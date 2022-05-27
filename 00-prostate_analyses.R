#############################
#     ANALYSES PROSTATE     #
#############################

library(dplyr)
library(stringr)
library(psy)
library(boot)
library(lme4)
library(pROC)
library(xlsx)


source("src/01-prostate_functions.R")


#=======================================
#DATA MANAGEMENT

roc <- read.csv2("data/sextant_20170407.csv") 
roc[roc$taille_IRM=="?", "taille_IRM"] <- NA
roc$taille_IRM <- as.numeric(as.character(roc$taille_IRM))
roc$sextants <- as.character(roc$sextants)
roc[is.na(roc$microcancer), "microcancer"] <- 0
roc$lobeD <- ifelse(roc$sextant %in% c("AD","MD","BD"),1,0)
#roc$taille_IRM <- ifelse (roc$microcancer == 1, 4, roc$taille_IRM) #NB : les microcancers sont-ils vus a l'IRM?

#imputation taille IRM 
roc %<>% left_join(roc %>% group_by(patient) %>% summarise(tmax = max(taille_IRM,na.rm=T))) %>% mutate(taille_IRM=tmax) %>% select(-tmax)

#tableau sextant : roc 

#tableau lobe
roc_lobe <- data.frame(roc %>% group_by(patient,lobeD) %>% select(ADK_histo, grep("_DCE",colnames(roc)), taille_IRM) %>%
                         summarise_each(funs(max(.,na.rm=T))))

#tableau patient
roc_pat <- data.frame(roc %>% group_by(patient) %>% select(ADK_histo, grep("_DCE",colnames(roc)), taille_IRM) %>%
                        summarise_each(funs(max(.,na.rm=T))))
#=============================

#tableau roc
#en histo, 3 patients n'ont pas de cancer, 34 patients ont un cancer, dont 21 patients sur un seul lobe et 13 patients sur les 2 lobes 
tmp <- roc %>% group_by(patient,lobeD) %>% summarise(K=max(ADK_histo)) %>% group_by(patient) %>% summarise(tot=sum(K))
table(tmp$tot)

#85 sextants avec un K en histo
table(roc$ADK_histo)

# moyennes et ecart-type sur l'ensemble des sextants, interet?
mean(roc$AL_DCE0)
mean(roc$AL_DCE1)
mean(roc$RP_DCE0)
mean(roc$RP_DCE1)



#fonctions utilisees pour Sensibilite_Specificite
# get_threshold : 
#   pour chaque variable dit si likert depasse le seuil de 3 et seuil de 4
# compute_se_sp :
#   pour chaque variable depassement du seuil 0/1, calcul Se et Sp par rapport a ADK_histo



#============================================
#toutes tailles
#1-choix de la taille de la tumeur et 0/1 selon seuil
roc1 <- roc
roc2 <- roc_lobe
roc3 <- roc_pat

#2-Calcul Se_sp et comparaison Se DCE0/DCE1
df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1:3){
  
  .l <- lapply(c(3:4), function(seuil){
    res <- compute_se_sp(.roc=get(df[i,"df"]), seuil, unit=df[i,"unit"], R=1000, type="bca") 
  })
  .l <- do.call(rbind, .l)
  
  sesp <- if (i==1) .l else rbind(sesp, .l)
}

sesp <- sesp %>% select(-N)


#3-Cohen's Kappa

df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1:3){
  .l <- lapply(c(3:4), function(seuil){
    res <- get_kappa_CI_conc (.roc=get(df[i,"df"]), seuil, unit=df[i,"unit"], R=1000, type="bca") 
  })
  .l <- do.call(rbind, .l)
  kappa_tot <- if (i==1) .l else rbind(kappa_tot, .l)
}


#4- AUC estimation
#var <- colnames(roc1[-c(1:4,9)])

df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)
graph_only = TRUE
graph_only = FALSE

for (i in 1:3) {
#for (i in 1:2) {
  data_tmp <- get(df[i,"df"])
  unit <- df[i,"unit"]
  
  data_tmp <- get_threshold(data_tmp)
  data_tmp$patient <- as.character(data_tmp$patient)
  
  var <- colnames(data_tmp)[grepl("DCE", colnames(data_tmp))]
  
  #Calcul de l'objet roc
  .l <- lapply(var, function(x) get_rocobj(data=data_tmp, var1=x, unit=unit, method = "glm"))

  impair.vec <- seq_along(.l)[seq_along(.l) %% 2 == 1]
  
  if (graph_only == FALSE){
    #intervalle de confiance
    dfAUC <- data.frame(variable=var,
                        do.call(rbind,lapply(1:12, function(x) round(as.numeric(ci.auc(.l[[x]], method="bootstrap", boot.n=1000, boot.stratified=FALSE)),3))))
    
    dfAUC$auc.95CI <- paste0(dfAUC$X2, "[", dfAUC$X1, "-", dfAUC$X3, "]")
    dfAUC <- dfAUC[,c("variable","auc.95CI")]
    
    #pvalue
    pvalue <- sapply(impair.vec, function(x) roc.test(.l[[x]],.l[[x+1]], method="bootstrap", boot.n=1000)$p.value) #seq_along(.l)[seq_along(.l) %% 2 == 1] : 1  3  5  7  9 11
    pvalue <- round(pvalue, 3)
    
    dfAUC$pvalue <- ""
    dfAUC$pvalue[impair.vec] <- pvalue
    
    dfAUC$unit <- unit
    dfAUC <- dfAUC[ ,c("unit", "variable", "auc.95CI", "pvalue")]
    
    if (i == 1 ) tab <- dfAUC
    else tab <- rbind(tab, dfAUC)
  }
  
  
  #browser()
  #le graphe sextant 1 est donc sextant toutes tailles, lobe1=lobe toutes tailles
  pdf(paste0("data/AUC_curve",unit,"all_size", Sys.Date(), ".pdf"))
  par(mfrow=c(2,2))
  for (j in impair.vec){
    seuil <- ifelse(nchar(var[j]) == 7, "none", str_sub(var[j],-1,-1))
    plot(.l[[j]], col="black",  lty= 2, ylim = c(0,1), xlim = c(1,0),main = paste0(unit, "; ", "juge = ", str_sub(var[j], 1,2), "; threshold = ", seuil))
    par(new=TRUE)
    plot(.l[[j+1]], lty=2, col = "grey",main = NULL)
    legend(0.5,0.22,  bty ="n",legend = c("injected", "not injected"),lty = c(2,2), col = c( "black","grey"))
  }
  
 dev.off()
}

tabAUC <- tab

write.xlsx(x = "toutes tailles", file = paste0("results/analyses_toutes_tailles_",Sys.Date(),".xlsx"), sheetName = "title", row.names = FALSE)
write.xlsx(x = sesp, file = paste0("results/analyses_toutes_tailles_",Sys.Date(),".xlsx"), sheetName = "Se Sp", row.names = FALSE, append = TRUE)
write.xlsx(x = kappa_tot, file = paste0("results/analyses_toutes_tailles_",Sys.Date(),".xlsx"), sheetName = "Kappa", row.names = FALSE, append=TRUE)
write.xlsx(x = tabAUC, file = paste0("results/analyses_toutes_tailles_",Sys.Date(),".xlsx"), sheetName = "AUC", row.names = FALSE, append=TRUE)
#============================================
#taille >=5

#1-choix de la taille de la tumeur et 0/1 selon seuil
roc1 <- roc[roc$taille_IRM>=5 & !is.na(roc$taille_IRM), ]
roc2 <- roc_lobe[roc_lobe$taille_IRM>=5 & !is.na(roc_lobe$taille_IRM), ]
roc3 <- roc_pat[roc_pat$taille_IRM>=5 & !is.na(roc_pat$taille_IRM), ]

#2-Calcul Se_sp et comparaison Se DCE0/DCE1
df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1:3){
  
  .l <- lapply(c(3:4), function(seuil){
    res <- compute_se_sp(.roc=get(df[i,"df"]), seuil, unit=df[i,"unit"], R=1000, type="bca") 
  })
  .l <- do.call(rbind, .l)
  
  sesp <- if (i==1) .l else rbind(sesp, .l)
}

sesp <- sesp %>% select(-N)


#3-Cohen's Kappa

df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1:3){
  .l <- lapply(c(3:4), function(seuil){
    res <- get_kappa_CI_conc (.roc=get(df[i,"df"]), seuil, unit=df[i,"unit"], R=1000, type="bca") 
  })
  .l <- do.call(rbind, .l)
  kappa_tot <- if (i==1) .l else rbind(kappa_tot, .l)
}



#4- AUC estimation
#var <- colnames(roc1[-c(1:4,9)])

get_rocobj <- function (data, var1, unit){
  #browser()
  tryCatch({
    print(var1)
    data$vartmp <- data[,var1]
    gm1 <- glmer(ADK_histo ~ vartmp + (1 | patient), data = data,
                 family = binomial)
    p <- as.numeric(predict(gm1, type="response"))
    if (var1 %in% c("AL_DCE1", "AL_DCE0", "RP_DCE1", "RP_DCE0") & unit!="patient" ) rocobj <- roc(data$ADK_histo, p, smooth=TRUE)
    else rocobj <- roc(data$ADK_histo, p, smooth=FALSE)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1) {
  data_tmp <- get(df[i,"df"])
  unit <- df[i,"unit"]
  
  data_tmp <- get_threshold(data_tmp)
  data_tmp$patient <- as.character(data_tmp$patient)
  
  var <- colnames(data_tmp)[grepl("DCE", colnames(data_tmp))]
  
  #Calcul de l'objet roc
  .l <- lapply(var, function(x) get_rocobj(data=data_tmp, var1=x, unit=unit))
  
  #intervalle de confiance
  dfAUC <- data.frame(variable=var, 
                      do.call(rbind,lapply(1:12, function(x) round(as.numeric(ci.auc(.l[[x]], method="bootstrap", boot.n=1000, boot.stratified=FALSE)),3))))
  dfAUC$auc.95CI <- paste0(dfAUC$X2, "[", dfAUC$X1, "-", dfAUC$X3, "]")
  dfAUC <- dfAUC[,c("variable","auc.95CI")]
  
  #pvalue
  impair.vec <- seq_along(.l)[seq_along(.l) %% 2 == 1]
  pvalue <- sapply(impair.vec, function(x) roc.test(.l[[x]],.l[[x+1]], method="bootstrap", boot.n=1000)$p.value) #seq_along(.l)[seq_along(.l) %% 2 == 1] : 1  3  5  7  9 11
  pvalue <- round(pvalue, 3)
  dfAUC$pvalue <- NA
  dfAUC$pvalue[impair.vec] <- pvalue
  
  dfAUC$unit <- unit
  dfAUC <- dfAUC[ ,c("unit", "variable", "auc.95CI", "pvalue")]
  
  if (i == 1 ) tab <- dfAUC
  else tab <- rbind(tab, dfAUC)
  
  pdf(paste0("data/AUC_curve",unit,"2.pdf"))
  par(mfrow=c(2,2))
  for (i in 1:length(var)){
    plot(.l[[i]], main = paste0(var[i],unit))
  }
  dev.off()
}

tabAUC <- tab
write.xlsx(x = "sup 5", file = paste0("results/analyses_sup5_",Sys.Date(),".xlsx"), sheetName = "title", row.names = FALSE)
write.xlsx(x = sesp, file = paste0("results/analyses_sup5_",Sys.Date(),".xlsx"), sheetName = "Se Sp", row.names = FALSE, append = TRUE)
write.xlsx(x = kappa_tot, file = paste0("results/analyses_sup5_",Sys.Date(),".xlsx"), sheetName = "Kappa", row.names = FALSE, append=TRUE)
write.xlsx(x = tabAUC, file = paste0("results/analyses_sup5_",Sys.Date(),".xlsx"), sheetName = "AUC", row.names = FALSE, append=TRUE)

#============================================
#taille <10

#1-choix de la taille de la tumeur et 0/1 selon seuil
roc1 <- roc[roc$taille_IRM<10 & !is.na(roc$taille_IRM), ]
roc2 <- roc_lobe[roc_lobe$taille_IRM<10 & !is.na(roc_lobe$taille_IRM), ]
roc3 <- roc_pat[roc_pat$taille_IRM<10 & !is.na(roc_pat$taille_IRM), ]

#2-Calcul Se_sp et comparaison Se DCE0/DCE1
df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1){
  
  .l <- lapply(c(3:4), function(seuil){
    res <- compute_se_sp(.roc=get(df[i,"df"]), seuil, unit=df[i,"unit"], R=1000, type="bca") 
  })
  .l <- do.call(rbind, .l)
  
  sesp <- if (i==1) .l else rbind(sesp, .l)
}

sesp <- sesp %>% select(-N)


#3-Cohen's Kappa

df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1:3){
  .l <- lapply(c(3:4), function(seuil){
    res <- get_kappa_CI_conc (.roc=get(df[i,"df"]), seuil, unit=df[i,"unit"], R=1000, type="bca") 
  })
  .l <- do.call(rbind, .l)
  kappa_tot <- if (i==1) .l else rbind(kappa_tot, .l)
}



#4- AUC estimation
#var <- colnames(roc1[-c(1:4,9)])




df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1) {
  data_tmp <- get(df[i,"df"])
  unit <- df[i,"unit"]
  
  data_tmp <- get_threshold(data_tmp)
  data_tmp$patient <- as.character(data_tmp$patient)
  
  var <- colnames(data_tmp)[grepl("DCE", colnames(data_tmp))]
  
  #Calcul de l'objet roc
  .l <- lapply(var, function(x) get_rocobj(data=data_tmp, var1=x, unit=unit))
  #intervalle de confiance
  dfAUC <- data.frame(variable=var, 
                      do.call(rbind,lapply(1:12, function(x) {
                        if(is.null(.l[[x]])) res <- c(NA, NA, NA)
                        else  res <- round(as.numeric(ci.auc(.l[[x]], method="bootstrap", boot.n=1000, boot.stratified=FALSE)),3)
                        return (res)
                        })))
  
  dfAUC$auc.95CI <- paste0(dfAUC$X2, "[", dfAUC$X1, "-", dfAUC$X3, "]")
  dfAUC <- dfAUC[,c("variable","auc.95CI")]
  
  #pvalue
  impair.vec <- seq_along(.l)[seq_along(.l) %% 2 == 1]
  pvalue <- sapply(impair.vec, function(x){
    if(is.null(.l[[x]]) | is.null(.l[[x+1]])) res <- NA
    else res <- roc.test(.l[[x]],.l[[x+1]], method="bootstrap", boot.n=1000)$p.value #seq_along(.l)[seq_along(.l) %% 2 == 1] : 1  3  5  7  9 11
    return(res)
  }) 
  pvalue <- round(pvalue, 3)
  dfAUC$pvalue <- NA
  dfAUC$pvalue[impair.vec] <- pvalue
  
  dfAUC$unit <- unit
  dfAUC <- dfAUC[ ,c("unit", "variable", "auc.95CI", "pvalue")]
  
  if (i == 1 ) tab <- dfAUC
  else tab <- rbind(tab, dfAUC)
  
  pdf(paste0("data/AUC_curve",unit,"3.pdf"))
  par(mfrow=c(2,2))
  for (i in 1:length(var)){
    plot(.l[[i]], main = paste0(var[i],unit))
  }
  dev.off()
}

tabAUC <- tab

write.xlsx(x = "inf 10", file = paste0("results/analyses_inf10_",Sys.Date(),".xlsx"), sheetName = "Title", row.names = FALSE)
write.xlsx(x = sesp, file = paste0("results/analyses_inf10_",Sys.Date(),".xlsx"), sheetName = "Se Sp", row.names = FALSE, append = TRUE)
write.xlsx(x = kappa_tot, file = paste0("results/analyses_inf10_",Sys.Date(),".xlsx"), sheetName = "Kappa", row.names = FALSE, append=TRUE)
write.xlsx(x = tabAUC, file = paste0("results/analyses_inf10_",Sys.Date(),".xlsx"), sheetName = "AUC", row.names = FALSE, append=TRUE)

#============================================
#taille >=10

#1-choix de la taille de la tumeur et 0/1 selon seuil
roc1 <- roc[roc$taille_IRM>=10 & !is.na(roc$taille_IRM), ]
roc2 <- roc_lobe[roc_lobe$taille_IRM>=10 & !is.na(roc_lobe$taille_IRM), ]
roc3 <- roc_pat[roc_pat$taille_IRM>=10 & !is.na(roc_pat$taille_IRM), ]

#2-Calcul Se_sp et comparaison Se DCE0/DCE1
df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1:3){
  
  .l <- lapply(c(3:4), function(seuil){
    res <- compute_se_sp(.roc=get(df[i,"df"]), seuil, unit=df[i,"unit"], R=1000, type="bca") 
  })
  .l <- do.call(rbind, .l)
  
  sesp <- if (i==1) .l else rbind(sesp, .l)
}

sesp <- sesp %>% select(-N)


#3-Cohen's Kappa

df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1:3){
  .l <- lapply(c(3:4), function(seuil){
    res <- get_kappa_CI_conc (.roc=get(df[i,"df"]), seuil, unit=df[i,"unit"], R=1000, type="bca") 
  })
  .l <- do.call(rbind, .l)
  kappa_tot <- if (i==1) .l else rbind(kappa_tot, .l)
}

kappa_tot <- kappa_tot

#4- AUC estimation
#var <- colnames(roc1[-c(1:4,9)])

get_rocobj <- function (data, var1, unit){
  #browser()
  tryCatch({
    print(var1)
    data$vartmp <- data[,var1]
    gm1 <- glmer(ADK_histo ~ vartmp + (1 | patient), data = data,
                 family = binomial)
    # gm1 <- glmer(data$ADK_histo ~ data[ ,var1] + (1 | data$patient),
    #              family = binomial)
    p <- as.numeric(predict(gm1, type="response"))
    if (var1 %in% c("AL_DCE1", "AL_DCE0", "RP_DCE1", "RP_DCE0") & unit!="patient" ) rocobj <- roc(data$ADK_histo, p, smooth=TRUE)
    else rocobj <- roc(data$ADK_histo, p, smooth=FALSE)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


df <- data.frame(df = c("roc1", "roc2", "roc3"),unit=c("sextant","lobe","patient"), stringsAsFactors=FALSE)

for (i in 1) {
  data_tmp <- get(df[i,"df"])
  unit <- df[i,"unit"]
  
  data_tmp <- get_threshold(data_tmp)
  data_tmp$patient <- as.character(data_tmp$patient)
  
  var <- colnames(data_tmp)[grepl("DCE", colnames(data_tmp))]
  
  #Calcul de l'objet roc
  .l <- lapply(var, function(x) get_rocobj(data=data_tmp, var1=x, unit=unit))
  
  #intervalle de confiance
  dfAUC <- data.frame(variable=var, 
                      do.call(rbind,lapply(1:12, function(x) round(as.numeric(ci.auc(.l[[x]], method="bootstrap", boot.n=1000, boot.stratified=FALSE)),3))))
  
  dfAUC$auc.95CI <- paste0(dfAUC$X2, "[", dfAUC$X1, "-", dfAUC$X3, "]")
  dfAUC <- dfAUC[,c("variable","auc.95CI")]
  
  #pvalue
  impair.vec <- seq_along(.l)[seq_along(.l) %% 2 == 1]
  pvalue <- sapply(impair.vec, function(x) roc.test(.l[[x]],.l[[x+1]], method="bootstrap", boot.n=1000)$p.value) #seq_along(.l)[seq_along(.l) %% 2 == 1] : 1  3  5  7  9 11
  pvalue <- round(pvalue, 3)
  dfAUC$pvalue <- NA
  dfAUC$pvalue[impair.vec] <- pvalue
  
  dfAUC$unit <- unit
  dfAUC <- dfAUC[ ,c("unit", "variable", "auc.95CI", "pvalue")]
  
  if (i == 1 ) tab <- dfAUC
  else tab <- rbind(tab, dfAUC)
  
  pdf(paste0("data/AUC_curve",unit,"4.pdf"))
  par(mfrow=c(2,2))
  for (i in 1:length(var)){
    plot(.l[[i]], main = paste0(var[i],unit))
  }
  dev.off()
}

tabAUC <- tab


write.xlsx(x = "sup 10", file = paste0("results/analyses_sup10_",Sys.Date(),".xlsx"), sheetName = "Title", row.names = FALSE)
write.xlsx(x = sesp, file = paste0("results/analyses_sup10_",Sys.Date(),".xlsx"), sheetName = "Se Sp", row.names = FALSE, append = TRUE)
write.xlsx(x = kappa_tot, file = paste0("results/analyses_sup10_",Sys.Date(),".xlsx"), sheetName = "Kappa", row.names = FALSE, append=TRUE)
write.xlsx(x = tabAUC, file = paste0("results/analyses_sup10_",Sys.Date(),".xlsx"), sheetName = "AUC", row.names = FALSE, append=TRUE)


#----------------------
#Analyse complementaire

d <- read.csv2("data/BPC_vs_BPST.csv")

d$BPST_lgmax <- as.numeric(as.character(d$BPST_lgmax))
d$BPC_lgmax <- as.numeric(as.character(d$BPC_lgmax))

table(d$BPST_positif)
table(d$BPC_positif)
table(d$BPST_positif, d$BPC_positif)

round(prop.table(table(d$BPST_positif)),2)
round(prop.table(table(d$BPC_positif)),2)
fisher.test(d$BPST_positif, d$BPC_positif)
mcnemar.test(d$BPST_positif, d$BPC_positif)

wilcox.test(d$BPST_lgmax, d$BPC_lgmax)

xBPST<-d$BPST_lgmax
xBPC<-d$BPC_lgmax

keep <- c(xBPST,xBPC)
obs0<- mean(xBPST)
obs1<- mean(xBPC)
#diff.obs <- abs(obs0 - obs1)
diff.obs <- obs0 - obs1

.gpe <- rep(c("gBPST","gBPC"),c(length(xBPST),length(xBPC)))

perm.test<- function(keepIT=keep,.gpe=.gpe){
  mixgpe <- sample(.gpe,replace = FALSE)
  g0 <- keepIT[mixgpe=="gBPST"]
  g1 <- keepIT[mixgpe=="gBPC"]
  m0<-mean(g0)
  m1<-mean(g1)
  #diff<-abs(m0-m1)
  diff<-m0-m1
  return(diff)
}

many.samp <- replicate (100000, perm.test(keep,.gpe))

p.val <-length(many.samp[abs(many.samp)>= abs(diff.obs)]) / length(many.samp)

# #pour tracer courbe
hist(many.samp,main=paste0("Difference de moyenne ",x))
abline(v=diff.obs,lwd=2,col=2)

cat(paste0("moyenne longueur carotte BPST = ", round(obs0, 2), "\nmoyenne longueur carotte BPC = ", round(obs1, 2), "\npval test de permutation = ", p.val ))
cat(sd(xBPST), sd(xBPC))  
