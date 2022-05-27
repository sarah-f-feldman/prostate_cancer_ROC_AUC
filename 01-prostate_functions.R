###########################
#   Prostate functions    #
###########################



get_threshold <- function(.roc){
  roc2 <- .roc 
  #ajout variable 0/1 pour seuil 3
  for (i in colnames(roc2)[str_sub(colnames(roc2),-5,-5)=="_"]){
    roc2[ ,paste0(i,"_3")] <- ifelse(roc2[ ,i] >= 3, 1, 0)
  }
  #ajout variable 0/1 pour seuil 4
  for (i in colnames(roc2)[str_sub(colnames(roc2),-5,-5)=="_"]){
    roc2[ ,paste0(i,"_4")] <- ifelse(roc2[ ,i] >= 4, 1, 0)
  }
  return(roc2)
}

compute_se_sp <- function(.roc, seuil, unit, R=1000, type="bca") {
  roc2 <- data.frame(get_threshold(.roc))
  
  vary <- data.frame(DCE1=paste0(c("AL","RP","AL","RP"),rep(c("_DCE1_3", "_DCE1_4"), each=2)), DCE0=paste0(c("AL","RP","AL","RP"),rep(c("_DCE0_3", "_DCE0_4"), each=2)))
  #vary <- apply(vary,2, as.character)
  #roc2 <- data.frame(.roc)
  
  #selection des colonnes "dépassement du seuil oui ou non"
  col_to_test <- colnames(roc2)[grep(seuil, colnames(roc2))] #"RP_DCE0_4" a tiret en 2e position en partant de la fin
  #pour chacune de ces colonnes, calcul se sp en comparant avec ADK_histo
  se_sp <- lapply(col_to_test, function(x){
    print(x)
    tmp <- roc2[,c("ADK_histo",x)]
    tab <- table(tmp[ ,x],tmp$ADK_histo,deparse.level = 2)
    #print(tab)
    #se : proba (signe/maladie)= proba (AL_DCE1_3=1 sachant ADK=1)
    #se <- paste0(round(tab[2,2]/sum(tab[,2]),2)*100,"%")
    #browser()
    se <- paste0(round(tab[rownames(tab)==1,colnames(tab)[colnames(tab)==1]]/sum(tab[,colnames(tab)[colnames(tab)==1]]),2)*100,"%")
    #sp : proba(pas de signe sachant pas malade = proba (AL_DCE1_3=0 sachant ADK=0))
    
    se_CI <- BootseCi(tmp, "se", R, type)
    #if(x=="RP_DCE0_3")browser()

    if(any(rownames(tab)==0) == FALSE) sp_CI <- NA 
    #if(any(colnames(tab)==0) == FALSE) sp <- "No M- with histo" 
    #if(any(colnames(tab)==0) == FALSE) browser()
    if((any(colnames(tab)==0)== FALSE & any(rownames(tab)==0)== FALSE) == TRUE) sp_CI <- NA 
    else {
      #sp <- paste0(round(tab[1,1]/sum(tab[,1]),2)*100,"%")
      sp <- paste0(round(tab[rownames(tab)==0,colnames(tab)[colnames(tab)==0]]/sum(tab[,colnames(tab)[colnames(tab)==0]]),2)*100,"%")
      sp_CI <- BootseCi(tmp, mes="sp" , R, type)
      #browser()
    }
    N <- nrow(tmp)
    #return(c(se,sp))
    #browser()
    nM0S0 <- tab[rownames(tab)==0, colnames(tab)==0]
    nM0S1 <- tab[rownames(tab)==1, colnames(tab)==0]
    nM1S0 <- tab[rownames(tab)==0, colnames(tab)==1]
    nM1S1 <- tab[rownames(tab)==1, colnames(tab)==1]
    for (i in c("nM0S0",  "nM0S1", "nM1S0", "nM1S1")){
    #lapply(c(nM0S0,  nM0S1, nM1S0, nM1S1), function(i){
      obj <- get(i)
      if (length(obj)==0) obj <- 0
      assign(i, obj)
    }
    #browser()
    if (str_sub(x, 7, 7)=="0") {
      #fisher.se <- ""
      mcnemar.se <- ""
      #fisher.sp <- ""
      mcnemar.sp <- ""
    } else {
      tab.se <- data.frame(t(test_SE_DCE(juge_DCE1_seuil=x, juge_DCE0_seuil=as.character(vary[vary$DCE1==x, 2]), .roc)), stringsAsFactors = FALSE)
      #fisher.se <- tab.se$fisher
      mcnemar.se <- tab.se$mcnemar
      tab.sp <- data.frame(t(test_SP_DCE(juge_DCE1_seuil=x, juge_DCE0_seuil=as.character(vary[vary$DCE1==x, 2]), .roc)), stringsAsFactors = FALSE)
      #fisher.sp <- tab.sp$fisher
      mcnemar.sp <- tab.sp$mcnemar
      #browser()
    }
    
    #if(x=="AL_DCE0_3" & unit =="patient")browser()
    #.df <- cbind(N, se_CI, sp_CI, nM0S0, nM0S1, nM1S0, nM1S1) ; colnames(.df) <- c("N","se_CI","sp_CI", "M-S-", "M-S+", "M+S-", "M+S+") #obligatoire pour que les NA ait un colnames aussi, sinon rbind impossible
    #.df <- cbind(N, se_CI, fisher, mcnemar, sp_CI, nM0S0, nM0S1, nM1S0, nM1S1) ; colnames(.df) <- c("N","se_CI","fisher","mcnemar", "sp_CI", "nM0S0", "nM0S1", "nM1S0", "nM1S1") #obligatoire pour que les NA ait un colnames aussi, sinon rbind impossible
    #if (str_sub(x, 7, 7)=="0")browser()
    #.df <- cbind(nM0S0, nM0S1, nM1S0, nM1S1, N, se_CI, sp_CI, fisher.se, mcnemar.se, fisher.sp, mcnemar.sp) ; colnames(.df) <- c("VN", "FP", "FN", "VP", "N", "se_CI", "sp_CI", "fisher.se", "mcnemar.se", "fisher.sp", "mcnemar.sp") #obligatoire pour que les NA ait un colnames aussi, sinon rbind impossible
    .df <- cbind(nM0S0, nM0S1, nM1S0, nM1S1, N, se_CI, sp_CI, mcnemar.se, mcnemar.sp) ; colnames(.df) <- c("VN", "FP", "FN", "VP", "N", "se_CI", "sp_CI", "mcnemar.se", "mcnemar.sp") #obligatoire pour que les NA ait un colnames aussi, sinon rbind impossible
    .df <- data.frame(.df,stringsAsFactors = FALSE)
    return(.df)
  })
  #if (unit=="patient")browser()
  se_sp <- do.call(rbind, se_sp)
  #chaque liste/ligne est une colonne "dépassement du seuil oui ou non" qui nous donne dans le nom de la colonne DCE, seuil et juge.
  #sesp <- data.frame(mesure = col_to_test, N = se_sp[ , 1], Se = se_sp[ , 2], Sp = se_sp[ , 3])
  sesp <- data.frame(mesure = col_to_test, se_sp)
  sesp$unit <- unit
  sesp$DCE <- str_sub(sesp$mesure, 7, 7)
  sesp$seuil <- str_sub(sesp$mesure, -1, -1)
  sesp$juge <- str_sub(sesp$mesure, 1, 2)
  
  #browser()
  #sesp <- sesp[ ,c(4,6,7,5,2:3)]
  sesp <- sesp[ ,c(11:14,2:10)]
  rownames(sesp) <- NULL
  return(sesp)
}

# rocbis <- get_threshold(roc_lobe)
# sesp_lobe <- compute_se_sp (data.frame(rocbis), seuil=3, unit="lobe")
# sesp_lobe <- lapply(3:4, function(x) compute_se_sp (roc_lobe, seuil=x, unit="lobe"))
# sesp_lobe <- do.call(rbind,sesp_lobe)

#sextant
# #seuil
# roc2 <- roc 
# #ajout variable 0/1 pour seuil 3
# for (i in colnames(roc2)[5:8]){
#   roc2[ ,paste0(i,"_3")] <- ifelse(roc2[ ,i] >= 3, 1, 0)
# }
# #ajout variable 0/1 pour seuil 4
# for (i in colnames(roc2)[5:8]){
#   roc2[ ,paste0(i,"_4")] <- ifelse(roc2[ ,i] >= 4, 1, 0)
# }
# 
# #calcul se sp
# #selection des colonnes "dépassement du seuil oui ou non"
# col_to_test <- colnames(roc2)[c(grep("DCE0_",colnames(roc2)),grep("DCE1_",colnames(roc2)))]
# #pour chacune de ces colonnes, calcul se sp en comparant avec ADK_histo
# se_sp <- lapply(col_to_test, function(x){
#          tmp <- roc2[,c("ADK_histo",x)]
#          tab <- table(tmp[ ,x],tmp$ADK_histo,deparse.level = 2)
#          #se : proba (signe/maladie)= proba (AL_DCE1_3=1 sachant ADK=1)
#          se <- paste0(round(tab[2,2]/sum(tab[,2]),2)*100,"%")
#          #sp : proba(pas de signe sachant pas malade = proba (AL_DCE1_3=0 sachant ADK=0))
#          sp <- paste0(round(tab[1,1]/sum(tab[,1]),2)*100,"%")
#          return(c(se,sp))
#        })
# se_sp <- do.call(rbind, se_sp)
# #chaque liste/ligne est une colonne "dépassement du seuil oui ou non" qui nous donne dans le nom de la colonne DCE, seuil et juge.
# sesp_sextant <- data.frame(mesure = col_to_test, Se = se_sp[,1], Sp = se_sp[,2])
# sesp_sextant$DCE <- str_sub(sesp_sextant$mesure, 7, 7)
# sesp_sextant$seuil <- str_sub(sesp_sextant$mesure, -1, -1)
# sesp_sextant$juge <- str_sub(sesp_sextant$mesure, 1, 2)
# 
# sesp_sextant <- sesp_sextant[ ,c(1,4:6,2:3)]

# #lobe
# roc2 <- roc_lobe 
# #ajout variable 0/1 pour seuil 3
# for (i in colnames(roc2)[str_sub(colnames(roc2),-5,-5)=="_"]){
#   roc2[ ,paste0(i,"_3")] <- ifelse(roc2[ ,i] >= 3, 1, 0)
# }
# #ajout variable 0/1 pour seuil 4
# for (i in colnames(roc2)[str_sub(colnames(roc2),-5,-5)=="_"]){
#   roc2[ ,paste0(i,"_4")] <- ifelse(roc2[ ,i] >= 4, 1, 0)
# }
# 
# #calcul se sp
# #selection des colonnes "dépassement du seuil oui ou non"
# col_to_test <- colnames(roc2)[str_sub(colnames(roc2),-2,-2)=="_"] #"RP_DCE0_4" a tiret en 2e position en partant de la fin
# #pour chacune de ces colonnes, calcul se sp en comparant avec ADK_histo
# se_sp <- lapply(col_to_test, function(x){
#   tmp <- roc2[,c("ADK_histo",x)]
#   tab <- table(tmp[ ,x],tmp$ADK_histo,deparse.level = 2)
#   #se : proba (signe/maladie)= proba (AL_DCE1_3=1 sachant ADK=1)
#   se <- paste0(round(tab[2,2]/sum(tab[,2]),2)*100,"%")
#   #sp : proba(pas de signe sachant pas malade = proba (AL_DCE1_3=0 sachant ADK=0))
#   sp <- paste0(round(tab[1,1]/sum(tab[,1]),2)*100,"%")
#   return(c(se,sp))
# })
# se_sp <- do.call(rbind, se_sp)
# #chaque liste/ligne est une colonne "dépassement du seuil oui ou non" qui nous donne dans le nom de la colonne DCE, seuil et juge.
# sesp_lobe <- data.frame(mesure = col_to_test, Se = se_sp[,1], Sp = se_sp[,2])
# sesp_lobe$DCE <- str_sub(sesp_lobe$mesure, 7, 7)
# sesp_lobe$seuil <- str_sub(sesp_lobe$mesure, -1, -1)
# sesp_lobe$juge <- str_sub(sesp_lobe$mesure, 1, 2)
# 
# sesp_lobe <- sesp_lobe[ ,c(1,4:6,2:3)]



se.boot<- function(data, .mes, indices){ 
  .dat<- data[indices,]
  tab <- table(.dat[ ,2],.dat$ADK_histo,deparse.level = 2)
  if(.mes=="se") {
    if (any(colnames(tab)==1) == FALSE | any(rownames(tab)==1) == FALSE) mesure <- NA
    else mesure <- tab[rownames(tab)==1,colnames(tab)[colnames(tab)==1]]/sum(tab[,colnames(tab)[colnames(tab)==1]])
  }
  if(.mes=="sp") {
    if (any(colnames(tab)==0) == FALSE | any(rownames(tab)==0) == FALSE) mesure <- NA
    else mesure <- tab[rownames(tab)==0,colnames(tab)[colnames(tab)==0]]/sum(tab[,colnames(tab)[colnames(tab)==0]])
  }
  return(mesure)
}

bootse<- function (tmp, mes, R)  {    #x = "1":6 ou x="tot"
  .df <- tmp
  .res<- boot(data=.df, statistic = se.boot , .mes=mes, R=R)
  return(.res)
}



BootseCi <- function(tmp, mes, R, type)  { #type="bca" ou "perc"
  .bootres <- bootse (tmp=tmp, mes, R=R)
  if (all(na.omit(.bootres$t)==na.omit(.bootres$t)[1])) return(paste0(as.numeric (.bootres$t0), "[ND]"))
  #si bug, décommenter ligne ci dessous
  #if (all(na.omit(.bootres$t)==na.omit(.bootres$t)[1])) return(NA)
  
  
  .n <- length (.bootres$t0) #donne le nombre de resultat boot realise : 1 pour internet, 1 pour telephone
  .list.ci <- lapply(1:.n, function(x) boot.ci(.bootres,index=x,type=type)) #fct boot.ci : intervalle de confiance pour chaque boot
  if (type=="perc") type2 <- "percent" 
  else type2 <- type
  .res <- data.frame (t (sapply (.list.ci, function (x) x[[type2]][4:5]))) #selectionne les valeur de IC
  rownames (.res) <- names (.bootres$t0)
  #print(.res)
  #print(str(.res))
  colnames (.res) <- c ("CI_L", "CI_U")
  .res$est <- as.numeric (.bootres$t0)
  .res$n <- nrow(tmp)
  .res <- .res[, c (4,3, 1, 2)]
  .ans <- round (.res, 2) #fait un arrondi sur chaque valeur
  #if(mes=="se").ans <- data.frame (N=.res$n, Se_CI=paste0 (.ans$est*100, "% [", .ans$CI_L*100, "%-", .ans$CI_U*100, "%]")) #met en forme les valeurs
  if(mes=="se").ans <- data.frame (Se_CI=paste0 (.ans$est*100, "% [", .ans$CI_L*100, "%-", .ans$CI_U*100, "%]")) #met en forme les valeurs
  if(mes=="sp").ans <- data.frame (Sp_CI=paste0 (.ans$est*100, "% [", .ans$CI_L*100, "%-", .ans$CI_U*100, "%]")) #met en forme les valeurs
  .ans <- as.vector(.ans)
  return (.ans)
}

#comparaison de sensibilité, DCE vs pas de DCE
test_SE_DCE <- function(juge_DCE1_seuil, juge_DCE0_seuil, data){
  rocbis <- data
  rocbis <- get_threshold(rocbis)
  comp <- rocbis[rocbis$ADK_histo==1, c(juge_DCE1_seuil, juge_DCE0_seuil)] 
  #FAUX : ne donne pas concordance entre DCE0 et DCE1
  # comp <- data.frame(signe=c(comp[,1], comp[,2]))
  # comp$DCE <- c(rep(1, nrow(rocbis[ rocbis$ADK_histo==1, ])), rep(0, nrow(rocbis[ rocbis$ADK_histo==1, ])))
  if (sum(dim(table(comp))) < 4) {
    #print("pb")
    #pvalft <- NA
    pvalmt <- NA
  } else {
    # ft <- fisher.test(comp$signe,comp$DCE)
    # pvalft <- round(ft$p.value, 3)
    mt <- mcnemar.test(table(comp))
    pvalmt <- round(mt$p.value,3)
  }
  #df <- c(juge = str_sub(juge_DCE1_seuil, 1, 2), seuil = str_sub(juge_DCE1_seuil, -1, -1), fisher = pvalft, mcnemar = pvalmt)
  df <- c(juge = str_sub(juge_DCE1_seuil, 1, 2), seuil = str_sub(juge_DCE1_seuil, -1, -1), mcnemar = pvalmt)
  return(df)
}

#comparaison de spécificité, DCE vs pas de DCE
test_SP_DCE <- function(juge_DCE1_seuil, juge_DCE0_seuil, data){
  rocbis <- data
  rocbis <- get_threshold(rocbis)
  comp <- rocbis[rocbis$ADK_histo==0, c(juge_DCE1_seuil, juge_DCE0_seuil)] 
  #browser()
  #ARCHIFAUX pour le mcnemar, ne donne pas la concordance entre DCE0 et 1, juste une répartition des VN et FP  
  # comp <- data.frame(signe=c(comp[,1], comp[,2]))
  # comp$DCE <- c(rep(1, nrow(rocbis[ rocbis$ADK_histo==0, ])), rep(0, nrow(rocbis[ rocbis$ADK_histo==0, ])))
  
  if (sum(dim(table(comp))) < 4) {
    #pvalft <- NA
    pvalmt <- NA
  } else {
    # ft <- fisher.test(comp$signe,comp$DCE)
    # pvalft <- round(ft$p.value, 3)
    mt <- mcnemar.test(table(comp))
    pvalmt <- round(mt$p.value,3)
  }
  #df <- c(juge = str_sub(juge_DCE1_seuil, 1, 2), seuil = str_sub(juge_DCE1_seuil, -1, -1), fisher = pvalft, mcnemar = pvalmt)
  df <- c(juge = str_sub(juge_DCE1_seuil, 1, 2), seuil = str_sub(juge_DCE1_seuil, -1, -1), mcnemar = pvalmt)
  return(df)
}


#KAPPA

get_kappa_CI_conc <- function (.roc, seuil, unit, R=1000, type="bca"){
  roc2 <- data.frame(get_threshold(.roc))
  for (i in c(0,1)){ #pour DCE0 puis pour DCE1
    #browser()
    AL <- paste0("AL_DCE",i,"_",seuil)
    RP <- paste0("RP_DCE",i,"_",seuil)
    print(paste(seuil,unit,i))
    #if(seuil==3 & unit=="lobe" & i==1) browser()
    tmp <- roc2[ ,c(AL, RP)]
    kappa <- ckappa.cmarge(tmp) 
    K_CI <- if (is.na(kappa$kappa))  paste0("NA [NA]") else bootkapa.ci (tmp, R, type="bca")
    conc_marge <- kappa$conc.marge
    df <- data.frame (juge1 = AL, juge2 = RP, kappa = K_CI, concordance_marge = conc_marge)
    colnames(df) <- c("juge1","juge2","Kappa[95%CI]","concordance des marges")
    res <- if (i==0) df else rbind(res, df)
  }
  res$unit <- unit
  res$seuil <- seuil
  res$DCE <- str_sub(res$juge1, -3,-3)
  res <- res[ ,c(5:7,3,4)]
  #browser()
  return(res)
}

ckappa.cmarge <- function (r) {
  r <- na.omit(r)
  r1 <- r[, 1]
  r2 <- r[, 2]
  n1 <- as.character(r1)
  n2 <- as.character(r2)
  lev <- levels(as.factor(c(n1, n2)))
  p <- length(lev)
  tab <- matrix(nrow = p, ncol = p)
  dimnames(tab) <- list(levels(as.factor(c(n1, n2))), levels(as.factor(c(n1, 
                                                                         n2))))
  dim1 <- dimnames(tab)[[1]]
  dim2 <- dimnames(tab)[[2]]
  tabi <- table(n1, n2)
  dimi1 <- dimnames(tabi)[[1]]
  dimi2 <- dimnames(tabi)[[2]]
  for (i in 1:p) for (j in 1:p) {
    if ((sum(dim1[i] == dimi1) == 1) & (sum(dim2[j] == dimi2) == 
                                        1)) 
      tab[i, j] <- tabi[dim1[i], dim2[j]]
    else tab[i, j] <- 0
  }
  tsum <- sum(tab)
  ttab <- tab/tsum
  tm1 <- apply(ttab, 1, sum)
  tm2 <- apply(ttab, 2, sum)
  agreeP <- sum(diag(ttab))
  chanceP <- sum(tm1 * tm2)
  kappa2 <- (agreeP - chanceP)/(1 - chanceP)
  
  PM0n1 <- tm2[1]
  PM1n1 <- tm2[2]
  PM0n2 <- tm1[1]
  PM1n2 <- tm1[2]
  cond <- (as.numeric(PM0n1-PM1n1)>0) == (as.numeric(PM0n2-PM1n2)>0)
  
  
  result <- list(table = tab, kappa = kappa2, conc.marge = cond)
  result
}

#ckappa.boot <- function(data,x) {ckappa(data[x,])[[2]]}
ckappa.boot <- function(data, indices) {ckappa(data[indices, ])[[2]]}

bootkapa.ci <- function(tmp, R=1000, type) {
  .bootres <- boot(data=tmp, ckappa.boot, R) #function(data,x) {ckappa(data[x,])[[2]]}
  .n <- length (.bootres$t0)
  if (all(na.omit(.bootres$t)==na.omit(.bootres$t)[1])) return(paste0(as.numeric (.bootres$t0), "[ND]"))
  .list.ci <- lapply(1:.n, function(x) boot.ci(.bootres,index=x,type=type))
  if (type=="perc") type2 <- "percent" 
  else type2 <- type
  .res <- data.frame (t (sapply (.list.ci, function (x) x[[type2]][4:5]))) #selectionne les valeur de IC
  rownames (.res) <- names (.bootres$t0)
  colnames (.res) <- c ("CI_L", "CI_U")
  .res$est <- as.numeric (.bootres$t0)
  .res$n <- nrow(tmp)
  .res <- .res[, c (4,3, 1, 2)]
  .ans <- round (.res, 2) #fait un arrondi sur chaque valeur
  .ans <- data.frame (Kappa_CI=paste0 (.ans$est, " [", .ans$CI_L, " - ", .ans$CI_U, "]"))
  .ans <- as.vector(.ans)
  return(.ans)
}

# get_rocobj <- function (data, var1, unit){
#   #browser()
#   #tryCatch({
#   #if(var1=="AL_DCE1" & unit=="lobe") browser() 
#     print(paste(var1, unit))
#     data$vartmp <- data[,var1]
#     gm1 <- glmer(ADK_histo ~ vartmp + (1 | patient), data = data,
#                  family = binomial)
#     # gm1 <- glmer(data$ADK_histo ~ data[ ,var1] + (1 | data$patient),
#     #              family = binomial)
#     p <- as.numeric(predict(gm1, type="response"))
#     if (var1 %in% c("AL_DCE1", "AL_DCE0", "RP_DCE1", "RP_DCE0") & unit!="patient" ) rocobj <- roc(data$ADK_histo, p, smooth=TRUE)
#     else rocobj <- roc(data$ADK_histo, p, smooth=FALSE)
#   #}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
get_rocobj <- function (data, var1, unit, method = "glmer"){
  #browser()
  tryCatch({
    print(var1)
    data$vartmp <- data[,var1]
    if (method == "glmer") gm1 <- glmer(ADK_histo ~ vartmp + (1 | patient), data = data, family = binomial)
    if (method == "glm") gm1 <- glm(ADK_histo ~ vartmp, data = data, family = binomial)
    # gm1 <- glmer(data$ADK_histo ~ data[ ,var1] + (1 | data$patient),
    #              family = binomial)
    p <- as.numeric(predict(gm1, type="response"))
    if (var1 %in% c("AL_DCE1", "AL_DCE0", "RP_DCE1", "RP_DCE0") & unit!="patient" ) rocobj <- roc(data$ADK_histo, p, smooth=FALSE)
    else rocobj <- roc(data$ADK_histo, p, smooth=FALSE)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
