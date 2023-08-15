#### Phenomic and Genomic Prediction of Maize Agronomic Traits Across and Within Environments ####
# Authors: Aaron J. DeSalvio and Alper Adak #

#### BGLR Prediction Script 1 ####
# This script is intended for use with traits that possess NIR data for all four kernels (CS11_WS, CS11_WW, CS12_WS, and CS12_WW). #
# These traits are: ASI, DTA, DTS, EH, KW, PH, and Yield #

setwd("Set your directory here")

library(dplyr)
library(plyr)
library(AGHmatrix)
library(gplots)
library(dplyr)
library(reshape)
library(reshape2)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(data.table)
library(BGLR)
library(caret)
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

traitnames <- c("ASI", "DTA", "DTS", "EH", "KW", "PH", "Yield")
models <- c("Eta1", "Eta2", "Eta3")

for (tr in 1:length(traitnames)) {
  
  if (tr == 1) {
    Pheno <- read.csv(file = paste("Pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr == 2) {
    Pheno <- read.csv(file = paste("Pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr == 3) {
    Pheno <- read.csv(file = paste("Pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr == 4) {
    Pheno <- read.csv(file = paste("Pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr == 5) {
    Pheno <- read.csv(file = paste("Pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr == 6) {
    Pheno <- read.csv(file = paste("Pheno_", traitnames[tr], ".csv", sep = ""))
  } else {
    Pheno <- read.csv(file = paste("Pheno_", traitnames[tr], ".csv", sep = ""))
  }
  
  
  Pheno$Env<-as.factor(Pheno$Env)
  ZE<-model.matrix(~Pheno$Env-1)
  ZE[1:10,1:4]
  dim(ZE)
  Pheno$Hybrid<-as.factor(Pheno$Hybrid)
  ZVar<-model.matrix(~Pheno$Hybrid-1)
  dim(ZVar)
  
  
  #### Genomic effect kernel ####
  G <-as.matrix(read.csv("G.csv",row.names = 1))
  ZG<-ZVar%*%G%*%t(ZVar)
  ZG[1:5,1:5]
  
  
  #### NIR kernel with centering and scaling ####
  NIR_CS11_WS <-as.matrix(read.csv("BLUP_CS11_WS.csv", row.names = 1))
  NIR_CS11_WS_sc<-scale(NIR_CS11_WS,center=TRUE,scale=TRUE)
  
  NIR_CS11_WW <-as.matrix(read.csv("BLUP_CS11_WW.csv", row.names = 1))
  NIR_CS11_WW_sc<-scale(NIR_CS11_WW,center=TRUE,scale=TRUE)
  
  NIR_CS12_WS <-as.matrix(read.csv("BLUP_CS12_WS.csv", row.names = 1))
  NIR_CS12_WS_sc<-scale(NIR_CS12_WS,center=TRUE,scale=TRUE)
  
  NIR_CS12_WW <-as.matrix(read.csv("BLUP_CS12_WW.csv", row.names = 1))
  NIR_CS12_WW_sc<-scale(NIR_CS12_WW,center=TRUE,scale=TRUE)
  
  v <-  as.matrix(rbind(NIR_CS11_WS_sc,NIR_CS11_WW_sc,NIR_CS12_WS_sc,NIR_CS12_WW_sc))
  ZN <-tcrossprod(v)/ncol(v)
  
  
  #### Interaction kernels ####
  ZEZE<-tcrossprod(ZE)
  ZGZE<-ZG*ZEZE # Interaction kernel between genomic kernel and environment kernel
  ZNZE<-ZN*ZEZE # Interaction kernel between phenomic kernel and environment kernel
  
  
  Eta1<-list(E=list(X=ZE,model="BRR"),
             G=list(K=ZG,model="RKHS"),
             EG=list(K=ZGZE,model="RKHS")) 
  # y = ZE + ZG + error
  
  Eta2<-list(E=list(X=ZE,model="BRR"),
             N=list(K=ZN,model="RKHS"),
             EN=list(K=ZNZE,model="RKHS")) 
  # y = ZE + ZN + error
  
  Eta3<-list(E=list(X=ZE,model="BRR"),
             G=list(K=ZG,model="RKHS"),
             N=list(K=ZN,model="RKHS"),
             EG=list(K=ZGZE,model="RKHS"),
             EN=list(K=ZNZE,model="RKHS"))
  # y = ZE + ZA + ZD + ZE*ZA + ZD*ZE + error
  
  
  #### Cross-Validation ####
  hybrid<-as.character(unique(Pheno$Hybrid))
  Phenotype_data1<-Pheno
  
  
  set.seed(123)
  cycles<- 10
  CV2  <- list()
  CV1  <- list()
  CV0  <- list()
  CV00 <- list()
  
  
  for (ETAnum in 1:length(models)) {  
    if (ETAnum == 1) {
      ETAX <- Eta1
    } else if (ETAnum == 2) {
      ETAX <- Eta2
    } else {
      ETAX <- Eta3
    }
    
    
    for (i in 1:25) {
      folds <- sample(rep(c(1:5), times=29))
      df <- as.data.frame(cbind(hybrid, folds))
      
      x <- (i - 1) * 5
      
      for(r in 1:max(folds))
      {
        test_geno <- subset(df, grepl(r, folds))$hybrid
        train_geno <- subset(df, !grepl(r, folds))$hybrid
        
        CV_Data_1_2<-Phenotype_data1
        CV_Data_1_2$Y2<-NA
        CV_Data_1_2$Y2[CV_Data_1_2$Hybrid%in%train_geno]<-CV_Data_1_2$Blup[CV_Data_1_2$Hybrid%in%train_geno]
        
        y_t<-as.numeric(CV_Data_1_2$Y2)
        fit<-BGLR(y=y_t,ETA=ETAX,nIter=5000,burnIn=1000, thin=10)
        CV_Data_1_2$yhat <- fit$yHat
        
        # CV1 #
        df_test <- subset(CV_Data_1_2, CV_Data_1_2$Hybrid %in% test_geno)
        CV1[[(r+x)]] <- as.data.frame(df_test %>% group_by(Env) %>% dplyr::summarize(cor=cor(Blup, yhat,use = "complete.obs")))
        
        # CV2 #
        df_train <- subset(CV_Data_1_2, CV_Data_1_2$Hybrid %in% train_geno)
        CV2[[(r+x)]] <- as.data.frame(df_train %>% group_by(Env) %>% dplyr::summarize(cor=cor(Blup, yhat,use = "complete.obs")))
        
        CV_Data_00_0<-CV_Data_1_2
        CV_Data_00_0$Y3<-NA
        CV_Data_00_0$Y3[CV_Data_00_0$Hybrid%in%train_geno]<-CV_Data_00_0$Blup[CV_Data_00_0$Hybrid%in%train_geno]
        CV_Data_00_0$Y3[grepl("WS",CV_Data_00_0$Env )]<- NA # Train-test split based on WW vs. WS. To split by year, 
        # specify "2012" instead of "WS".
        y_t1<-as.numeric(CV_Data_00_0$Y3)
        fit2<-BGLR(y=y_t1,ETA=ETAX,nIter=5000,burnIn=1000, thin=10)
        CV_Data_00_0$yhat2 <- fit2$yHat
        
        # CV0 #
        df_train1 <- subset(CV_Data_00_0, CV_Data_00_0$Hybrid %in% train_geno)
        CV0[[(r+x)]] <- as.data.frame(df_train1 %>% group_by(Env) %>% dplyr::summarize(cor=cor(Blup, yhat2,use = "complete.obs")))
        
        # CV00 #
        df_test1 <- subset(CV_Data_00_0, CV_Data_00_0$Hybrid %in% test_geno)
        CV00[[(r+x)]] <- as.data.frame(df_test1 %>% group_by(Env) %>% dplyr::summarize(cor=cor(Blup, yhat2,use = "complete.obs")))
        
      }
      
      if (i == 25) {
        CV00out <- plyr::ldply(CV00, data.frame)
        CV0out <- plyr::ldply(CV0, data.frame)
        CV1out <- plyr::ldply(CV1, data.frame)
        CV2out <- plyr::ldply(CV2, data.frame)
        
        write.csv(CV00out, file = paste("ACC_", traitnames[tr], "_CV00_", models[ETAnum], ".csv", sep = ""), row.names = F)
        write.csv(CV0out, file = paste("ACC_", traitnames[tr], "_CV0_", models[ETAnum], ".csv", sep = ""), row.names = F)
        write.csv(CV1out, file = paste("ACC_", traitnames[tr], "_CV1_", models[ETAnum], ".csv", sep = ""), row.names = F)
        write.csv(CV2out, file = paste("ACC_", traitnames[tr], "_CV2_", models[ETAnum], ".csv", sep = ""), row.names = F)
      }
    }
  }
}

################################################################################

#### BGLR Prediction Script 2 ####

####
# This script is intended for use with traits that possess NIR data for ONLY THREE of four kernels (CS11_WW, CS12_WS, and CS12_WW). #
# The traits Phosphorus, Protein, and Starch do not have recorded NIR data for CS11_ws #
####

setwd("Set your directory here")

library(dplyr)
library(plyr)
library(AGHmatrix)
library(gplots)
library(dplyr)
library(reshape)
library(reshape2)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(data.table)
library(BGLR)


traitnames <- c("Phosphorus", "Protein", "Starch")
models <- c("Eta1", "Eta2", "Eta3")

for (tr in 1:length(traitnames)) {
  
  if (tr == 1) {
    Pheno <- read.csv(file = paste("Pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr == 2) {
    Pheno <- read.csv(file = paste("Pheno_", traitnames[tr], ".csv", sep = ""))
  } else {
    Pheno <- read.csv(file = paste("Pheno_", traitnames[tr], ".csv", sep = ""))
  }
  
  
  Pheno$Env<-as.factor(Pheno$Env)
  ZE<-model.matrix(~Pheno$Env-1)
  ZE[1:10,1:3]
  dim(ZE)
  Pheno$Hybrid<-as.factor(Pheno$Hybrid)
  ZVar<-model.matrix(~Pheno$Hybrid-1)
  dim(ZVar)
  
  
  #### Genomic effect kernel ####
  G <-as.matrix(read.csv("G.csv",row.names = 1))
  ZG<-ZVar%*%G%*%t(ZVar)
  ZG[1:5,1:5]
  
  #### NIR kernel with centering and scaling ####
  #NIR_CS11_WS <-as.matrix(read.csv("BLUP_CS11_WS.csv", row.names = 1))
  #NIR_CS11_WS_sc<-scale(NIR_CS11_WS,center=TRUE,scale=TRUE)
  
  NIR_CS11_WW <-as.matrix(read.csv("BLUP_CS11_WW.csv", row.names = 1))
  NIR_CS11_WW_sc<-scale(NIR_CS11_WW,center=TRUE,scale=TRUE)
  
  NIR_CS12_WS <-as.matrix(read.csv("BLUP_CS12_WS.csv", row.names = 1))
  NIR_CS12_WS_sc<-scale(NIR_CS12_WS,center=TRUE,scale=TRUE)
  
  NIR_CS12_WW <-as.matrix(read.csv("BLUP_CS12_WW.csv", row.names = 1))
  NIR_CS12_WW_sc<-scale(NIR_CS12_WW,center=TRUE,scale=TRUE)
  
  v <-  as.matrix(rbind(NIR_CS11_WW_sc,NIR_CS12_WS_sc,NIR_CS12_WW_sc))
  ZN <-tcrossprod(v)/ncol(v)
  
  
  #### Interaction kernels ####
  ZEZE<-tcrossprod(ZE)
  ZGZE<-ZG*ZEZE # Interaction kernel between genomic kernel and environment kernel
  ZNZE<-ZN*ZEZE # Interaction kernel between phenomic kernel and environment kernel
  
  
  Eta1<-list(E=list(X=ZE,model="BRR"),
             G=list(K=ZG,model="RKHS"),
             EG=list(K=ZGZE,model="RKHS")) 
  # y = ZE + ZG + error
  
  Eta2<-list(E=list(X=ZE,model="BRR"),
             N=list(K=ZN,model="RKHS"),
             EN=list(K=ZNZE,model="RKHS")) 
  # y = ZE + ZN + error
  
  Eta3<-list(E=list(X=ZE,model="BRR"),
             G=list(K=ZG,model="RKHS"),
             N=list(K=ZN,model="RKHS"),
             EG=list(K=ZGZE,model="RKHS"),
             EN=list(K=ZNZE,model="RKHS"))
  # y = ZE + ZA + ZD + ZE*ZA + ZD*ZE + error
  
  
  #### Cross-Validation ####
  hybrid<-as.character(unique(Pheno$Hybrid))
  Phenotype_data1<-Pheno
  
  #CV
  set.seed(123)
  cycles<- 10
  CV2  <- list()
  CV1  <- list()
  CV0  <- list()
  CV00 <- list()
  
  
  for (ETAnum in 1:length(models)) {  
    if (ETAnum == 1) {
      ETAX <- Eta1
    } else if (ETAnum == 2) {
      ETAX <- Eta2
    } else {
      ETAX <- Eta3
    }
    
    
    for (i in 1:25) {
      folds <- sample(rep(c(1:5), times=29))
      df <- as.data.frame(cbind(hybrid, folds))
      
      x <- (i - 1) * 5
      
      for(r in 1:max(folds))
      {
        test_geno <- subset(df, grepl(r, folds))$hybrid
        train_geno <- subset(df, !grepl(r, folds))$hybrid
        
        CV_Data_1_2<-Phenotype_data1
        CV_Data_1_2$Y2<-NA
        CV_Data_1_2$Y2[CV_Data_1_2$Hybrid%in%train_geno]<-CV_Data_1_2$Blup[CV_Data_1_2$Hybrid%in%train_geno]
        
        y_t<-as.numeric(CV_Data_1_2$Y2)
        fit<-BGLR(y=y_t,ETA=ETAX,nIter=5000,burnIn=1000, thin=10)
        CV_Data_1_2$yhat <- fit$yHat
        
        # CV1 #
        df_test <- subset(CV_Data_1_2, CV_Data_1_2$Hybrid %in% test_geno)
        CV1[[(r+x)]] <- as.data.frame(df_test %>% group_by(Env) %>% dplyr::summarize(cor=cor(Blup, yhat,use = "complete.obs")))
        
        # CV2 #
        df_train <- subset(CV_Data_1_2, CV_Data_1_2$Hybrid %in% train_geno)
        CV2[[(r+x)]] <- as.data.frame(df_train %>% group_by(Env) %>% dplyr::summarize(cor=cor(Blup, yhat,use = "complete.obs")))
        
        CV_Data_00_0<-CV_Data_1_2
        CV_Data_00_0$Y3<-NA
        CV_Data_00_0$Y3[CV_Data_00_0$Hybrid%in%train_geno]<-CV_Data_00_0$Blup[CV_Data_00_0$Hybrid%in%train_geno]
        CV_Data_00_0$Y3[grepl("WS",CV_Data_00_0$Env )]<- NA # Train-test split based on WW vs. WS. To split by year, 
        # specify "2012" instead of "WS".
        y_t1<-as.numeric(CV_Data_00_0$Y3)
        fit2<-BGLR(y=y_t1,ETA=ETAX,nIter=5000,burnIn=1000, thin=10)
        CV_Data_00_0$yhat2 <- fit2$yHat
        
        # CV0 #
        df_train1 <- subset(CV_Data_00_0, CV_Data_00_0$Hybrid %in% train_geno)
        CV0[[(r+x)]] <- as.data.frame(df_train1 %>% group_by(Env) %>% dplyr::summarize(cor=cor(Blup, yhat2,use = "complete.obs")))
        
        # CV00 #
        df_test1 <- subset(CV_Data_00_0, CV_Data_00_0$Hybrid %in% test_geno)
        CV00[[(r+x)]] <- as.data.frame(df_test1 %>% group_by(Env) %>% dplyr::summarize(cor=cor(Blup, yhat2,use = "complete.obs")))
        
      }
      
      if (i == 25) {
        CV00out <- plyr::ldply(CV00, data.frame)
        CV0out <- plyr::ldply(CV0, data.frame)
        CV1out <- plyr::ldply(CV1, data.frame)
        CV2out <- plyr::ldply(CV2, data.frame)
        
        write.csv(CV00out, file = paste("ACC_", traitnames[tr], "_CV00_", models[ETAnum], ".csv", sep = ""), row.names = F)
        write.csv(CV0out, file = paste("ACC_", traitnames[tr], "_CV0_", models[ETAnum], ".csv", sep = ""), row.names = F)
        write.csv(CV1out, file = paste("ACC_", traitnames[tr], "_CV1_", models[ETAnum], ".csv", sep = ""), row.names = F)
        write.csv(CV2out, file = paste("ACC_", traitnames[tr], "_CV2_", models[ETAnum], ".csv", sep = ""), row.names = F)
      }
    }
  }
}

################################################################################

#### Within-Environment Prediction ####

setwd("Set your directory here")

####
# To take advantage of the maximum number of hybrids available in each environment,
# within-environment prediction is shown below for prediction of agronomic traits using 
# as many as possible of the 199 hybrids in CS11_WS, 270 hybrids in CS11_WW, 319 hybrids 
# in CS12_WS, and 220 hybrids in CS12_WW. Note that 9 hybrids were not genotyped; these are
# CML382, Mp04:97, Mp07:121, Mp07:153, Mp07:158, NC312, NC408, S2B73, and W401.
####

ab <- read.table("SNP60000.hmp.txt", header = T)
ab[1:10,1:20]
myGAPIT <- GAPIT(G=ab, output.numerical=TRUE) # Hybrids: 346; SNPs: 61408
myGD <- myGAPIT$GD 
myGM <- myGAPIT$GM
myGD[1:5,1:5]
GBS <- myGD[,-1]; dim(GBS) # Numerical genotype-by-sequencing (GBS) data

#### Genomic and NIRS Matrices ####

# CS11_WS #
ws11 <-as.matrix(read.csv("PEDIGREE_CS11_WS.csv", row.names = 1)) # Import NIRS data for CS11_WS
ws11<-scale(ws11,center=TRUE,scale=TRUE) # Scale and center NIRS bands
c1 <- intersect(rownames(ws11), rownames(GBS))
GBS_ws11 <- subset(GBS, rownames(GBS) %in% c1) # Find subset of hybrids in NIRS data that match those in the GBS data
ws11 <- subset(ws11, rownames(ws11) %in% c1)
dim(ws11)
dim(GBS_ws11)
NIR_ws11 <- tcrossprod(ws11)/ncol(ws11); dim(NIR_ws11)
G_ws11 <- Gmatrix(SNPmatrix=as.matrix(GBS_ws11), missingValue="NA", maf=0.05, method="VanRaden"); dim(G_ws11)

# CS11_WW #
ww11 <-as.matrix(read.csv("PEDIGREE_CS11_WW.csv", row.names = 1)) # Import NIRS data for CS11_WW
ww11<-scale(ww11,center=TRUE,scale=TRUE) # Scale and center NIRS bands
c2 <- intersect(rownames(ww11), rownames(GBS))
GBS_ww11 <- subset(GBS, rownames(GBS) %in% c2) # Find subset of hybrids in NIRS data that match those in the GBS data
ww11 <- subset(ww11, rownames(ww11) %in% c2)
dim(ww11)
dim(GBS_ww11)
NIR_ww11 <- tcrossprod(ww11)/ncol(ww11)
G_ww11 <- Gmatrix(SNPmatrix=as.matrix(GBS_ww11), missingValue="NA", maf=0.05, method="VanRaden") 

# CS12_WS #
ws12 <-as.matrix(read.csv("PEDIGREE_CS12_WS.csv", row.names = 1)) # Import NIRS data for CS12_WS
ws12<-scale(ws12,center=TRUE,scale=TRUE) # Scale and center NIRS bands
c3 <- intersect(rownames(ws12), rownames(GBS))
GBS_ws12 <- subset(GBS, rownames(GBS) %in% c3) # Find subset of hybrids in NIRS data that match those in the GBS data
ws12 <- subset(ws12, rownames(ws12) %in% c3)
dim(ws12)
dim(GBS_ws12)
NIR_ws12 <- tcrossprod(ws12)/ncol(ws12)
G_ws12 <- Gmatrix(SNPmatrix=as.matrix(GBS_ws12), missingValue="NA", maf=0.05, method="VanRaden") 

# CS12_WW #
ww12 <-as.matrix(read.csv("PEDIGREE_CS12_WW.csv", row.names = 1)) # Import NIRS data for CS12_WW
ww12<-scale(ww12,center=TRUE,scale=TRUE) # Scale and center NIRS bands
c4 <- intersect(rownames(ww12), rownames(GBS))
GBS_ww12 <- subset(GBS, rownames(GBS) %in% c4) # Find subset of hybrids in NIRS data that match those in the GBS data
ww12 <- subset(ww12, rownames(ww12) %in% c4)
dim(ww12)
dim(GBS_ww12)
NIR_ww12 <- tcrossprod(ww12)/ncol(ww12)
G_ww12 <- Gmatrix(SNPmatrix=as.matrix(GBS_ww12), missingValue="NA", maf=0.05, method="VanRaden") 

#### Import Phenotypic Data (BLUPs) for Agronomic Traits ####
Pheno <- read.csv("BLUPs_Wide.csv")

Pheno_ws11 <- subset(Pheno, Pheno$Hybrid %in% c1 ,  select= grep("2011.WS", names(Pheno))) 
Pheno_ws11$Hybrid <- c1
Pheno_ww11 <- subset(Pheno, Pheno$Hybrid %in% c2 ,  select= grep("2011.WW", names(Pheno)))
Pheno_ww11$Hybrid <- c2
Pheno_ws12 <- subset(Pheno, Pheno$Hybrid %in% c3 ,  select= grep("2012.WS", names(Pheno))) 
Pheno_ws12$Hybrid <- c3
Pheno_ww12 <- subset(Pheno, Pheno$Hybrid %in% c4 ,  select= grep("2012.WW", names(Pheno))) 
Pheno_ww12$Hybrid <- c4

# Ensure all columns are in the same order for pheno files #
Pheno_ww11 <- Pheno_ww11[, c(1, 2, 3, 4, 6, 7, 11, 5, 8, 9, 10, 12)]
Pheno_ws12 <- Pheno_ws12[, c(1, 2, 3, 4, 6, 7, 11, 5, 8, 9, 10, 12)]
Pheno_ww12 <- Pheno_ww12[, c(1, 2, 3, 4, 6, 7, 11, 5, 8, 9, 10, 12)]

#### Create Model Kernels ####
# CS11_WS #
Eta1_ws11<-list(A=list(K=G_ws11,model="RKHS")) 
# y = G + error

Eta2_ws11<-list(N=list(K=NIR_ws11,model="RKHS")) 
# y = N + error

Eta3_ws11<-list(A=list(K=G_ws11,model="RKHS"),
                N=list(K=NIR_ws11,model="RKHS")) 
# y = G + N + error



# CS11_WW #
Eta1_ww11<-list(A=list(K=G_ww11,model="RKHS")) 
# y = G + error

Eta2_ww11<-list(N=list(K=NIR_ww11,model="RKHS")) 
# y = N + error

Eta3_ww11<-list(A=list(K=G_ww11,model="RKHS"),
                 N=list(K=NIR_ww11,model="RKHS")) 
# y = G + N + error



# CS12_WS #
Eta1_ws12<-list(A=list(K=G_ws12,model="RKHS")) 
# y = G + error

Eta2_ws12<-list(N=list(K=NIR_ws12,model="RKHS")) 
# y = N + error

Eta3_ws12<-list(A=list(K=G_ws12,model="RKHS"),
                 N=list(K=NIR_ws12,model="RKHS")) 
# y = G + N + error



# CS12_WW #
Eta1_ww12<-list(A=list(K=G_ww12,model="RKHS")) 
# y = G + error

Eta2_ww12<-list(N=list(K=NIR_ww12,model="RKHS")) 
# y = N + error

Eta3_ww12<-list(A=list(K=G_ww12,model="RKHS"),
                 N=list(K=NIR_ww12,model="RKHS")) 
# y = G + N + error



ETA_names <- c("Eta1_CS11_WS", "Eta2_CS11_WS", "Eta3_CS11_WS",
               "Eta1_CS11_WW", "Eta2_CS11_WW", "Eta3_CS11_WW",
               "Eta1_CS12_WS", "Eta2_CS12_WS", "Eta3_CS12_WS",
               "Eta1_CS12_WW", "Eta2_CS12_WW", "Eta3_CS12_WW")
ETA_list <- list(Eta1_ws11, Eta2_ws11, Eta3_ws11,
                 Eta1_ww11, Eta2_ww11, Eta3_ww11,
                 Eta1_ws12, Eta2_ws12, Eta3_ws12,
                 Eta1_ww12, Eta2_ww12, Eta3_ww12) # Create list of ETAs

pheno_data_list <- list(Pheno_ws11, Pheno_ww11, Pheno_ws12, Pheno_ww12)

#### Within-Environment Prediction ####
set.seed(789)
CV1 <- list()
CV2 <- list()

for (ETAnum in 1:length(ETA_list)){
  
  ETAX <- ETA_list[[ETAnum]] # ETA will change as ETAnum goes from 1 to 12
  
  if (ETAnum <= 3) {
    pheno_data <- pheno_data_list[[1]]
  } else if (ETAnum <= 6) {
    pheno_data <- pheno_data_list[[2]]
  } else if (ETAnum <= 9) {
    pheno_data <- pheno_data_list[[3]]
  } else if (ETAnum <= 12) {
    pheno_data <- pheno_data_list[[4]]
  }
  
  # Begin looping though pheno_data columns
  data_cols <- ncol(pheno_data)-1 # Number of columns with phenotypic data, used later for indexing
  total_cols <- as.numeric(ncol(pheno_data)) # Value of the last column
  
  for (col in 1:data_cols) {
    trait_n <- pheno_data[,c(col, total_cols)]
    hybrid <- trait_n$Hybrid
    
    for (rep_num in 1:25) {
      
      x <- (rep_num - 1) * 5 # Allows saving in lists as rep_num increases
      
      indices <- 1:nrow(trait_n) # Create a vector of indices from 1 to the number of observations
      k <- 5 # Set the number of folds
      obs_per_fold <- ceiling(nrow(trait_n) / k) # Calculate the number of observations per fold
      permuted_indices <- sample(indices) # Randomly permute the indices
      folds <- sample(rep(1:k, each = obs_per_fold, length.out = nrow(trait_n)))[permuted_indices] # Split the permuted 
      # indices into k approximately equal-sized subsets
      df <- as.data.frame(cbind(hybrid, folds))
      
      for(fold_num in 1:max(folds)) {
        
        test_geno <- subset(df, grepl(fold_num, folds))$hybrid
        train_geno <- subset(df, !grepl(fold_num, folds))$hybrid
        
        CV_Data_1_2 <- trait_n
        CV_Data_1_2$Y2<-NA
        CV_Data_1_2$Y2[CV_Data_1_2$Hybrid%in%train_geno] <- CV_Data_1_2[,1][CV_Data_1_2$Hybrid%in%train_geno]
        
        y_t<-as.numeric(CV_Data_1_2$Y2)
        fit<-BGLR(y=y_t,ETA=ETAX,nIter=5000,burnIn=1000, thin=10)
        CV_Data_1_2$yhat <- fit$yHat
        
        # CV1 #
        df_test <- subset(CV_Data_1_2, CV_Data_1_2$Hybrid %in% test_geno)
        CV1[[(fold_num+x)]] <- as.data.frame(df_test %>% dplyr::summarize(cor=cor(df_test[,1], yhat,use = "complete.obs")))
        
        
        # CV2 #
        df_train <- subset(CV_Data_1_2, CV_Data_1_2$Hybrid %in% train_geno)
        CV2[[(fold_num+x)]] <- as.data.frame(df_train %>% dplyr::summarize(cor=cor(df_train[,1], yhat,use = "complete.obs")))
        
      }
      
      if (rep_num == 25) {
        CV1out <- plyr::ldply(CV1, data.frame)
        CV2out <- plyr::ldply(CV2, data.frame)
        
        write.csv(CV1out, file = paste("ACC_", colnames(trait_n[1]), "_CV1_", ETA_names[ETAnum], ".csv", sep = ""), row.names = F)
        write.csv(CV2out, file = paste("ACC_", colnames(trait_n[1]), "_CV2_", ETA_names[ETAnum], ".csv", sep = ""), row.names = F)
      }
    }
  }
}



save.image("Save RData for future loading")


