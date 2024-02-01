rm(list=ls())
library(nlme)
library(sjstats)
library(lme4)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(gplots)
library(plyr)
library(viridis)
library(BGLR)
library(lessR)
library(devtools)
library(webr)


#### Multi-environment ANOVA for NIRS bands ####
df <- data.table::fread("df1.csv", header = T) %>% as.data.frame()

df$Range <- as.factor(df$Range)
df$Row <- as.factor(df$Row)
df$Rep <- as.factor(df$Rep)
df$Pedigree <- as.factor(df$Pedigree)
df$Env <- as.factor(df$Env)

varcomp <- list()
NIR <- list()

for (i in 81:3192)
{
  
  anova <- lmer(df[,i] ~
                  (1|Pedigree)+
                  (1|Env)+
                  (1|Pedigree:Env)+
                  (1|Range)+
                  (1|Row)+
                  (1|Rep:Env), df )
  anova
  rmse<- sjstats::rmse(anova)
  
  R<- MuMIn::r.squaredGLMM(anova)[,2]
  R
  VC<-as.data.frame(print(VarCorr(anova ), comp=c("Variance")))
  VC$Percent<-VC$vcov/sum(VC$vcov)
  VC
  heritability <- VC[2,6] / (VC[2,6] + (VC[1,6]/4) +   (VC[7,6]/8)     )
  round(heritability,3)  
  VC[,7] <-rmse
  VC[,8] <-heritability
  VC[,9] <-R
  VC[,10] <- colnames(df)[i]
  names(VC)[7:10] <- c("RMSE","Repeatability","Rsquared","Trait")
  
  varcomp[[i]]<- list(VC)
  
  E <- ranef(anova)$'Env'
  E[,2] <- rownames(E)
  names(E)[2] <- "Env"
  
  G <- ranef(anova)$'Pedigree'
  G[,2] <- rownames(G)
  names(G)[2] <- "Pedigree"
  
  GE <- coef(anova)$'Pedigree:Env'
  GE[,2] <- rownames(GE)
  GE$Pedigree <- lapply(strsplit(as.character(GE$V2), "\\:"), "[", 1)
  GE$Env <- lapply(strsplit(as.character(GE$V2), "\\:"), "[", 2)
  head(GE)
  GE <- as.data.frame(lapply(GE, unlist))
  head(GE)
  names(GE)[1:2] <- c( paste("X", colnames(df)[i], sep="" )  ,  "interaction" )
  head(GE)
  
  head(GE)
  head(E)
  head(G)
  
  GE <- left_join(GE,E, by="Env")
  GE <- left_join(GE,G, by="Pedigree")
  GE[,1] <- GE[,1]+GE[,5]+GE[,6]
  
  
  NIR[[i]]<- list(GE)
  
}


#### Save the variance components of each subset ####
varcomp <- varcomp[!is.na(varcomp)]
varcomp <- ldply(varcomp, data.frame)
write.csv(varcomp, "Var.comp.csv",row.names = F)
head(varcomp)
unique(varcomp$grp)

varcomp <- read.csv("Var.comp.csv")

varcomp$grp<- factor(varcomp$grp, levels = c( "Pedigree",
                                              "Env",
                                              "Pedigree:Env",
                                              "Range",
                                              "Row",
                                              "Rep:Env",
                                              "Residual"))


varcomp$Trait <- as.numeric(varcomp$Trait)

p <- varcomp %>% 
  ggplot( aes(x= Trait, y=Percent, fill=grp, text=grp)) +
  geom_area( ) +
  geom_line(aes(x= Trait, y=Repeatability))+
  scale_y_continuous("Explained percent variation by component (%)")+
  scale_x_continuous(bquote("NIRS Bands"~(cm^-1)),breaks = c(4000,5000,6000,7000,8000,9000,10000) ,guide = guide_axis(n.dodge=1))+
  theme_bw() +
  theme(legend.position=c(.8,.4),
        legend.background =  element_rect(fill = alpha("white", 0.7) , size=0.1, linetype="solid"))+
  guides(fill=guide_legend(title="Variance \ncomponent"))+ 
  ggtitle("A) Explained percent variation by component")

p
jpeg("Var Comp (NIR).jpeg", width = 8,height =5,units = "in", res=600)
p
dev.off()


#### Save the genotypic value results of each subset ####
NIR <- NIR[!is.na(NIR)]
NIR[sapply(NIR, is.null)] <- NULL

NIR <- Reduce(function(x, y) merge(x, y, by="interaction", all=TRUE), NIR)
dim(NIR)

write.csv(NIR, "NIR.csv",row.names = F)

NIR <- data.table::fread("NIR.stacked.txt") %>% as.data.frame()
NIR[1:5,1:5]

NIR$Band <- gsub('X', '',NIR$Band)
NIR$Band <- as.numeric(NIR$Band)
NIR$Year[NIR$Year == 2011] <- "CS11"
NIR$Year[NIR$Year == 2012] <- "CS12"

p2 <- ggplot(data=NIR, aes(x=Band, y=Data, group=Pedigree, color=interaction(Year,Trial))) +
  geom_line( ) +
  scale_y_continuous("Genotypic values of NIR bands")+
  scale_x_continuous(bquote("NIRS Bands"~(cm^-1)), breaks = c(4000,5000,6000,7000,8000,9000,10000),guide = guide_axis(n.dodge=1) )+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position=c(.9,.35),
        legend.background =  element_rect(fill = alpha("azure", 0.8) , size=0.1, linetype="solid"))+
  facet_grid(Year~Trial)+
  guides(color = FALSE)+
  ggtitle("B) Genotypic values of pedigree")
p2

jpeg("BLUP (NIR).jpeg", width = 8,height =5,units = "in", res=600)
p2
dev.off()


jpeg("Var.Comp.and.NIR_BLUP.jpeg",width = 6,height =8,units = "in", res=600)
grid.arrange(p, p2, ncol=1)
dev.off()


#### HapMap data ####
rm(list=ls())
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
myG <- data.table::fread("SNP60000.hmp.txt", head = F)
myG[1:5,1:15]

myGAPIT <- GAPIT(G=myG, 
                 output.numerical=TRUE,
                 file.output = F,
                 SNP.impute = "Major",
                 SNP.MAF = 0.05)
myGD= myGAPIT$GD
myGM= myGAPIT$GM
myGD[1:15,1:10]
names(myGD)[1] <- "Pedigree"
myGD <- myGD[order(myGD$Pedigree),]


#### Multi-environment ANOVA for KW and GY ####
rm(list=ls())


df <- data.table::fread("df1.csv", header = T) %>% as.data.frame()

df$Range <- as.factor(df$Range)
df$Row <- as.factor(df$Row)
df$Rep <- as.factor(df$Rep)
df$Pedigree <- as.factor(df$Pedigree)
df$Env <- as.factor(df$Env)

varcomp <- list()
NIR <- list()

anova <- lmer(df$kernel500Wt ~
                (1|Pedigree)+
                (1|Env)+
                (1|Pedigree:Env)+
                (1|Range)+
                (1|Row)+
                (1|Rep:Env), df )
anova
rmse<- sjstats::rmse(anova)

R<- MuMIn::r.squaredGLMM(anova)[,2]
R
VC<-as.data.frame(print(VarCorr(anova ), comp=c("Variance")))
VC$Percent<-VC$vcov/sum(VC$vcov)
VC
heritability <- VC[2,6] / (VC[2,6] + (VC[1,6]/4) +   (VC[7,6]/8)     )
round(heritability,3)  
VC[,7] <-rmse
VC[,8] <-heritability
VC[,9] <-R
VC[,10] <- "KW"
names(VC)[7:10] <- c("RMSE","Repeatability","Rsquared","Trait")
write.csv(VC,"varcomp.KW.csv")

E <- ranef(anova)$'Env'
E[,2] <- rownames(E)
names(E)[2] <- "Env"

G <- ranef(anova)$'Pedigree'
G[,2] <- rownames(G)
names(G)[2] <- "Pedigree"

GE <- coef(anova)$'Pedigree:Env'
GE[,2] <- rownames(GE)
GE$Pedigree <- lapply(strsplit(as.character(GE$V2), "\\:"), "[", 1)
GE$Env <- lapply(strsplit(as.character(GE$V2), "\\:"), "[", 2)
head(GE)
GE <- as.data.frame(lapply(GE, unlist))
head(GE)


head(GE)
head(E)
head(G)

GE <- left_join(GE,E, by="Env")
GE <- left_join(GE,G, by="Pedigree")
GE[,1] <- GE[,1]+GE[,5]+GE[,6]
write.csv(GE, "KW.BLUPs.csv")


df <- data.table::fread("df1.csv", header = T) %>% as.data.frame()

df$Range <- as.factor(df$Range)
df$Row <- as.factor(df$Row)
df$Rep <- as.factor(df$Rep)
df$Pedigree <- as.factor(df$Pedigree)
df$Env <- as.factor(df$Env)

anova <- lmer(df$yield.buac ~
                (1|Pedigree)+
                (1|Env)+
                (1|Pedigree:Env)+
                (1|Range)+
                (1|Row)+
                (1|Rep:Env), df )
anova
rmse<- sjstats::rmse(anova)

R<- MuMIn::r.squaredGLMM(anova)[,2]
R
VC<-as.data.frame(print(VarCorr(anova ), comp=c("Variance")))
VC$Percent<-VC$vcov/sum(VC$vcov)
VC
heritability <- VC[2,6] / (VC[2,6] + (VC[1,6]/4) +   (VC[7,6]/8)     )
round(heritability,3)  
VC[,7] <-rmse
VC[,8] <-heritability
VC[,9] <-R
VC[,10] <- "GY"
names(VC)[7:10] <- c("RMSE","Repeatability","Rsquared","Trait")
write.csv(VC,"varcomp.GY.csv")

E <- ranef(anova)$'Env'
E[,2] <- rownames(E)
names(E)[2] <- "Env"

G <- ranef(anova)$'Pedigree'
G[,2] <- rownames(G)
names(G)[2] <- "Pedigree"

GE <- coef(anova)$'Pedigree:Env'
GE[,2] <- rownames(GE)
GE$Pedigree <- lapply(strsplit(as.character(GE$V2), "\\:"), "[", 1)
GE$Env <- lapply(strsplit(as.character(GE$V2), "\\:"), "[", 2)
head(GE)
GE <- as.data.frame(lapply(GE, unlist))
head(GE)


head(GE)
head(E)
head(G)

GE <- left_join(GE,E, by="Env")
GE <- left_join(GE,G, by="Pedigree")
GE[,1] <- GE[,1]+GE[,5]+GE[,6]
write.csv(GE, "GY.BLUPs.csv")




VC.KW <- read.csv("varcomp.KW.csv")  
VC.GY <- read.csv("varcomp.GY.csv")   

VC <- rbind(VC.KW,VC.GY)

VC$grp<- factor(VC$grp, levels = c( "Pedigree","Env","Pedigree:Env", "Range","Row", "Rep:Env","Residual"))

p1 <-  ggplot(data=VC, aes(x= Trait, y=Percent)) +
  geom_col(mapping=aes(fill=grp), width= 0.7) +
  geom_point(mapping=aes(y =Repeatability,shape="Repeatability" ), fill="white", color="black",size=4, alpha=0.5) +
  geom_point(mapping=aes(y =Rsquared,shape="RSquared"),size=4,  alpha=1) +
  scale_shape_manual(values=c(1,2), name = NULL) +
  theme_bw() +
  theme(legend.position=c(.75,.4),
        legend.background =  element_rect(fill = alpha("azure", 0.5) , size=0.1, linetype="solid")) +
  guides(fill=guide_legend(title="Variance \ncomponent"))+ 
  ggtitle("A) Explained percent variation by component")

p1
jpeg("Var.Comp.GY.KW.jpeg",width = 5,height =6,units = "in", res=600)
p1
dev.off()

Blup.KW <- read.csv("KW.BLUPs.csv")
Blup.KW <- Blup.KW[,c(2,4,5,6)]
Blup.KW$Trait <- "KW (g)"

Blup.GY <- read.csv("GY.BLUPs.csv")
Blup.GY <- Blup.GY[,c(2,4,5,6)]
Blup.GY$Trait <- "GY (t/ha)"

Blup.GY$Data <- Blup.GY$Data*0.0673

Blups <- rbind(Blup.KW,Blup.GY)
head(Blups)

p2 <- ggplot(Blups, aes(x=as.factor(Year), y=Data, fill=Trial)) + 
  geom_boxplot(position=position_dodge(1),outlier.shape = NA)+
  facet_wrap(~Trait,scales = "free_y")+
  theme_bw()+
  scale_y_continuous("Genotypic value")+ 
  scale_x_discrete("Year")+
  ggtitle("B) Genotypic values of pedigree")+
  theme(legend.position=c(.9,.2))
p2

jpeg("Var.Comp.and.NIR_BLUPKW.GY.jpeg",width = 11,height =6,units = "in", res=600)
grid.arrange(p1, p2, ncol=2)
dev.off()


#### Correlation between NIRS bands and traits ####
cor <- data.table::fread("Correlation.between.band.and.trait.txt") %>% as.data.frame()
cor$Year <- ifelse(cor$Year == 2011, "CS11", "CS12")

p1 <-ggplot(cor, aes(Band , Correlation, color= interaction(Year,Trial) )) +
  geom_line(size=1)+ 
  theme_bw()+
  scale_x_continuous(bquote("NIRS Bands"~(cm^-1)), guide = guide_axis(n.dodge=1))+
  scale_y_continuous("Pearson correlation")+
  facet_grid(~Trait)+ 
  guides(color=guide_legend(title="Year.Trial"))+
  theme_bw() +
  theme(legend.position=c(.65,.2),
        legend.background =  element_rect(fill = alpha("lightblue", 0.4) , size=0.1, linetype="solid")) 

p1
jpeg("Correlation.jpeg",width = 8,height =5,units = "in", res=600)
p1
dev.off()


#### Data preparation before genomic/phenomic prediction ####
NIR <-  data.table::fread("NIR.csv") %>% as.data.frame()
NIR[1:10,1:10]
unique(NIR$Env)

Geno <-  data.table::fread("GAPIT.Genotype.Numerical.txt") %>% as.data.frame()
Geno[1:5,1:5]

sh <- intersect(Geno$taxa, NIR$Pedigree)


NIR <- NIR[NIR$Pedigree %in% sh,] 
Geno <- Geno[Geno$taxa %in% sh,]

length(unique(NIR$Pedigree))
length(Geno$taxa)

a <- NIR[NIR$Env=="CS11_WS",]
b <- NIR[NIR$Env=="CS11_WW",]
c <- NIR[NIR$Env=="CS12_WS",]
d <- NIR[NIR$Env=="CS12_WW",]


ch <- Reduce(intersect, list(Geno$taxa, a$Pedigree,b$Pedigree ,c$Pedigree, d$Pedigree)); length(ch)


#### Prediction ####
NIR[1:5,1:10]
Geno[1:5,1:5]
Geno <- Geno[order(Geno$taxa),]

# Genomic relationship matrix #
rownames(Geno) <- Geno$taxa
Geno <- Geno[,-1] 
Geno[1:5,1:5]

Geno.a <- scale(subset(Geno, row.names(Geno) %in% a$Pedigree),center=TRUE,scale=TRUE)
Geno.b <- scale(subset(Geno, row.names(Geno) %in% b$Pedigree),center=TRUE,scale=TRUE)
Geno.c <- scale(subset(Geno, row.names(Geno) %in% c$Pedigree),center=TRUE,scale=TRUE)
Geno.d <- scale(subset(Geno, row.names(Geno) %in% d$Pedigree),center=TRUE,scale=TRUE)

Geno.all <- rbind(Geno.a,Geno.b,Geno.c,Geno.d)
ZG <- tcrossprod(as.matrix(Geno.all))/ncol(as.matrix(Geno.all))
dim(ZG)

# Phenomic relationship matrix #
rownames(a) <- a$Pedigree
rownames(b) <- b$Pedigree
rownames(c) <- c$Pedigree
rownames(d) <- d$Pedigree

NIR.a <- scale(a[,-c(1:5)],center=TRUE,scale=TRUE)
NIR.b <- scale(b[,-c(1:5)],center=TRUE,scale=TRUE)
NIR.c <- scale(c[,-c(1:5)],center=TRUE,scale=TRUE)
NIR.d <- scale(d[,-c(1:5)],center=TRUE,scale=TRUE)


NIR.all <- rbind(NIR.a,NIR.b,NIR.c,NIR.d)
ZP <- tcrossprod(as.matrix(NIR.all))/ncol(as.matrix(NIR.all))
dim(ZP)


#### Lasso for NIRS band variable importance calculations ####
NIR[1:5,1:8]
library(caret)

control <- trainControl(method= "cv",
                        number=5,
                        verboseIter = T)


lasso_caret<- train(KW~.,
                    a[,-c(1,2,3,5)], 
                    method = "glmnet",
                    preProc = c("center", "scale"),
                    trControl=control,
                    tuneGrid = expand.grid(alpha = 1, lambda = 0))
lasso_caret

ib <- varImp(lasso_caret)$importance 

write.csv(ib, "CS11_WS_NIRbands.varImp.KW.csv")



lasso_caret<- train(KW~.,
                    NIR[,-c(1,2,3,5)], 
                    method = "glmnet",
                    preProc = c("center", "scale"),
                    trControl=control,
                    tuneGrid = expand.grid(alpha = 1, lambda = 0))

lasso_caret

ib <- varImp(lasso_caret)$importance 
head(ib)
write.csv(ib, "All.env_NIRbands.varImp.KW.csv")

varImp <- read.csv("All.varImp.csv")

q <- read.csv("All.env_NIRbands.varImp.KW.csv")
q$Trait <-"KW"
q$Env   <-"Combined"
w <- read.csv("CS11_WS_NIRbands.varImp.KW.csv")
w$Trait <-"KW"
w$Env   <-"2011.WS"
e <- read.csv("CS11_WW_NIRbands.varImp.KW.csv")
e$Trait <-"KW"
e$Env   <-"2011.WW"
r <- read.csv("CS12_WS_NIRbands.varImp.KW.csv")
r$Trait <-"KW"
r$Env   <-"2012.WS"
t <- read.csv("CS12_WW_NIRbands.varImp.KW.csv")
t$Trait <-"KW"
t$Env   <-"2012.WW"
y <- read.csv("All.env_NIRbands.varImp.GY.csv")
y$Trait <-"GY"
y$Env   <-"Combined"
u <- read.csv("CS11_WS_NIRbands.varImp.GY.csv")
u$Trait <- "GY"
u$Env   <- "2011.WS"
i <- read.csv("CS11_WW_NIRbands.varImp.GY.csv")
i$Trait <- "GY"
i$Env   <- "2011.WW"
o <- read.csv("CS12_WS_NIRbands.varImp.GY.csv")
o$Trait <- "GY"
o$Env   <- "2012.WS"
p <- read.csv("CS12_WW_NIRbands.varImp.GY.csv")
p$Trait <- "GY"
p$Env   <- "2012.WW"

varImp <- rbind(q,w,e,r,t,y,u,i,o,p)
names(varImp)[1] <- "Band"
head(varImp)
varImp$Band <- gsub("X", "" , varImp$Band)
varImp$Band <- as.numeric(varImp$Band)

p <- ggplot(varImp, aes(Band, Env, fill = Overall )) +
  geom_tile() +
  facet_grid(~Trait)+
  scale_fill_gradientn("Variable \nimportance", colors = hcl.colors(20, "RdYlGn"))
p

jpeg("varImp.jpeg",width = 8,height =3,units = "in", res=600)
p
dev.off()

write.csv(varImp, "all.varImp.csv")


#### Phenomic relationship matrix with subsetted bands based on lasso variable importance ####
varImp <- read.csv("all.varImp.csv")
head(varImp)
varImp.KW <- varImp[varImp$Trait=="KW",]
varImp.GY <- varImp[varImp$Trait=="GY",]

subset.bands <- subset(varImp,varImp.GY$Overall>0, )
subset.bands <- unique(subset.bands$Band)

NIR.all <- rbind(NIR.a,NIR.b,NIR.c,NIR.d)
NIR.all.subset <- NIR.all[,subset.bands]


ZP1 <- tcrossprod(as.matrix(NIR.all.subset))/ncol(as.matrix(NIR.all.subset))
dim(ZP1)


#### Environment kernel ####
Pheno <-rbind(a[,c(1,2,5)],b[,c(1,2,5)],c[,c(1,2,5)],d[,c(1,2,5)]) ## This is for GY
names(Pheno)[1] <- "hybrid"
head(Pheno)

Pheno$Env <- as.factor(Pheno$Env)
ZE<-model.matrix(~Pheno$Env-1)
ZE[1:10,1:4]
dim(ZE)


ZEZE   <- tcrossprod(ZE)
ZGZE   <- ZG*ZEZE
ZPZE   <- ZP*ZEZE
ZP1ZE   <- ZP1*ZEZE

Eta1<-list(E=list(X=ZE,model="BRR"),
           G=list(K=ZG,model="RKHS"),
           EG=list(K=ZGZE,model="RKHS")
)

Eta2<-list(E=list(X=ZE,model="BRR"),
           P=list(K=ZP,model="RKHS"),
           EP=list(K=ZPZE,model="RKHS")
)

Eta3<-list(E=list(X=ZE,model="BRR"),
           P1=list(K=ZP1,model="RKHS"),
           EP1=list(K=ZP1ZE,model="RKHS")
)

Eta4<-list(E=list(X=ZE,model="BRR"),
           G=list(K=ZG,model="RKHS"),
           P=list(K=ZP,model="RKHS"),
           EG=list(K=ZGZE,model="RKHS"),
           EP=list(K=ZPZE,model="RKHS")
) 

Eta5<-list(E=list(X=ZE,model="BRR"),
           G=list(K=ZG,model="RKHS"),
           P1=list(K=ZP1,model="RKHS"),
           EG=list(K=ZGZE,model="RKHS"),
           EP1=list(K=ZP1ZE,model="RKHS")
) 


#### Cross-Validation ####
hybrid<-ch ### These are the common hybrids determined above
length(ch)

CV2  <- list()
CV1  <- list()

CV0_CS11_WS  <- list()
CV00_CS11_WS <- list()

CV0_CS11_WW  <- list()
CV00_CS11_WW <- list()

CV0_CS12_WS  <- list()
CV00_CS12_WS <- list()

CV0_CS12_WW  <- list()
CV00_CS12_WW <- list()

Models <- list(Eta1, Eta2, Eta3, Eta4, Eta5)

library(parallel)
library(doParallel)
library(MASS)

numCores <- detectCores()
numCores
registerDoParallel(numCores-1)

library(BGLR)

foreach(MODEL = 1:length(Models), .packages = c("BGLR", "dplyr")) %dopar% {
  
  
  #for (MODEL in 1:length(Models))   {
  
  for (rep_num in 1:20)
  {
    
    x <- (rep_num - 1) * 5 # Allows saving in lists as rep_num increases
    
    k <- 5 # Set the number of folds
    obs_per_fold <- ceiling(length(hybrid)/k) # Calculate the number of observations per fold
    folds <-sample(rep(1:k, each = obs_per_fold, length.out = length(hybrid)))
    # Partition indices into k approximately equal-sized subsets
    df <- as.data.frame(cbind(hybrid, folds))
    
    for(fold_num in 1:max(folds))
    {
      
      #test_geno <- subset(df, grepl(fold_num, folds))$hybrid
      train_geno <- subset(df, !grepl(fold_num, folds))$hybrid
      
      
      yield <- Pheno 
      yield$Y2<-NA
      yield$Y2[yield$hybrid%in%train_geno] <- yield[,3][yield$hybrid%in%train_geno]
      y_t<-as.numeric(yield$Y2)
      
      fit<-BGLR(y=y_t,ETA=Models[[MODEL]],nIter=5000,burnIn=1000, thin=10)
      yield$yhat <- fit$yHat
      # CV1_M1 #
      df_test <- yield[!yield$hybrid %in% train_geno,]
      CV1[[(fold_num+x)]] <- df_test %>% group_by(Env) %>% dplyr::summarize(cor=cor(GY, yhat,use = "complete.obs")) %>% as.data.frame()
      # CV2_M1 #
      df_train <- yield[yield$hybrid %in% train_geno,]
      CV2[[(fold_num+x)]] <- df_train %>% group_by(Env) %>% dplyr::summarize(cor=cor(GY, yhat,use = "complete.obs")) %>% as.data.frame()
      
      ####################################
      ##### CS11_WS is untested env ######
      ####################################
      
      yield <- Pheno 
      yield$Y2<-NA
      yield$Y2[yield$hybrid%in%train_geno] <- yield[,3][yield$hybrid%in%train_geno]
      yield$Y2[yield$Env ==  'CS11_WS' ] <- NA 
      y_t<-as.numeric(yield$Y2)
      
      fit<-BGLR(y=y_t,ETA=Models[[MODEL]],nIter=5000,burnIn=1000, thin=10)
      yield$yhat <- fit$yHat
      # CV00 #
      df_test <- yield[!yield$hybrid %in% train_geno,]
      CV00_CS11_WS[[(fold_num+x)]] <- df_test %>% group_by(Env) %>% dplyr::summarize(cor=cor(GY, yhat,use = "complete.obs")) %>% as.data.frame()
      # CV0 #
      df_train <- yield[yield$hybrid %in% train_geno,]
      CV0_CS11_WS[[(fold_num+x)]] <- df_train %>% group_by(Env) %>% dplyr::summarize(cor=cor(GY, yhat,use = "complete.obs")) %>% as.data.frame()
      
      ####################################
      ##### CS11_WW is untested env ######
      ####################################
      
      yield <- Pheno 
      yield$Y2<-NA
      yield$Y2[yield$hybrid%in%train_geno] <- yield[,3][yield$hybrid%in%train_geno]
      yield$Y2[yield$year ==  'CS11_WW' ] <- NA 
      y_t<-as.numeric(yield$Y2)
      
      fit<-BGLR(y=y_t,ETA=Models[[MODEL]],nIter=5000,burnIn=1000, thin=10)
      yield$yhat <- fit$yHat
      # CV00 #
      df_test <- yield[!yield$hybrid %in% train_geno,]
      CV00_CS11_WW[[(fold_num+x)]] <- df_test %>% group_by(Env) %>% dplyr::summarize(cor=cor(GY, yhat,use = "complete.obs")) %>% as.data.frame()
      # CV0 #
      df_train <- yield[yield$hybrid %in% train_geno,]
      CV0_CS11_WW[[(fold_num+x)]] <- df_train %>% group_by(Env) %>% dplyr::summarize(cor=cor(GY, yhat,use = "complete.obs")) %>% as.data.frame()
      
      ####################################
      ##### CS12_WS is untested env ######
      ####################################
      
      yield <- Pheno 
      yield$Y2<-NA
      yield$Y2[yield$hybrid%in%train_geno] <- yield[,3][yield$hybrid%in%train_geno]
      yield$Y2[yield$Env == 'CS12_WS' ] <- NA 
      y_t<-as.numeric(yield$Y2)
      
      fit<-BGLR(y=y_t,ETA=Models[[MODEL]],nIter=5000,burnIn=1000, thin=10)
      yield$yhat <- fit$yHat
      # CV00 #
      df_test <- yield[!yield$hybrid %in% train_geno,]
      CV00_CS12_WS[[(fold_num+x)]] <- df_test %>% group_by(Env) %>% dplyr::summarize(cor=cor(GY, yhat,use = "complete.obs")) %>% as.data.frame()
      # CV0 #
      df_train <- yield[yield$hybrid %in% train_geno,]
      CV0_CS12_WS[[(fold_num+x)]] <- df_train %>% group_by(Env) %>% dplyr::summarize(cor=cor(GY, yhat,use = "complete.obs")) %>% as.data.frame()
      
      ####################################
      ##### CS12_WW is untested env ######
      ####################################
      
      yield <- Pheno 
      yield$Y2<-NA
      yield$Y2[yield$hybrid%in%train_geno] <- yield[,3][yield$hybrid%in%train_geno]
      yield$Y2[yield$year == 'CS12_WW' ] <- NA 
      y_t<-as.numeric(yield$Y2)
      
      fit<-BGLR(y=y_t,ETA=Models[[MODEL]],nIter=5000,burnIn=1000, thin=10)
      yield$yhat <- fit$yHat
      # CV00 #
      df_test <- yield[!yield$hybrid %in% train_geno,]
      CV00_CS12_WW[[(fold_num+x)]] <- df_test %>% group_by(Env) %>% dplyr::summarize(cor=cor(GY, yhat,use = "complete.obs")) %>% as.data.frame()
      # CV0 #
      df_train <- yield[yield$hybrid %in% train_geno,]
      CV0_CS12_WW[[(fold_num+x)]] <- df_train %>% group_by(Env) %>% dplyr::summarize(cor=cor(GY, yhat,use = "complete.obs")) %>% as.data.frame()
    }
    
    if (rep_num == 20){
      CV1out <- plyr::ldply(CV1, data.frame)
      CV2out <- plyr::ldply(CV2, data.frame)
      write.csv(CV1out,file = paste("CV1_", "GY_", "Eta", MODEL, ".csv", sep = ""), row.names = F)
      write.csv(CV2out,file = paste("CV2_", "GY_", "Eta", MODEL, ".csv", sep = ""), row.names = F)
      
      CV00out <- plyr::ldply(CV00_CS11_WS, data.frame)
      CV0out  <- plyr::ldply(CV0_CS11_WS, data.frame)
      write.csv(CV00out,file = paste("CV00_", "GY_", 'CS11_WS' , "_Eta", MODEL, ".csv", sep = ""), row.names = F)
      write.csv(CV0out,file = paste("CV0_", "GY_", 'CS11_WS' ,"_Eta", MODEL, ".csv", sep = ""), row.names = F)
      
      CV00out <- plyr::ldply(CV00_CS11_WW, data.frame)
      CV0out  <- plyr::ldply(CV0_CS11_WW, data.frame)
      write.csv(CV00out,file = paste("CV00_", "GY_", 'CS11_WW' , "_Eta", MODEL, ".csv", sep = ""), row.names = F)
      write.csv(CV0out,file = paste("CV0_", "GY_", 'CS11_WW' ,"_Eta", MODEL, ".csv", sep = ""), row.names = F)
      
      CV00out <- plyr::ldply(CV00_CS12_WS, data.frame)
      CV0out  <- plyr::ldply(CV0_CS12_WS, data.frame)
      write.csv(CV00out,file = paste("CV00_", "GY_", 'CS12_WS' , "_Eta", MODEL, ".csv", sep = ""), row.names = F)
      write.csv(CV0out,file = paste("CV0_", "GY_", 'CS12_WS' ,"_Eta", MODEL, ".csv", sep = ""), row.names = F)
      
      CV00out <- plyr::ldply(CV00_CS12_WW, data.frame)
      CV0out  <- plyr::ldply(CV0_CS12_WW, data.frame)
      write.csv(CV00out,file = paste("CV00_", "GY_", 'CS12_WW' , "_Eta", MODEL, ".csv", sep = ""), row.names = F)
      write.csv(CV0out,file = paste("CV0_", "GY_", 'CS12_WW' ,"_Eta", MODEL, ".csv", sep = ""), row.names = F)
      
    }
  }
}

#### Saving prediction abilities ####
library(readr)
list_csv_files <- list.files(path = "path_to_CV00_CV0_files")
df2 <- readr::read_csv(list_csv_files, id = "file_name") %>% as.data.frame()
df2
write.csv(df2, "Pred.ability.CV00.CV0.csv")


library(readr)
list_csv_files <- list.files(path = "path_to_CV1_CV2_files")
df2 <- readr::read_csv(list_csv_files, id = "file_name") %>% as.data.frame()
df2
write.csv(df2, "Pred.ability.CV1.CV2.csv")

df <- read.csv("Pred.ability.CV1.CV2.CV0.CV00.csv") 
head(df)
table(df$Env)

df<- df[df$Select=="1", ]

df <- as.data.frame(df %>%  dplyr::group_by(Model,CV,Trait) %>% 
                      dplyr::summarise(M = mean(cor, na.rm=TRUE),
                                       SD = sd(cor, na.rm=TRUE)))
head(df)
df$CV <- factor(df$CV, levels =  c("CV2", "CV1", "CV0", "CV00"))

p <- ggplot(df  , aes(Model, y=M, fill=Trait)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=round(M,2)  ), hjust=3, color="white",
            position = position_dodge(0.9), angle = 90,size=3.5)+
  geom_errorbar(aes(ymin=M, ymax=M+SD), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  facet_grid(~CV)+
  scale_y_continuous("Prediction ability")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p

jpeg("Pred.ability.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()