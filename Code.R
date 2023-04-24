library(caTools)
library(ROSE)
library(cvTools)
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(survminer)
library(corrplot)
library(data.table)
library(riskRegression)
library(randtests)
library(DescTools)
library(lmtest)
library(survey)
library(rsample)
library(pec)
library(pROC)
library(ROCR)
library(Metrics)
library(foreign)
library(car)
library(coxme)
library(caret)
library(cluster)
library(factoextra)
library(NbClust)
library(clValid)
library(extrafont)
loadfonts(device = "win", quiet = TRUE)
op <- par(no.readonly = TRUE)
par(family = "Times New Roman")
windowsFonts(A = windowsFont("Times New Roman"))

#CLUSTER
datos<-read.csv(file=file.choose(),header=T,sep=';')#datos_cluster_ALL.csv
attach(datos)
head(datos)

#"1. Take a random sample for select the optimal number of clusters"
set.seed(101) 
sample <- sample.int(n = nrow(datos), size = floor(.1*nrow(datos)), replace = F)
k_clust <- datos[sample, ]

variables<-k_clust[,c(14:18)]


fviz_pca_ind(prcomp(variables), geom = "point")
fviz_nbclust(variables,kmeans,method = "silhouette")
fviz_nbclust(variables,kmeans,method = "wss")+geom_vline(xintercept = 3,linetype=2)


#"2. Build the k cluster in sample data"
kmeans.2 <- kmeans(variables, centers = 2)

#"3. Build the cluster with k groups in all data"
kmeans_2_all <- kmeans(datos[,c(9,10,21,22)], centers = 2)
datos<- cbind(datos, cluster2 = kmeans_2_all$cluster)

#4. "Analize the cluster conformation"
kruskal.test(Geology ~ cluster2, data = datos)
kruskal.test(Geomorphology ~ cluster2, data = datos)
kruskal.test(Slope_Rast ~ cluster2, data = datos)
kruskal.test(STREAM_DIST ~ cluster2, data = datos)

#MODEL: Cox
#1. "Select the data for training and testing"

datos_1 <- datos[datos$cluster2 == "1", c(1:23)]
set.seed(101) 
sample <- sample.int(n = nrow(datos_1), size = floor(.8*nrow(datos_1)), replace = F)
train <- datos_1[sample, ]
test  <- datos_1[-sample, ]

#2. "Data balance"

rl.both <- ovun.sample(Censure~Land_date+Geology+Geomorphology+STREAM_DIST
                       +Slope_Rast,
                       data=train,method = "both", p=0.5,seed=123,N=nrow(train))$data

#3. "Prepare the data for training"

Datos_training<-rl.both
attach(Datos_training)
head(Datos_training)

#4. "k-Folds CV"

k <- 10 #the number of folds
folds <- cvFolds(NROW(Datos_training), K=k)

#5. "The Cox model: training"
AUC_temp1 <-matrix(0,nrow = 5, ncol = 6)
n=0
p2 =c()
for(i in 1:k){
  train <- Datos_training[folds$subsets[folds$which != i], ] #Set the training set
  validation <- Datos_training[folds$subsets[folds$which == i], ] #Set the validation set
  #model_cox <- coxph(Surv(Land_date, Censure) ~STREAM_DIST+Slope_Rast+Geology+Geomorphology,ties="breslow", x = TRUE, y = TRUE, data=train)
  
  model_cox <- coxph(Surv(Land_date, Censure) ~Geomorphology+STREAM_DIST,ties="breslow", x = TRUE, y = TRUE, data=train)
  
  mcox_pred <- predictCox(model_cox, newdata=validation, times = 471)
  mcox_sup<-mcox_pred$survival #PROBABILIDADES DE SUPERVIVENCIA
  cox_predict<-1-mcox_sup
  for (j in 1:9){
    z<-seq(0.1,0.9, by=0.1)
    cox_predict2 <- ifelse(cox_predict > 0.5,1,0)
    p2 <- cbind(p2,cox_predict2 )
    n = n+1
    AUC_temp1[n] <- as.numeric(auc(validation$Censure,cox_predict2))
    print(c(n,as.numeric(auc(validation$Censure,cox_predict2))))
  }
}

model_cox
max(cox_predict)
apply(AUC_temp1, 1, FUN=mean)
max(apply(AUC_temp1, 1, FUN=mean))
cox.zph(model_cox)


#6. "RMSE"
Censure_validation=validation$Censure
RMSE_Cox<-rmse(Censure_validation, cox_predict)
RMSE_Cox

#7. "Confusion matrix"
valpredf_Cox<-as.factor(cox_predict2)
val_obsf_Cox<-as.factor(Censure_validation)
matriz_cox<-confusionMatrix(val_obsf_Cox, valpredf_Cox)
matriz_cox

#8. "AUC-ROC CURVE"

ROC_COX=roc(Censure_validation, cox_predict2,family="A", plot = TRUE, legacy.axes = TRUE,
            percent = TRUE, xlab = "False positives percentage",
            ylab = "True positives percentage", col = "#377eb8", lwd = 2,
            print.auc = TRUE, main="Cluster 1 Cox model: Training")


#9. Cox: testing"

rl.both_TEST <- ovun.sample(Censure~Land_date+Geology+Geomorphology+STREAM_DIST+Slope_Rast,data=test,method = "both", p=0.5,seed=123,N=nrow(test))$data

sum(rl.both_TEST$Censure)


Datos_testing<-rl.both_TEST
attach(Datos_testing)

k <- 10
folds <- cvFolds(NROW(Datos_testing), K=k)

AUC_temp4 <-matrix(0,nrow = 9, ncol = 10)
n=0

for(i in 1:k){
  train_TEST <- Datos_testing[folds$subsets[folds$which != i], ] 
  validation_TEST <- Datos_testing[folds$subsets[folds$which == i], ]
  model_cox_strat_TEST <- coxph(Surv(Land_date, Censure) ~ Geomorphology+STREAM_DIST,ties="breslow", x = TRUE, y = TRUE, data=train_TEST)
  mcox_strat_pred_TEST <- predictCox(model_cox_strat_TEST, newdata=validation_TEST, times = 471)
  mcox_strat_sup_TEST<-mcox_strat_pred_TEST$survival #PROBABILIDADES DE SUPERVIVENCIA
  cox_strat_predict_TEST<-1-mcox_strat_sup_TEST
  for (j in 1:9){
    z<-seq(0.1,0.9, by=0.1)
    cox_strat_predict2_TEST <- ifelse(cox_strat_predict_TEST > 0.5,1,0)
    p2 <- cbind(p2,log_predict2)
    n = n+1
    AUC_temp4[n] <- as.numeric(auc(validation_TEST$Censure,cox_strat_predict2_TEST))
    print(c(n,as.numeric(auc(validation$Censure,log_predict2))))
  }
}

model_cox_strat_TEST
apply(AUC_temp4, 1, FUN=mean)
max(apply(AUC_temp7, 1, FUN=mean))
cox.zph(model_cox_strat_TEST)
base_cox_cs <- basehaz(model_cox_strat_TEST)

Censure_validation_TEST=validation_TEST$Censure

roc(Censure_validation_TEST, family="A",cox_strat_predict_TEST, plot = TRUE, legacy.axes = TRUE,
    percent = TRUE, xlab = "False positives percentage",
    ylab = "True positives percentage", col = "#377eb8", lwd = 2,
    print.auc = TRUE, main="Cluster 1 Cox model: Testing")

RMSE_TEST<-rmse(Censure_validation_TEST, cox_strat_predict_TEST)
valpredf_TEST<-as.factor(cox_strat_predict2_TEST)
val_obsf_TEST<-as.factor(Censure_validation_TEST)
matriz_cox_TEST<-confusionMatrix(val_obsf_TEST, valpredf_TEST)
matriz_cox_TEST

#10. "COX PREDICTION cluster 1"
model_cox <- coxph(Surv(Land_date, Censure) ~STREAM_DIST+Geomorphology,ties="breslow", x = TRUE, y = TRUE, data=Datos_training) #Get your new linear model (just fit on the train data)
model_cox
mcox_pred <- predictCox(model_cox, newdata=datos_1, times = 471)
mcox_sup<-mcox_pred$survival 
cox_predict_ALL<-1-mcox_sup
for (j in 1:9){
  z<-seq(0.1,0.9, by=0.1)
  cox_predict2_ALL <- ifelse(cox_predict_ALL > 0.5,1,0)
}

datos_1<- cbind(datos_1, cox_predict_ALL)
write.csv(datos_1, file="datos_1.csv",row.names=FALSE)

#10. MODEL: Stratified Cox

model_cox <- coxph(Surv(Land_date, Censure) ~Geology+STREAM_DIST+Slope_Rast+strata(Geomorphology),ties="breslow", x = TRUE, y = TRUE, data=train)

mcox_pred <- predictCox(model_cox, newdata=validation, times = 471)
mcox_sup<-mcox_pred$survival
cox_predict<-1-mcox_sup
for (j in 1:9){
  z<-seq(0.1,0.9, by=0.1)
  cox_predict2 <- ifelse(cox_predict > 0.5,1,0)
  p2 <- cbind(p2,cox_predict2 )
  n = n+1
  AUC_temp1[n] <- as.numeric(auc(validation$Censure,cox_predict2))
  print(c(n,as.numeric(auc(validation$Censure,cox_predict2))))
}


model_cox
max(cox_predict)
apply(AUC_temp1, 1, FUN=mean)
max(apply(AUC_temp1, 1, FUN=mean))
cox.zph(model_cox)


#"RMSE"
Censure_validation=validation$Censure
RMSE_Cox<-rmse(Censure_validation, cox_predict)
RMSE_Cox

#"Confusion matrix"
valpredf_Cox<-as.factor(cox_predict2)
val_obsf_Cox<-as.factor(Censure_validation)
matriz_cox<-confusionMatrix(val_obsf_Cox, valpredf_Cox)
matriz_cox

#"AUC-ROC CURVE"

ROC_COX=roc(Censure_validation, cox_predict2,family="A", plot = TRUE, legacy.axes = TRUE,
            percent = TRUE, xlab = "False positives percentage",
            ylab = "True positives percentage", col = "#377eb8", lwd = 2,
            print.auc = TRUE, main="Cluster 1 stratified Cox model: Validation")


#"Stratified Cox: testing"

rl.both_TEST <- ovun.sample(Censure~Land_date+Geology+Geomorphology+STREAM_DIST+Slope_Rast,data=test,method = "both", p=0.5,seed=123,N=nrow(test))$data
sum(rl.both_TEST$Censure)
Datos_testing<-rl.both_TEST
attach(Datos_testing)
k <- 10
folds <- cvFolds(NROW(Datos_testing), K=k)
AUC_temp4 <-matrix(0,nrow = 9, ncol = 10)
n=0
for(i in 1:k){
  train_TEST <- Datos_testing[folds$subsets[folds$which != i], ] 
  validation_TEST <- Datos_testing[folds$subsets[folds$which == i], ]
  model_cox_strat_TEST <- coxph(Surv(Land_date, Censure) ~ Geology+STREAM_DIST+Slope_Rast+strata(Geomorphology),ties="breslow", x = TRUE, y = TRUE, data=train_TEST)
  mcox_strat_pred_TEST <- predictCox(model_cox_strat_TEST, newdata=validation_TEST, times = 471)
  mcox_strat_sup_TEST<-mcox_strat_pred_TEST$survival #PROBABILIDADES DE SUPERVIVENCIA
  cox_strat_predict_TEST<-1-mcox_strat_sup_TEST
  for (j in 1:9){
    z<-seq(0.1,0.9, by=0.1)
    cox_strat_predict2_TEST <- ifelse(cox_strat_predict_TEST > 0.5,1,0)
    p2 <- cbind(p2,log_predict2)
    n = n+1
    AUC_temp4[n] <- as.numeric(auc(validation_TEST$Censure,cox_strat_predict2_TEST))
    print(c(n,as.numeric(auc(validation$Censure,log_predict2))))
  }
}

model_cox_strat_TEST
apply(AUC_temp4, 1, FUN=mean)
max(apply(AUC_temp7, 1, FUN=mean)) # gana >0.5
cox.zph(model_cox_strat_TEST)
base_cox_cs <- basehaz(model_cox_strat_TEST)

Censure_validation_TEST=validation_TEST$Censure

roc(Censure_validation_TEST, family="A",cox_strat_predict_TEST, plot = TRUE, legacy.axes = TRUE,
    percent = TRUE, xlab = "False positives percentage",
    ylab = "True positives percentage", col = "#377eb8", lwd = 2,
    print.auc = TRUE, main="Cluster 1 stratified Cox model: Testing")

RMSE_TEST<-rmse(Censure_validation_TEST, cox_strat_predict_TEST)
valpredf_TEST<-as.factor(cox_strat_predict2_TEST)
val_obsf_TEST<-as.factor(Censure_validation_TEST)
matriz_cox_TEST<-confusionMatrix(val_obsf_TEST, valpredf_TEST)
matriz_cox_TEST

#"Stratified Cox PREDICTION cluster 1"
model_cox <- coxph(Surv(Land_date, Censure) ~Geology+STREAM_DIST+Slope_Rast+strata(Geomorphology),ties="breslow", x = TRUE, y = TRUE, data=Datos_training) #Get your new linear model (just fit on the train data)
model_cox
mcox_pred <- predictCox(model_cox, newdata=datos_1, times = 471)
mcox_sup<-mcox_pred$survival #PROBABILIDADES DE SUPERVIVENCIA
cox_predict_ALL<-1-mcox_sup
for (j in 1:9){
  z<-seq(0.1,0.9, by=0.1)
  cox_predict2_ALL <- ifelse(cox_predict_ALL > 0.5,1,0)
}

max(cox_predict_ALL)
datos<- cbind(datos_1, cluster2 = kmeans_2_all$cluster)

datos_1<- cbind(datos_1, cox_predict_ALL)
write.csv(datos_1, file="datos_1.csv",row.names=FALSE)

