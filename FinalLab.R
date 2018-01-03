library(dplyr)
library(ggplot2)
library(glmnet)
library(olsrr)
library(boot)
load('data/fMRIdata.RData')

#data preparation
#leave 20% for validation
partition_num = 1400
#prepare
fit_feat<-apply(fit_feat, 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
training_data <- fit_feat[1:partition_num,]
training_label <- resp_dat[1:partition_num,]
validation_data <- fit_feat[(partition_num+1):1750,]
validation_label <- resp_dat[(partition_num+1):1750,]

# Question 2

# Stepwise Forward Regression
# this may take a while
# you may like to skip running this part
Data <- data.frame(training_data)
Response <- data.frame(training_label)
Data_y1 <- cbind(Data,Response[,1])
colnames(Data_y1)[10922] <- 'y'
model <- lm(y ~ ., data = Data_y1)
k <- ols_step_forward(model, penter = 0.2, details = TRUE)
k <- ols_stepaic_forward(model, details = TRUE)

#CV
best_lambda <- matrix(0,20,1)
MSE <- matrix(0,20,1)
for (voxel in c(1:20)){
  cvfit = cv.glmnet(training_data, training_label[,voxel])
  lasso.cv.pred <- predict(cvfit, s = cvfit$lambda.min, training_data)
  best_lambda[voxel,1] <- cvfit$lambda.min
  MSE[voxel,1] <- mean((lasso.cv.pred - training_label[,voxel])^2)
}
mean(MSE)

# ESCV
nfold <- 10
index_all <- sample(1:1400,1400)
ESCV_estimate <- matrix(0,nfold,1400)
ESCV_ES_result <- matrix(0,100,20)
voxel <- 1 #internation

lambda_pool <- 10^seq(-1,-3, length = 100)
for (voxel in c(1:20)){
  i = 1
  for (lambda in lambda_pool){
    
    for (fold in c(1:nfold)){
      index_validation <- index_all[c((1+(fold-1)*1400/nfold):(140+(fold-1)*1400/nfold))]
      index_train <-index_all[-c((1+(fold-1)*1400/nfold):(140+(fold-1)*1400/nfold))]
      CV_train_data <- training_data[index_train,]
      CV_train_label <- training_label[index_train,voxel]
      CV_validation_data <- training_data[index_validation,]
      CV_validation_label <- training_label[index_validation,voxel]
      
      ESCV_lasso_mod <-glmnet(CV_train_data,CV_train_label,alpha = 1, lambda = lambda)
      ESCV_lasso_pred <- predict(ESCV_lasso_mod, s = lambda, newx = training_data)
      ESCV_estimate[fold,] <- ESCV_lasso_pred
    }
    
    ESCV_y_bar <- colMeans(ESCV_estimate)
    ESCV_y_var <- sum((ESCV_estimate[1,]-ESCV_y_bar)**2)
    for (fold in c(2:nfold)){
      ESCV_y_var <- ESCV_y_var + sum((ESCV_estimate[fold,]-ESCV_y_bar)**2)
    }
    ESCV_y_var <- ESCV_y_var/nfold
    ESCV_ES <- ESCV_y_var/sum(ESCV_y_bar**2)
    ESCV_ES_result[i] <- ESCV_ES
    i = i + 1
  }
}


#AIC AICc BIC
lambda_pool <- 10^seq(-1,-3, length = 50)
AIC_result <- matrix(0,50,20)
AICc_result <- matrix(0,50,20)
BIC_result <- matrix(0,50,20)
i = 1
n <- 1400
for (lambda in lambda_pool){
  for (voxel in c(1:20)){
    AIC_lasso_mod = glmnet(training_data, training_label[,voxel], alpha = 1, lambda = lambda)
    AIC_lasso_predict <- predict(AIC_lasso_mod, s = lambda, validation_data)
    k <- AIC_lasso_mod$beta@p[2]
    AIC_result[i,voxel] <- 2 * k + n * log(sum((AIC_lasso_predict-validation_label[,voxel])^2)/n, base = exp(1))
    AICc_result[i,voxel] <- AICc_result[i,voxel] + (2 * k * (k + 1))/(n - k - 1)
    BIC_result[i,voxel] <- n * log(sum((AIC_lasso_predict-validation_label[,voxel])^2)/n, base = exp(1)) + k * log(n, base = exp(1))
    print(c(voxel,i))
  }
  i = i + 1
}

#Figures
#to save your time, Please load the results
load('extra/Results/AIC.Rdata')
load('extra/Results/AICc.Rdata')
load('extra/Results/BIC.Rdata')
load('extra/Results/ESCV.Rdata')

#AIC
lambda <- 10^seq(-1,-3, length = 50)
AIC <- data.frame(AIC_result)
voxel <- 1
AICvoxel <- cbind(AIC[,voxel],lambda)
AICvoxel <- cbind(AICvoxel,rep(voxel,50))
AICtotal <- AICvoxel
for (voxel in c(2:20)){
  AICvoxel <- cbind(AIC[,voxel],lambda)
  AICvoxel <- cbind(AICvoxel,rep(voxel,50))
  AICtotal <- rbind(AICtotal,AICvoxel)
}
AICtotal <- data.frame(AICtotal)
colnames(AICtotal) <- c('AIC','lambda','Voxel')

pAIC <- ggplot(AICtotal, aes(x = log10(lambda), y = AIC))
pAIC <- pAIC + geom_line() +
  facet_wrap(~Voxel, ncol = 5)
pAIC + theme_minimal()+
  theme(text = element_text(color = "turquoise"))

#AICc
AICc <- data.frame(AICc_result)
voxel <- 1
AICcvoxel <- cbind(AICc[,voxel],lambda)
AICcvoxel <- cbind(AICcvoxel,rep(voxel,50))
AICctotal <- AICcvoxel
for (voxel in c(2:20)){
  AICcvoxel <- cbind(AICc[,voxel],lambda)
  AICcvoxel <- cbind(AICcvoxel,rep(voxel,50))
  AICctotal <- rbind(AICctotal,AICcvoxel)
}
AICctotal <- data.frame(AICctotal)
colnames(AICctotal) <- c('AICc','lambda','Voxel')

pAICc <- ggplot(AICctotal, aes(x = log(lambda), y = AICc))
pAICc <- pAICc + geom_line() +
  facet_wrap(~Voxel, ncol = 5)
pAICc + theme_minimal()+
  theme(text = element_text(color = "turquoise"))


#BIC
BIC <- data.frame(BIC_result)
voxel <- 1
BICvoxel <- cbind(BIC[,voxel],lambda)
BICvoxel <- cbind(BICvoxel,rep(voxel,50))
BICtotal <- BICvoxel
for (voxel in c(2:20)){
  BICvoxel <- cbind(BIC[,voxel],lambda)
  BICvoxel <- cbind(BICvoxel,rep(voxel,50))
  BICtotal <- rbind(BICtotal,BICvoxel)
}
BICtotal <- data.frame(BICtotal)
colnames(BICtotal) <- c('BIC','lambda','Voxel')

pBIC <- ggplot(BICtotal, aes(x = log(lambda), y = BIC))
pBIC <- pBIC + geom_line() +
  facet_wrap(~Voxel, ncol = 5)
pBIC + theme_minimal()+
  theme(text = element_text(color = "turquoise"))

#CV

#ESCV
lambda <- 10^seq(-0.5,-2, length = 100)
ESCV <- data.frame(ESCV)
voxel <- 1
ESCVvoxel <- cbind(ESCV[,voxel],lambda)
ESCVvoxel <- cbind(ESCVvoxel,rep(voxel,100))
ESCVtotal <- ESCVvoxel
for (voxel in c(2:19)){
  ESCVvoxel <- cbind(ESCV[,voxel],lambda)
  ESCVvoxel <- cbind(ESCVvoxel,rep(voxel,100))
  ESCVtotal <- rbind(ESCVtotal,ESCVvoxel)
}
ESCVtotal <- data.frame(ESCVtotal)
colnames(ESCVtotal) <- c('ESCV','lambda','Voxel')

pESCV <- ggplot(ESCVtotal, aes(x = lambda, y = ESCV))
pESCV <- pESCV + geom_line() +
  facet_wrap(~Voxel, ncol = 5)
pESCV + theme_minimal()+
  theme(text = element_text(color = "turquoise"))

ESCV1 <- ESCVtotal[ESCVtotal$Voxel == 20,]
ESCV1 <- ESCV1[ESCV1$ESCV < 1,]
p1 <- ggplot(ESCV1,aes(x = lambda, y = ESCV))
p1 <- p1 +geom_line()
p1

#Question 3

#Correlation between fitted values and observed values
load('extra/Results/ESCV.rdata')
lambda_pool <- 10^seq(-0.5,-2, length = 100)

cor_result <- array(0,dim=c(20,1))
lambda_best <- array(0,dim=c(20,1))
for (voxel in c(1:20)){
  lambda_min <- lambda_pool[which.min(ESCV[,voxel])]
  lambda_best[voxel] <- lambda_min
  lasso.mod <- glmnet(training_data, training_label[,voxel], alpha = 1, lambda = lambda_min)
  lasso.pred <- predict(lasso.mod, s = lambda_min, newx = validation_data)
  cor_result[voxel] <- cor(lasso.pred, validation_label[,voxel])
}


#Question 4
Stability <- array(0,dim=c(10,2,2))
for (t in c(1:10)){
  nfold <- 10
  index_all <- sample(1:1400,1400)
  training_data <- fit_feat[index_all,]
  training_label <- resp_dat[index_all,]
  validation_data <- fit_feat[-c(index_all),]
  validation_label <- resp_dat[-c(index_all),]
  ESCV_estimate <- matrix(0,nfold,1400)
  ESCV_ES_result <- matrix(0,20,20)
  lambda_pool <- 10^seq(-1,-2, length = 20)
  for (voxel in c(1,9)){
    i = 1
    for (lambda in lambda_pool){
      
      for (fold in c(1:nfold)){
        index_validation <- index_all[c((1+(fold-1)*1400/nfold):(140+(fold-1)*1400/nfold))]
        index_train <-index_all[-c((1+(fold-1)*1400/nfold):(140+(fold-1)*1400/nfold))]
        CV_train_data <- fit_feat[index_train,]
        CV_train_label <- resp_dat[index_train,voxel]
        CV_validation_data <- fit_feat[index_validation,]
        CV_validation_label <- resp_dat[index_validation,voxel]
        
        ESCV_lasso_mod <-glmnet(CV_train_data,CV_train_label,alpha = 1, lambda = lambda)
        ESCV_lasso_pred <- predict(ESCV_lasso_mod, s = lambda, newx = fit_feat[index_all,])
        ESCV_estimate[fold,] <- ESCV_lasso_pred
      }
      
      ESCV_y_bar <- colMeans(ESCV_estimate)
      ESCV_y_var <- sum((ESCV_estimate[1,]-ESCV_y_bar)**2)
      for (fold in c(2:nfold)){
        ESCV_y_var <- ESCV_y_var + sum((ESCV_estimate[fold,]-ESCV_y_bar)**2)
      }
      ESCV_y_var <- ESCV_y_var/nfold
      ESCV_ES <- ESCV_y_var/sum(ESCV_y_bar**2)
      ESCV_ES_result[i,voxel] <- ESCV_ES
      i = i + 1
      print(c(voxel,i))
    }
  }
  
  voxel <- 9
  voxel9 <- ESCV_ES_result[,voxel]
  lambda_voxel_9 <- lambda_pool[which.min(voxel9)]
  lasso.mod <- glmnet(training_data, training_label[,voxel], alpha = 1, lambda = lambda_voxel_9)
  lasso.pred <- predict(lasso.mod, s = lambda_voxel_9, newx = validation_data)
  Stability[t,2,1] <- cor(lasso.pred,validation_label[,voxel])
  Stability[t,2,2] <- mean((lasso.pred-validation_label[,voxel])^2)
  voxel <- 1
  voxel1 <- ESCV_ES_result[,voxel]
  lambda_voxel_1 <- lambda_pool[which.min(voxel1)]
  
  lasso.mod <- glmnet(training_data, training_label[,voxel], alpha = 1, lambda = lambda_voxel_1)
  lasso.pred <- predict(lasso.mod, s = lambda_voxel_1, newx = validation_data)
  Stability[t,1,1] <- cor(lasso.pred,validation_label[,voxel])
  Stability[t,1,2] <- mean((lasso.pred-validation_label[,voxel])^2)
}

#stability of CV
stabilityCV <- array(0,dim=c(10,2,2))
for (t in c(1:10)){
  best_lambda <- matrix(0,10,1)
  MSE <- matrix(0,10,1)
  index_all <- sample(1:1400,1400)
  CV_training <- fit_feat[index_all,]
  CV_training_l <- resp_dat[index_all,]
  CV_validation <- fit_feat[-c(index_all),]
  CV_validation_l <- resp_dat[-c(index_all),]
  for (voxel in c(1,9)){
    cvfit = cv.glmnet(CV_training, CV_training_l[,voxel])
    lasso.cv.pred <- predict(cvfit, s = cvfit$lambda.min, CV_training)
    best_lambda[voxel,1] <- cvfit$lambda.min
    MSE[voxel,1] <- mean((lasso.cv.pred-CV_training_l[,voxel])^2)
  }
  voxel <- 9
  lambda_voxel_9 <- best_lambda[voxel]
  lasso.mod <- glmnet(CV_training, CV_training_l[,voxel], alpha = 1, lambda = lambda_voxel_9)
  lasso.pred <- predict(lasso.mod, s = lambda_voxel_9, newx = CV_validation)
  stabilityCV[t,2,1] <- cor(lasso.pred,CV_validation_l[,voxel])
  stabilityCV[t,2,2] <- mean((lasso.pred-CV_validation_l[,voxel])^2)
  
  voxel <- 1
  lambda_voxel_1 <- best_lambda[voxel]
  lasso.mod <- glmnet(CV_training, CV_training_l[,voxel], alpha = 1, lambda = lambda_voxel_1)
  lasso.pred <- predict(lasso.mod, s = lambda_voxel_1, newx = CV_validation)
  stabilityCV[t,1,1] <- cor(lasso.pred,CV_validation_l[,voxel])
  stabilityCV[t,1,2] <- mean((lasso.pred-CV_validation_l[,voxel])^2)
}
voxel1 <- stabilityCV[,1,]
voxel9 <- stabilityCV[,2,]
apply(voxel1,2,sd)
apply(voxel9,2,sd)

#Question 5
load('extra/Results/ESCV.rdata')
fit_stim <- read.csv('data/fit_stim.csv')
real_wav <- read.csv('data/real_wav.csv')

lambda_pool <- 10^seq(-0.5,-2, length = 100)
voxel <- 13
Voxel9 <- ESCV[,voxel]
lambda_voxel_9 <- lambda_pool[which.min(Voxel9)]

lasso.mod <- glmnet(training_data, training_label[,voxel], alpha = 1, lambda = lambda_voxel_9)
lasso.predict <- predict(lasso.mod, s = lambda_voxel_9, newx = validation_data)
cor(lasso.predict, validation_label[,voxel])
coef_voxel9 <- coef(lasso.mod)

coef_non0_voxel9 <- coef(lasso.mod)[which(coef(lasso.mod) != 0)]
coef_non0_voxel9_X <- which(coef(lasso.mod) != 0)-1

real_wav_voxel9 <- real_wav[,coef_non0_voxel9_X]
real_wav_voxel9_sum <- apply(real_wav_voxel9,1,sum)
real_wav_voxel9_location <- which(real_wav_voxel9_sum !=0)
x <- rep(1:128,128)
y <- rep(1:128, each = 128)
value <- rep(0, 128*128)
image <- cbind(x,y)
image <- cbind(image,value)
image <- data.frame(image)
colnames(image) <- c('x','y','sensitiveness')
for (i in real_wav_voxel9_location){
  image[i,3] <- real_wav_voxel9_sum[i]
}

ggplot(image) + geom_point(aes(x = x, y = y, color = sensitiveness)) +
  ggtitle("The response to images of Voxel 13")

colnames(image) <- c('x','y','region')
for (i in real_wav_voxel9_location){
  image[i,3] <- 1
}

ggplot(image) + geom_point(aes(x = x, y = y, color = region)) +
  ggtitle("The response region to images of Voxel 13")

#Importance
importance_order <- order(coef_non0_voxel9)
importance <- array(0,dim=c(length(coef_non0_voxel9),1))
j=0
for (i in importance_order){
  importance[length(coef_non0_voxel9)-j] <- coef_non0_voxel9_X[i]
  j = j + 1
}

#bootstrp
# Bootstrap 95% CI for regression coefficients 


Data <- data.frame(fit_feat)
Response <- data.frame(resp_dat)
Data_y1 <- cbind(Data,Response[,1])
colnames(Data_y1)[10922] <- 'y'

bs <- function(formula, data, indices) {
  d <- data[indices,] 
  fit <- lm(formula, data=d)
  return(coef(fit)) 
} 
# bootstrapping with 1000 replications 
results <- boot(data=Data_y1, statistic=bs, 
                R=1000, formula=y ~.)

# view results
results
plot(results, index=1) 
plot(results, index=2) 
plot(results, index=3) 

# get 95% confidence intervals 
boot.ci(results, type="bca", index=1) 
boot.ci(results, type="bca", index=2) 
boot.ci(results, type="bca", index=3) 

#Question 6
val_feat <- apply(val_feat, 2, function(y) (y-mean(y)) / sd(y) ^ as.logical(sd(y)))
prediction <- array(0,dim=c(120,20))
for (voxel in c(1:20)){
  lambda_min <- lambda_pool[which.min(ESCV[,voxel])]
  lambda_best[voxel] <- lambda_min
  lasso.mod <- glmnet(training_data, training_label[,voxel], alpha = 1, lambda = lambda_min)
  lasso.pred <- predict(lasso.mod, s = lambda_min, newx = val_feat)
  prediction[,voxel] <- lasso.pred
}

write.table(prediction, file="predv1_Hongxu.txt", row.names=FALSE, col.names=FALSE)
