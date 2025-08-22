# Real data analysis for estimating treatment effects

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())
library(sigmoid)
library(MASS)

esti.INT <- function(rct){
  n0 <- nrow(rct)
  n0_treated <- sum(rct$X)
  n0_control <- n0 - n0_treated
  rct_treated <- rct[(rct$X==1), -1]
  rct_control <- rct[(rct$X==0), -1]
  
  propensity_model <- glm(X~., family = 'binomial', data = rct[, -2]) # glm(X~age+height+Weight+BMI+Marriage+Education+Gender+Job, family = "binomial", data = rct)
  propensity <- propensity_model$fitted.values
  
  reg1 <- glm(Outcome~., family = 'binomial', data = rct_treated)  # glm(Outcome~age+height+Weight+BMI+Marriage+Education+Gender+Job, family = "binomial", data = rct_treated)
  reg0 <- glm(Outcome~., family = 'binomial', data = rct_control)  # glm(Outcome~age+height+Weight+BMI+Marriage+Education+Gender+Job, family = "binomial", data = rct_control)
  m1 <- predict(reg1, newdata = rct[, -c(1,2)], type = "response") # predict(reg1, newdata = rct[,c('age', 'height', 'Weight', 'BMI', 'Marriage', 'Education', 'Gender', 'Job')], type = "response")
  m0 <- predict(reg0, newdata = rct[, -c(1,2)], type = "response") # predict(reg0, newdata = rct[,c('age', 'height', 'Weight', 'BMI', 'Marriage', 'Education', 'Gender', 'Job')], type = "response")
  
  hattau_1 <- mean(rct$X*(rct$Outcome - m1)/propensity + m1)
  hattau_0 <- mean((1 - rct$X)*(rct$Outcome - m0)/(1 - propensity) + m0)
  phi <- (rct$X*(rct$Outcome - m1)/propensity + m1) - 
    ((1 - rct$X)*(rct$Outcome - m0)/(1 - propensity) + m0)
  eta <- (1 - rct$X)*(rct$Outcome - m0)/(1 - propensity) + m0 - hattau_0
  
  internal <- hattau_1 - hattau_0
  inter_variance <- var(phi)/n0
  pvalue_int <- 1 - pnorm(internal/sqrt(inter_variance))
  c(internal, pvalue_int)
}

esti.ADF <- function(rct, tiltabeta, sigma, n1, C){
  n0 <- nrow(rct)
  n0_treated <- sum(rct$X)
  n0_control <- n0 - n0_treated
  rho <- n1/n0
  rct_treated <- rct[(rct$X==1), -1]
  rct_control <- rct[(rct$X==0), -1]
  
  propensity_model <- glm(X~., family = 'binomial', data = rct[, -2]) # glm(X~age+height+Weight+BMI+Marriage+Education+Gender+Job, family = "binomial", data = rct)
  propensity <- propensity_model$fitted.values
  
  reg1 <- glm(Outcome~., family = 'binomial', data = rct_treated)  # glm(Outcome~age+height+Weight+BMI+Marriage+Education+Gender+Job, family = "binomial", data = rct_treated)
  reg0 <- glm(Outcome~., family = 'binomial', data = rct_control)  # glm(Outcome~age+height+Weight+BMI+Marriage+Education+Gender+Job, family = "binomial", data = rct_control)
  m1 <- predict(reg1, newdata = rct[, -c(1,2)], type = "response") # predict(reg1, newdata = rct[,c('age', 'height', 'Weight', 'BMI', 'Marriage', 'Education', 'Gender', 'Job')], type = "response")
  m0 <- predict(reg0, newdata = rct[, -c(1,2)], type = "response") # predict(reg0, newdata = rct[,c('age', 'height', 'Weight', 'BMI', 'Marriage', 'Education', 'Gender', 'Job')], type = "response")
  
  hattau_1 <- mean(rct$X*(rct$Outcome - m1)/propensity + m1)
  hattau_0 <- mean((1 - rct$X)*(rct$Outcome - m0)/(1 - propensity) + m0)
  phi <- (rct$X*(rct$Outcome - m1)/propensity + m1) - 
    ((1 - rct$X)*(rct$Outcome - m0)/(1 - propensity) + m0)
  eta <- (1 - rct$X)*(rct$Outcome - m0)/(1 - propensity) + m0 - hattau_0
  
  S00 <- t(phi)%*%phi/n0^2
  S01 <- t(phi)%*%eta/n0^2
  S11 <- (sigma / rho + t(eta)%*%eta/n0)/n0
  
  wb <- pmax(1 - C * (hattau_0 - tildatau_0)^4 * sqrt(n0), 0)
  
  hattau <- hattau_1 - hattau_0 + S01 %*% wb %*% ginv(sqrt(wb) %*% S11 %*% sqrt(wb) + (1 - wb) * S11)%*%
    (tildatau_0 - hattau_0)
  
  variance <- t(phi)%*%phi/n0^{2} - S01 %*% wb %*% 
                     ginv(sqrt(wb) %*% S11 %*% sqrt(wb) + (1 - wb) * S11) %*% wb %*%
                     t(S01)
  
  stderror <- sqrt(variance)
  c(hattau, stderror)
}

crossvalidation <- function(dat,tildabeta,sigma,kfold,n0,n1){
  C.list <- 1:10
  shuffle.dat <- dat[sample(nrow(dat)),]
  subset <- cut(seq(1: nrow(dat)), breaks=kfold, labels=FALSE)
  
  cv.error <- rep(NA,length(C.list))
  for (j in 1:length(C.list)) {
    error <- rep(NA, kfold)
    for (i in 1:kfold) {
      #Segement your data by 'subset' using the which() function 
      testIndexes <- which(subset==i,arr.ind=TRUE)
      testData <- shuffle.dat[testIndexes, ]
      trainData <- shuffle.dat[-testIndexes, ]
      #Use the test and train data partitions however you desire...
      
      result_ADF <- esti.ADF(trainData,tildabeta, sigma,n1, C.list[j])
      htau_ADF <- result_ADF[1]
      
      result_INT <- esti.INT(testData)
      test_htau <- result_INT[1]
      error[i] <- mean((htau_ADF - test_htau)^2)
    }
    cv.error[j] <- mean(error)
  }
  C.optimal <- C.list[which.min(cv.error)]
}

boot.ADF.pval <- function(hattau1, b_try, S00, S01, S11, n0, n1, C, NB = 500){
  rho <- n1/n0
  
  Boot <- matrix(NA, length(hattau1), NB)
  
  for (j in 1:NB) {
    sumdata_tmp <- mvrnorm(1, c(hattau1, b_try), rbind(cbind(S00, - S01), 
                                                       cbind(t(- S01), S11)))
    tau_tmp <- sumdata_tmp[1:length(hattau1)]
    b_tmp <- sumdata_tmp[- (1:length(hattau1))]
    
    wb <- pmax(1 - C * b_tmp^4 * sqrt(n0), 0)
    
    hattau <- tau_tmp + S01 %*% wb %*% ginv(sqrt(wb) %*% S11 %*% sqrt(wb) + (1 - wb) * S11)%*%
      b_tmp
    
    Boot[, j] <- hattau - hattau1
  }
  
  mean(hattau1 < Boot)
}


obs_control <- read.csv("./HP_OBS_control.csv")
rct_treated <- read.csv("./HP_RCT_treated.csv")
rct_control <- read.csv("./HP_RCT_control.csv")


n0_treated <- length(rct_treated$Outcome)
n0_control <- length(rct_control$Outcome)
n0 <- n0_treated + n0_control
n1 <- length(obs_control$Outcome)
rho <- n1/n0


rct_treated$Education <- factor(rct_treated$Education)
rct_treated$Degree.of.stomach.pain.PRE <- factor(rct_treated$Degree.of.stomach.pain.PRE)
rct_control$Education <- factor(rct_control$Education)
rct_control$Degree.of.stomach.pain.PRE <- factor(rct_control$Degree.of.stomach.pain.PRE)

rct_treated$X <- rep(1, length(rct_treated$X))
rct_control$X <- rep(0, length(rct_control$X))
rct <- rbind(rct_treated, rct_control)
rct_treated <- rct[1:n0_treated, -1]
rct_control <- rct[-(1:n0_treated), -1]

## summary statistics
tildatau_0 <- mean(obs_control$Outcome)
sigma <- var(obs_control$Outcome) # mean((obs_control$Outcome - tildatau_0)^2) #

propensity_model <- glm(X~., family = 'binomial', data = rct[, -2]) # glm(X~age+height+Weight+BMI+Marriage+Education+Gender+Job, family = "binomial", data = rct)
propensity <- propensity_model$fitted.values

reg1 <- glm(Outcome~., family = 'binomial', data = rct_treated)  # glm(Outcome~age+height+Weight+BMI+Marriage+Education+Gender+Job, family = "binomial", data = rct_treated)
reg0 <- glm(Outcome~., family = 'binomial', data = rct_control)  # glm(Outcome~age+height+Weight+BMI+Marriage+Education+Gender+Job, family = "binomial", data = rct_control)
m1 <- predict(reg1, newdata = rct[, -c(1,2)], type = "response") # predict(reg1, newdata = rct[,c('age', 'height', 'Weight', 'BMI', 'Marriage', 'Education', 'Gender', 'Job')], type = "response")
m0 <- predict(reg0, newdata = rct[, -c(1,2)], type = "response") # predict(reg0, newdata = rct[,c('age', 'height', 'Weight', 'BMI', 'Marriage', 'Education', 'Gender', 'Job')], type = "response")

hattau_1 <- mean(rct$X*(rct$Outcome - m1)/propensity + m1)
hattau_0 <- mean((1 - rct$X)*(rct$Outcome - m0)/(1 - propensity) + m0)
phi <- (rct$X*(rct$Outcome - m1)/propensity + m1) - 
  ((1 - rct$X)*(rct$Outcome - m0)/(1 - propensity) + m0)
eta <- (1 - rct$X)*(rct$Outcome - m0)/(1 - propensity) + m0 - hattau_0

## Using only internal data
internal <- hattau_1 - hattau_0
inter_variance <- var(phi)/n0
pvalue_int <- 1 - pnorm(internal/sqrt(inter_variance))

## data-fused efficient estimator
EFF <- hattau_1 - hattau_0 - mean(phi*eta)/(sigma/rho + var(eta))*(hattau_0 - tildatau_0)
B <- var(phi) - mean(phi*eta)^2/(sigma/rho + var(eta))
EFF_variance <- B/n0
pvalue_EFF <- 1 - pnorm(EFF/sqrt(EFF_variance))

## crude estimator
PRM <- hattau_1 - tildatau_0
PRM_variance <- var(rct$X*(rct$Outcome - m1)/propensity + m1)/n0 + sigma/n1
pvalue_PRM <- 1 - pnorm(PRM/sqrt(PRM_variance))

## debiased estimator
C <- crossvalidation(rct,tildatau_0,sigma,kfold=3,n0,n1)
result_ADF <- esti.ADF(rct,tildatau_0,sigma, n1,C)
htau_ADF <- result_ADF[1]
pvalue_ADF <- 1 - pnorm(htau_ADF/result_ADF[2])

##bootstrapped ADF
set.seed(0)
KS <- 10
pval_boot <- c()
b_est <- tildatau_0 - hattau_0

hattau_1 - hattau_0 - mean(phi*eta)/(sigma/rho + var(eta))*(hattau_0 - tildatau_0)
S00 <- var(phi) / n0 
S11 <- (sigma/rho + var(eta)) / n0
S01 <- mean(phi*eta)/n0

p_h <- as.numeric(1 - pchisq(b_est^2 / S11, 1))

b_pos <- b_est
S_pos <- S11
S_pos <- S_pos
for (k in 1:KS) {
  b_try <- mvrnorm(1, b_pos, S_pos)
  b_try <- (p_h > 0.05 / log(n0)) * sqrt(pmax(b_pos^2 - S11, 0) / (b_pos^2 + S11)) * b_try + 
    (p_h <= 0.05 / log(n0)) * b_pos
  pval_boot[k] <- boot.ADF.pval(htau_ADF, b_try, S00, S01, S11, n0, n1, C, NB = 500)
}

max(pval_boot)
