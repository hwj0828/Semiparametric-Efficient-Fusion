# Simulation for Scenario II with biased summary statistics
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

library(parallel)
library(boot)
library(MASS)
library(gim)
library(expm)
library(glmnet)

# parallel replicate
# refer to https://rdrr.io/github/grayclhn/dbframe-R-library/man/RepParallel.html
# parallel replicate

source("mclapply.hack.R")

parReplicate <- function(n, expr, simplify = "array",...) {
  answer <-
    mclapply(integer(n), eval.parent(substitute(function(...) expr)),...)
  if (!identical(simplify, FALSE) && length(answer)) 
    return(simplify2array(answer, higher = (simplify == "array")))
  else return(answer)
}


esti.INT <- function(x,y){
  n0 <- length(y)
  result <- lm(y~x)
  design <- cbind(rep(1,n0), x)
  hattau1 <- result$coef[-1]   ## use only internal data
  # variance <- diag(solve(t(design)%*%design)%*%t(design)%*%
  #         diag(result$residual)^2%*%design%*%solve(t(design)%*%design))
  variance <- (mean(result$residual^2)*diag(solve(t(design)%*%design)))[-1]
  stderror <- sqrt(variance)
  c(hattau1, stderror)
}

esti.ORC <- function(x,y, tildabeta, V,n1){
  n0 <- length(y)
  n <- n0+n1
  rho <- n1/n0
  result <- lm(y~x)
  hattau1 <- result$coef[-1]
  
  tildabeta1 <- tildabeta[1]
  tildabeta2 <- tildabeta[2]
  reg1 <- lm(y~x[,1])
  hbeta1 <- reg1$coef[2]
  varX <- t(x)%*%x/n0
  phi <- (x*result$residuals)%*%solve(varX)
  eta1 <- x[,1]*reg1$residuals/varX[1,1]
  
  tmp <- mean(eta1^2)  ##Eeta1^2
  sigma2 <- mean(result$residuals^2)
  hattau <- hattau1 - (n1/(n0+n1)/tmp*sigma2*c(1/varX[1,1],0)*(hbeta1 - tildabeta1))
  variance <- diag(t(phi)%*%phi/n0 - n1/(n0+n1)*(t(phi)%*%eta1)%*%(t(eta1)%*%phi)/tmp/n0^2 )/n0
  stderror <- sqrt(variance)
  c(hattau, stderror)
}

esti.ADF <- function(x,y, tildabeta, V,n1, C){
  n0 <- length(y)
  n <- n0+n1
  rho <- n1/n0
  result <- lm(y~x)
  hattau1 <- result$coef[-1]
  
  tildabeta1 <- tildabeta[1]
  tildabeta2 <- tildabeta[2]
  reg1 <- lm(y~x[,1])
  hbeta1 <- reg1$coef[2]
  reg2 <- lm(y~x[,2])
  hbeta2 <- reg2$coef[2]
  hbeta <- c(hbeta1, hbeta2)
  
  varX <- t(x)%*%x/n0
  phi <- (x*result$residuals)%*%solve(varX)
  eta1 <- x[,1]*reg1$residuals/varX[1,1]
  eta2 <- x[,2]*reg2$residuals/varX[2,2]
  eta <- cbind(eta1, eta2)
  
  S00 <- t(phi)%*%phi/n0^2
  S01 <- t(phi)%*%eta/n0^2
  S11 <- (V / rho + t(eta)%*%eta/n0)/n0
  
  wb <- pmax(1 - C * (hbeta - tildabeta)^4 * sqrt(n0), 0)
  
  hattau <- hattau1 + S01 %*% diag(wb) %*% ginv(diag(sqrt(wb)) %*% S11 %*% diag(sqrt(wb)) + diag((1 - wb) * diag(S11)))%*%
    (tildabeta - hbeta)
  
  variance <- diag(t(phi)%*%phi/n0^{2} - S01 %*% diag(wb) %*% 
                     ginv(diag(sqrt(wb)) %*% S11 %*% diag(sqrt(wb)) + diag((1 - wb) * diag(S11))) %*% diag(wb) %*%
                     t(S01))
  
  stderror <- sqrt(variance)
  c(hattau, stderror)
}

crossvalidation <- function(x,y,tildabeta,V,kfold,n1){
  n0 <- length(y)
  C.list <- c(1 / (1:5), 1:5)
  cv.error <- rep(NA,length(C.list))
  S <- 50
  
  dat <- cbind(y,x)

  for (j in 1:length(C.list)) {
    error <- rep(NA, kfold*S)
    for (s in 1:S) {
      shuffle.dat <- dat[sample(nrow(dat)),]
      subset <- cut(seq(1: nrow(dat)), breaks=kfold, labels=FALSE)
      for (i in 1:kfold) {
        #Segement your data by 'subset' using the which() function 
        testIndexes <- which(subset==i,arr.ind=TRUE)
        testData <- shuffle.dat[testIndexes, ]
        trainData <- shuffle.dat[-testIndexes, ]
        #Use the test and train data partitions however you desire...
        
        result_ADF <- esti.ADF(x=trainData[,-1], y=trainData[,1], 
                               tildabeta, V,n1*(1 - 1/kfold), C.list[j])
        htau_ADF <- result_ADF[c(1,2)]
        
        result_INT <- esti.INT(testData[,-1], testData[,1])
        test_htau <- result_INT[c(1,2)]
        error[i+kfold*(s-1)] <- mean((htau_ADF - test_htau)^2)
      }
    }
    cv.error[j] <- mean(error)
  }
  C.optimal <- C.list[which.min(cv.error)]
}

esti.EFF <- function(x,y, tildabeta, V,n1){
  n0 <- length(y)
  n <- n0+n1
  rho <- n1/n0
  result <- lm(y~x)
  hattau1 <- result$coef[-1]
  
  tildabeta1 <- tildabeta[1]
  tildabeta2 <- tildabeta[2]
  reg1 <- lm(y~x[,1])
  hbeta1 <- reg1$coef[2]
  reg2 <- lm(y~x[,2])
  hbeta2 <- reg2$coef[2]
  varX <- t(x)%*%x/n0
  phi <- (x*result$residuals)%*%solve(varX)
  eta1 <- x[,1]*reg1$residuals/varX[1,1]
  eta2 <- x[,2]*reg2$residuals/varX[2,2]
  eta <- cbind(eta1, eta2)
  W <- t(eta)%*%eta/n0
  
  hattau <- hattau1 + (t(phi)%*%eta/n0)%*%solve(V/rho + W)%*%(
    c(tildabeta1 - hbeta1, tildabeta2 - hbeta2))
  variance <- diag(t(phi)%*%phi/n0 - (t(phi)%*%eta)%*%solve(t(eta)%*%eta)%*%(t(eta)%*%phi)/n0 )/n0
  stderror <- sqrt(variance)
  c(hattau, stderror)
}


esti.GIM <- function(x,y, tildabeta, V,n1){
  n0 <- length(y)
  dat0 <- as.data.frame(cbind(y,x))
  form1 <- 'y ~ x1'
  form2 <- 'y ~ x2'
  model <- list()
  ## partial information is available
  model[[1]] <- list(form = form1, 
                     info = data.frame(var = 'x1', 
                                       bet = tildabeta[1]))
  
  model[[2]] <- list(form = form2, 
                     info = data.frame(var = 'x2', 
                                       bet = tildabeta[2]))
  form <- 'y ~ x1+x2'
  nsample <- matrix(c(n1, n1, n1, n1),2, 2)
  skip_to_next <- FALSE
  tryCatch(fit <- gim(form, 'gaussian', dat0, model, nsample), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next){
    hattau <- c(NA, NA)
    variance <- c(NA, NA)
  } else{
    hattau <- fit$coef[-1]
    variance <- diag(fit$vcov)[-1]
  }
  stderror <- sqrt(variance)
  c(hattau, stderror)
}


simuOne <- function(n,n0,n1,kfold){
  tau <- c(1,1)
  kappa <- 0.6
  Sigma <- matrix(c(1,kappa, kappa, 1), nrow = 2, byrow = TRUE)
  rho <- n1/n0
  
  x <- mvrnorm(n, mu=c(0,0), Sigma = Sigma)
  y <- x%*%tau + rnorm(n,sd=2)
  dat <- as.data.frame(cbind(y,x))
  colnames(dat) <- c('y', 'x1', 'x2')
  
  
  dat0 <- dat[1:n0, ]
  dat1 <- dat[(n0+1):(n0+n1), ]
  dat1$x2 <- dat1$x2 + rnorm(n1,sd=1)
  
  x_inter <- as.matrix(dat0[,c(2,3)])
  y_inter <- dat0$y 
  
  x_exter <- as.matrix(dat1[,c(2,3)])
  y_exter <- dat1$y
  
  reg1 <- lm(y_exter~x_exter[,1])
  tildabeta1 <- reg1$coef[2]
  reg2 <- lm(y_exter~x_exter[,2])
  tildabeta2 <- reg2$coef[2]
  tildabeta <- c(tildabeta1, tildabeta2)
  varX_ext <- t(x_exter)%*%x_exter/n1
  tildaeta1 <- x_exter[,1]*reg1$residuals/varX_ext[1,1]
  tildaeta2 <- x_exter[,2]*reg2$residuals/varX_ext[2,2]
  tildaeta <- cbind(tildaeta1, tildaeta2)
  V <- t(tildaeta)%*%tildaeta/n1
  
  ## efficient estimator using only internal data
  result_INT <- esti.INT(x_inter, y_inter)
  htau_INT <- result_INT[c(1,2)]
  stderror_INT <- result_INT[c(3,4)]
  cv_INT1 <- (tau[1] >= htau_INT[1] - 1.96*stderror_INT[1]) & 
    (tau[1] <= htau_INT[1] + 1.96*stderror_INT[1])
  cv_INT2 <- (tau[2] >= htau_INT[2] - 1.96*stderror_INT[2]) & 
    (tau[2] <= htau_INT[2] + 1.96*stderror_INT[2])
  cv_INT <- c(cv_INT1, cv_INT2)
  
  ## oracle estimator
  result_ORC <- esti.ORC(x_inter, y_inter,tildabeta, V,n1)
  htau_ORC <- result_ORC[c(1,2)]
  stderror_ORC <- result_ORC[c(3,4)]
  cv_ORC1 <- (tau[1] >= htau_ORC[1] - 1.96*stderror_ORC[1]) & 
    (tau[1] <= htau_ORC[1] + 1.96*stderror_ORC[1])
  cv_ORC2 <- (tau[2] >= htau_ORC[2] - 1.96*stderror_ORC[2]) & 
    (tau[2] <= htau_ORC[2] + 1.96*stderror_ORC[2])
  cv_ORC <- c(cv_ORC1, cv_ORC2)
  
  ## debiased estimator
  C <- crossvalidation(x_inter,y_inter,tildabeta,V,kfold,n1)
  result_ADF <- esti.ADF(x_inter, y_inter, tildabeta, V,n1,C)
  htau_ADF <- result_ADF[c(1,2)]
  stderror_ADF <- result_ADF[c(3,4)]
  cv_ADF1 <- (tau[1] >= htau_ADF[1] - 1.96*stderror_ADF[1]) & 
    (tau[1] <= htau_ORC[1] + 1.96*stderror_ORC[1])
  cv_ADF2 <- (tau[2] >= htau_ADF[2] - 1.96*stderror_ADF[2]) & 
    (tau[2] <= htau_ADF[2] + 1.96*stderror_ADF[2])
  cv_ADF <- c(cv_ADF1, cv_ADF2)
  
  ## data-fused efficient estimator
  result_EFF <- esti.EFF(x_inter, y_inter,tildabeta, V,n1)
  htau_EFF <- result_EFF[c(1,2)]
  stderror_EFF <- result_EFF[c(3,4)]
  cv_EFF1 <- (tau[1] >= htau_EFF[1] - 1.96*stderror_EFF[1]) & 
    (tau[1] <= htau_EFF[1] + 1.96*stderror_EFF[1])
  cv_EFF2 <- (tau[2] >= htau_EFF[2] - 1.96*stderror_EFF[2]) & 
    (tau[2] <= htau_EFF[2] + 1.96*stderror_EFF[2])
  cv_EFF <- c(cv_EFF1, cv_EFF2)
  
  ## GIM
  result_GIM <- esti.GIM(x_inter, y_inter,tildabeta, V,n1)
  htau_GIM <- result_GIM[c(1,2)]
  stderror_GIM <- result_GIM[c(3,4)]
  cv_GIM1 <- (tau[1] >= htau_GIM[1] - 1.96*stderror_GIM[1]) & 
    (tau[1] <= htau_GIM[1] + 1.96*stderror_GIM[1])
  cv_GIM2 <- (tau[2] >= htau_GIM[2] - 1.96*stderror_GIM[2]) & 
    (tau[2] <= htau_GIM[2] + 1.96*stderror_GIM[2])
  cv_GIM <- c(cv_GIM1, cv_GIM2)
  
  list(htau_INT=htau_INT, htau_ORC=htau_ORC, htau_ADF=htau_ADF, 
       htau_EFF=htau_EFF, htau_GIM=htau_GIM,
       stderror_INT=stderror_INT, stderror_ORC=stderror_ORC, stderror_ADF=stderror_ADF,
       stderror_EFF=stderror_EFF, stderror_GIM=stderror_GIM,
       cv_INT=cv_INT, cv_ORC=cv_ORC, cv_ADF=cv_ADF, cv_EFF=cv_EFF, cv_GIM=cv_GIM)
}



nsims <- 1000
ncpus <- 20
n0 <- 500 #internal sample size in the maintext
n1 <- 2000 #external sample size in the main text
n <- n0 + n1
kfold <- 3

duration <- Sys.time()
duration

# set RNG seed for reproducibility 
RNGkind("L'Ecuyer-CMRG")
set.seed(777)

simu <- t(parReplicate(nsims, expr=simuOne(n, n0, n1, kfold = kfold), 
                       simplify=TRUE, mc.cores=ncpus, mc.set.seed=TRUE))

htau_INT <- do.call(rbind, simu[, 'htau_INT'])
htau_ORC <- do.call(rbind, simu[, 'htau_ORC'])
htau_ADF <- do.call(rbind, simu[, 'htau_ADF'])
htau_EFF <- do.call(rbind, simu[, 'htau_EFF'])
htau_GIM <- do.call(rbind, simu[, 'htau_GIM'])

stderror_INT <- do.call(rbind, simu[, 'stderror_INT'])
stderror_ORC <- do.call(rbind, simu[, 'stderror_ORC'])
stderror_ADF <- do.call(rbind, simu[, 'stderror_ADF'])
stderror_EFF <- do.call(rbind, simu[, 'stderror_EFF'])
stderror_GIM <- do.call(rbind, simu[, 'stderror_GIM'])

ASE_INT <- sqrt(apply(stderror_INT^2, 2, mean))
ASE_ORC <- sqrt(apply(stderror_ORC^2, 2, mean))
ASE_ADF <- sqrt(apply(stderror_ADF^2, 2, mean))
ASE_EFF <- sqrt(apply(stderror_EFF^2, 2, mean))
ASE_GIM <- sqrt(apply(stderror_GIM^2, 2, mean, na.rm=TRUE))

cv_INT <- do.call(rbind, simu[, 'cv_INT'])
cv_ORC <- do.call(rbind, simu[, 'cv_ORC'])
cv_ADF <- do.call(rbind, simu[, 'cv_ADF'])
cv_EFF <- do.call(rbind, simu[, 'cv_EFF'])
cv_GIM <- do.call(rbind, simu[, 'cv_GIM'])

cvp_INT <- apply(cv_INT, 2, mean)
cvp_ORC <- apply(cv_ORC, 2, mean)
cvp_ADF <- apply(cv_ADF, 2, mean)
cvp_EFF <- apply(cv_EFF, 2, mean)
cvp_GIM <- apply(cv_GIM, 2, mean, na.rm=TRUE)


duration <- Sys.time() - duration
duration

htau <- cbind(htau_INT, htau_ORC, htau_ADF, htau_EFF, htau_GIM)
MSE <- colMeans(sweep(htau, 2, rep(1,10))^2, na.rm=TRUE)

# save results
path.save <- paste0('Results/')
dir.create(path.save, recursive=TRUE)
path.image <- paste0(path.save,'/Simu2.RData')
path.out <- paste0(path.save,'/simu2.out')
save.image(file=path.image)


path.graph <- paste0(path.save,'/biased1.pdf')
pdf(file=path.graph, width=9, height=9)
op <- par(mfcol=c(1, 1), mai=c(0.6,0.5,0.6,0.4), 
          oma=c(2,2,2,2), cex=1)
boxplot(htau[,c(1,3,5,7,9)], outline=FALSE, xaxt='n', lwd=2, cex.axis=2.5,
        boxwex = 0.5)
title(main=expression(tau[1]),cex.main=4, line=2)
axis(1, at = 1:5, labels = c("INT", "ORC", "ADF", "EFF", "GIM"), line=1.5, tick=FALSE, cex.axis=2.5)
abline(h=1)
par(op)
dev.off()

path.graph <- paste0(path.save,'/biased2.pdf')
pdf(file=path.graph, width=9, height=9)
op <- par(mfcol=c(1, 1), mai=c(0.6,0.5,0.6,0.4), 
          oma=c(2,2,2,2), cex=1)
boxplot(htau[,c(2,4,6,8,10)], outline=FALSE, xaxt='n', lwd=2, cex.axis=2.5,
        boxwex = 0.5)
title(main=expression(tau[2]),cex.main=4, line=2)
axis(1, at = 1:5, labels = c("INT", "ORC", "ADF", "EFF", "GIM"), line=1.5, tick=FALSE, cex.axis=2.5)
abline(h=1)
par(op)
dev.off()


round(cbind(sqrt(MSE[seq(1, 10, 2)]), c(ASE_INT[1], ASE_ORC[1], ASE_ADF[1], ASE_EFF[1], ASE_GIM[1]), 
      c(cvp_INT[1], cvp_ORC[1], cvp_ADF[1], cvp_EFF[1], cvp_GIM[1]),
      sqrt(MSE[seq(2, 10, 2)]), c(ASE_INT[2], ASE_ORC[2], ASE_ADF[2], ASE_EFF[2], ASE_GIM[2]), 
      c(cvp_INT[2], cvp_ORC[2], cvp_ADF[2], cvp_EFF[2], cvp_GIM[2])) * 100, 2)

