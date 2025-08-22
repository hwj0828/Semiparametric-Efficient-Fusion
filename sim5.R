# Simulation for Scenario I
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

library(parallel)
library(boot)
library(MASS)


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


esti.INT <- function(x,y,A,n0,n1){
  n <- n0+n1
  ps <- glm(A~x, family = binomial(link =  "cauchit"))$fitted.values
  regressor <- cbind('Intercept' = rep(1,n0),'x' = x,
                     'A'=A, 'Ax'= A*x)
  oreg <- lm(y~0+regressor)$coef
  oreg1 <- c(oreg[1:2]+oreg[3:4], rep(0,2))
  oreg0 <- c(oreg[1:2], rep(0,2))
  tmp <- c(A*y/ps - (1 - A)*y/(1 - ps)
           -(A/ps - 1)*(regressor%*%oreg1) +
             ((1 - A)/(1 - ps) - 1)*(regressor%*%oreg0))
  ## efficient estimator using only internal data
  htau_INT <- mean(tmp)
  variance <- mean((tmp - htau_INT)^2)/n0
  stderror <- sqrt(variance)
  c(htau_INT, stderror)
}

esti.PRM <- function(x,y,A,hbeta1,n0,n1){
  n <- n0+n1
  ps <- glm(A~x, family = binomial(link =  "cauchit"))$fitted.values
  regressor <- cbind('Intercept' = rep(1,n0),'x' = x,
                     'A'=A, 'Ax'= A*x)
  oreg <- lm(y~0+regressor)$coef
  oreg1 <- c(oreg[1:2]+oreg[3:4], rep(0,2))
  oreg0 <- c(oreg[1:2], rep(0,2))
  tmp <- c(A*y/ps - (1 - A)*y/(1 - ps)
           -(A/ps - 1)*(regressor%*%oreg1) +
             ((1 - A)/(1 - ps) - 1)*(regressor%*%oreg0))
  residue <- c(y - cbind(rep(1,n0), x, A)%*%hbeta1)
  design <- cbind(rep(1,n0), x, A)
  htau_PRM <- mean(tmp - colMeans(c(tmp*residue)*design)%*%
                     solve(t(design)%*%diag(residue^2)%*%design/n0)%*%t(residue*design))
  variance <- (mean((tmp - htau_PRM)^2) + (n0/n1 - 1)*colMeans(c(tmp*residue)*design)%*%
                 solve(t(design)%*%diag(residue^2)%*%design/n0)%*%colMeans(c(tmp*residue)*design))/n0  
  stderror <- sqrt(variance)
  c(htau_PRM, stderror)
}

esti.EFF <- function(x,y,A,hbeta1,n0,n1){
  n <- n0+n1
  ps <- glm(A~x, family = binomial(link = "cauchit"))$fitted.values
  regressor <- cbind('Intercept' = rep(1,n0),'x' = x,
                     'A'=A, 'Ax'= A*x)
  oreg <- lm(y~0+regressor)$coef
  oreg1 <- c(oreg[1:2]+oreg[3:4], rep(0,2))
  oreg0 <- c(oreg[1:2], rep(0,2))
  tmp <- c(A*y/ps - (1 - A)*y/(1 - ps)
           -(A/ps - 1)*(regressor%*%oreg1) +
             ((1 - A)/(1 - ps) - 1)*(regressor%*%oreg0))
  lm_inter <- lm(y~x+A) 
  hbeta2 <- lm_inter$coef
  residue <- c(y - cbind(rep(1,n0), x, A)%*%hbeta2)
  design <- cbind(rep(1,n0), x, A)
  htau_EFF <- mean(tmp) - n1/n*colMeans(c(tmp*residue)*design)%*%
    solve(t(design)%*%diag(residue^2)%*%design/n0)%*%(t(design)%*%design/n0)%*%(hbeta2 - hbeta1)
  variance <- (mean((tmp - htau_EFF)^2) - n1/(n0+n1)*colMeans(c(tmp*residue)*design)%*%
                 solve(t(design)%*%diag(residue^2)%*%design/n0)%*%colMeans(c(tmp*residue)*design))/n0
  stderror <- sqrt(variance)
  c(htau_EFF, stderror)
}

esti.KNW <- function(x,y,A,beta,n0,n1){
  n <- n0+n1
  ps <- glm(A~x, family = binomial(link = "cauchit"))$fitted.values
  regressor <- cbind('Intercept' = rep(1,n0),'x' = x, 
                     'A'=A, 'Ax'= A*x)
  oreg <- lm(y~0+regressor)$coef
  oreg1 <- c(oreg[1:2]+oreg[3:4], rep(0,2))
  oreg0 <- c(oreg[1:2], rep(0,2))
  tmp <- c(A*y/ps - (1 - A)*y/(1 - ps)
           -(A/ps - 1)*(regressor%*%oreg1) +
             ((1 - A)/(1 - ps) - 1)*(regressor%*%oreg0))
  residue <- c(y - cbind(rep(1,n0), x, A)%*%beta)
  design <- cbind(rep(1,n0), x, A)
  htau_KNW <- mean(tmp - colMeans(c(tmp*residue)*design)%*%
                     solve(t(design)%*%diag(residue^2)%*%design/n0)%*%t(residue*design))
  variance <- (mean((tmp - htau_KNW)^2) - colMeans(c(tmp*residue)*design)%*%
                 solve(t(design)%*%diag(residue^2)%*%design/n0)%*%colMeans(c(tmp*residue)*design))/n0
  stderror <- sqrt(variance)
  c(htau_KNW, stderror)
}


simuOne <- function(n, n0, n1) {
  alpha <- c(1,1,1)
  tau <- 1
  beta <- c(1.2306, 0.6069, 0.5929)
  
  x <- rnorm(n, mean=0, sd=1)
  A <- rbinom(n,1, 1/(1+exp(-(1-x))))
  y1 <- cbind(1,x,x^2)%*%alpha + rnorm(n, mean=0, sd=2)
  y0 <- cbind(1,x,x^2)%*%(alpha - c(0,0,1)) + rnorm(n, mean=0, sd=1)
  y <- A*y1+(1-A)*y0
  
  x_inter <- x[1:n0]
  A_inter <- A[1:n0]
  y_inter <- y[1:n0]
  
  x_exter <- x[(n0+1):(n0+n1)]
  A_exter <- A[(n0+1):(n0+n1)]
  y_exter <- y[(n0+1):(n0+n1)]
  
  lm1 <- lm(y_exter~x_exter+A_exter)
  hbeta1 <- lm1$coef
  
  ## efficient estimator using only internal data
  result_INT <- esti.INT(x_inter, y_inter, A_inter,n0,n1)
  htau_INT <- result_INT[1]
  stderror_INT <- result_INT[2]
  cv_INT <- (tau >= htau_INT - 1.96*stderror_INT) & 
    (tau <= htau_INT + 1.96*stderror_INT)
  
  ## the crude estimator ignoring uncertainty of $\tilde{\beta}$
  result_PRM <- esti.PRM(x_inter, y_inter, A_inter, hbeta1,n0,n1)
  htau_PRM <- result_PRM[1]
  stderror_PRM <- result_PRM[2]
  cv_PRM <- (tau >= htau_PRM - 1.96*stderror_PRM) & 
    (tau <= htau_PRM + 1.96*stderror_PRM)
  
  ## data-fused efficient estimator
  result_EFF <- esti.EFF(x_inter, y_inter, A_inter, hbeta1,n0,n1)
  htau_EFF <- result_EFF[1]
  stderror_EFF <- result_EFF[2]
  cv_EFF <- (tau >= htau_EFF - 1.96*stderror_EFF) & 
    (tau <= htau_EFF + 1.96*stderror_EFF)
  
  ## the efficient estimator knowing the true value of $\beta$
  result_KNW <- esti.KNW(x_inter, y_inter, A_inter, beta,n0,n1)
  htau_KNW <- result_KNW[1]
  stderror_KNW <- result_KNW[2]
  cv_KNW <- (tau >= htau_KNW - 1.96*stderror_KNW) & 
    (tau <= htau_KNW + 1.96*stderror_KNW)
  
  list(htau_INT=htau_INT, htau_PRM=htau_PRM, htau_EFF=htau_EFF, htau_KNW=htau_KNW, 
       stderror_INT=stderror_INT, stderror_PRM=stderror_PRM, stderror_EFF=stderror_EFF, stderror_KNW=stderror_KNW,
       cv_INT=cv_INT, cv_PRM=cv_PRM, cv_EFF=cv_EFF, cv_KNW=cv_KNW)
}


nsims <- 1000
ncpus <- 20
n <- 3000 # = n0 + n1
n0 <- 1000 #internal sample size in the maintext
n1 <- 2000 #external sample size in the main text

duration <- Sys.time()
duration

# set RNG seed for reproducibility 
RNGkind("L'Ecuyer-CMRG")
set.seed(100)

simu <- t(parReplicate(nsims, expr=simuOne(n, n0, n1), 
                       simplify=TRUE, mc.cores=ncpus, mc.set.seed=TRUE))

htau_INT <- do.call(rbind, simu[, 'htau_INT'])
htau_PRM <- do.call(rbind, simu[, 'htau_PRM'])
htau_EFF <- do.call(rbind, simu[, 'htau_EFF'])
htau_KNW <- do.call(rbind, simu[, 'htau_KNW'])

stderror_INT <- do.call(rbind, simu[, 'stderror_INT'])
stderror_PRM <- do.call(rbind, simu[, 'stderror_PRM'])
stderror_EFF <- do.call(rbind, simu[, 'stderror_EFF'])
stderror_KNW <- do.call(rbind, simu[, 'stderror_KNW'])

ASE_INT <- sqrt(apply(stderror_INT^2, 2, mean))
ASE_PRM <- sqrt(apply(stderror_PRM^2, 2, mean))
ASE_EFF <- sqrt(apply(stderror_EFF^2, 2, mean))
ASE_KNW <- sqrt(apply(stderror_KNW^2, 2, mean))

cv_INT <- do.call(rbind, simu[, 'cv_INT'])
cv_PRM <- do.call(rbind, simu[, 'cv_PRM'])
cv_EFF <- do.call(rbind, simu[, 'cv_EFF'])
cv_KNW <- do.call(rbind, simu[, 'cv_KNW'])

cvp_INT <- apply(cv_INT, 2, mean)
cvp_PRM <- apply(cv_PRM, 2, mean)
cvp_EFF <- apply(cv_EFF, 2, mean)
cvp_KNW <- apply(cv_KNW, 2, mean)


duration <- Sys.time() - duration
duration

htau <- cbind(htau_INT, htau_PRM, htau_EFF, htau_KNW)

MSE <- colMeans(sweep(htau, 2, rep(1,4))^2)

RMSE <- sqrt(MSE)

# save results
path.save <- paste0('./Results/')
dir.create(path.save, recursive=TRUE)
path.image <- paste0(path.save,'/',n1,'Simu7.RData')
path.out <- paste0(path.save,'/simu.out')
save.image(file=path.image)


path.graph <- paste0(path.save,'/m=',n1,'model_misspecified2.pdf')
pdf(file=path.graph, width=9, height=9)
op <- par(mfcol=c(1, 1), mai=c(0.6,0.5,0.6,0.4), 
          oma=c(2,2,2,2), cex=1)
ylim.box=c(0.4, 1.6)
boxplot(htau, outline=FALSE, xaxt='n', lwd=2, cex.axis=2.5,
        ylim=ylim.box, boxwex = 0.5)
title(main=bquote("m ="~ .(n1)~"(Model misspecified)"),cex.main=3, line=2)
axis(1, at = 1:4, labels = c("INT", "PRM", "EFF", "KNW"), line=1.5, tick=FALSE, cex.axis=2.5)
abline(h=1)
par(op)
dev.off()
