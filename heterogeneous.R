# Sept 8, 2022
# Simulation for heterogeneous variance of Example 3 in supplementary material

rm(list = ls())
library(MASS)
theta <- c(1,1)
kappa <- 0.5
Sigma <- matrix(c(1,kappa, kappa, 1), nrow = 2, byrow = TRUE)
n = 10000
n0 <- 1000
n1 <- 2000 
rho <- n1/n0
set.seed(100)

iteration <- 1000
htheta <- matrix(NA, ncol = 4, nrow = iteration)

for (i in 1:iteration) {
  x <- mvrnorm(n, mu=c(0,0), Sigma = Sigma)
  sd <- sqrt(1 + abs(3*x[,2]+ x[,1])^2 ) 
  y <- x%*%theta + rnorm(n, sd = sd)
  dat <- as.data.frame(cbind(y,x))
  colnames(dat) <- c('y', 'x1', 'x2')
  
  
  dat0 <- dat[1:n0, ]
  dat1 <- dat[(n0+1):(n0+n1), ]
  
  sd <- sd[1:n0]
  x_inter <- as.matrix(dat0[,c(2,3)])
  y_inter <- dat0$y 
  
  x_exter <- as.matrix(dat1[,c(2,3)])
  y_exter <- dat1$y 
  
  result_inter <- lm(y_inter~x_inter, weights = 1/sd^2)
  ## use only internal data
  hattheta1 <- result_inter$coef[-1]   
  
  ## data-fused efficient estimator
  reg1 <- lm(y_exter~x_exter[,1])
  tildabeta1 <- reg1$coef[2]

  reg3 <- lm(y_inter~x_inter[,1])
  hbeta1 <- reg3$coef[2]

  varX <- t(x_exter)%*%x_exter/n1
  tildaeta1 <- x_exter[,1]*reg1$residuals/varX[1,1]
  V <- mean(tildaeta1^2)

  varX_in <- t(x_inter)%*%x_inter/n0
  eta1 <- x_inter[,1]*reg3$residuals/varX_in[1,1]
  W <- t(eta1)%*%eta1/n0
  

  phi <- (x_inter/sd^2*result_inter$residuals)%*%solve(t(x_inter)%*%(x_inter/sd^2)/n0)
  hattheta2 <- hattheta1 - (t(phi)%*%eta1/n0)%*%solve(V/rho + W)%*%(hbeta1 - tildabeta1)

  
  htheta[i,] <- c(as.numeric(hattheta1), as.numeric(hattheta2))
}

MSE <- colMeans(sweep(htheta, 2, rep(1,4))^2)
RMSE <- sqrt(MSE)

path.save <- paste0('./Results/')
dir.create(path.save, showWarnings=FALSE)
path.image <- paste0(path.save,'/heteroSimu.RData')
save.image(file=path.image)

path.graph <- paste0(path.save,'/heterogeneous.pdf')
pdf(file=path.graph, width=9, height=9)
op <- par(mfcol=c(1, 1), mai=c(0.6,0.5,0.6,0.4), 
          oma=c(2,2,2,2), cex=1)
local.box <- c(0.5,1.5,3.5,4.5)
# col.box <- c('white', 'grey', 'white', 'grey')
boxplot(cbind(htheta[,c(1,3)],htheta[,c(2,4)]), at=local.box, outline=FALSE, xaxt='n', 
        lwd=2, cex.axis=2.5, boxwex=0.5
        # , ylim=c(1.8,2.2)
        )
title("",cex.main=3, line=2)
abline(h=1)
axis(1, at = c(0.5,1.5, 3.5,4.5), labels = c("INT","EFF", "INT", "EFF"), line=1, tick=FALSE, cex.axis=2.5)
mtext(side=1, at=c(1,4), text=c(expression(tau[1]),expression(tau[2])), cex=2.5, line = 4)
par(op)
dev.off()

