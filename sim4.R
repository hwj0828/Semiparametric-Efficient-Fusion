# Simulation for Scenario I with Neural Networks and Chernozhukov's Sample Splitting
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

library(parallel)
library(boot)
library(MASS)
library(torch)

# Configure PyTorch for multiprocessing compatibility
torch_set_num_threads(1)  # Prevent thread conflicts
if (torch_is_installed()) {
  # Set sharing strategy for multiprocessing
  options(torch.serialization_version = 2)
}

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

# Define neural network architecture
propensity_net <- nn_module(
  "PropensityNet",
  initialize = function() {
    self$layer1 <- nn_linear(1, 8)
    self$layer2 <- nn_linear(8, 8)
    self$layer3 <- nn_linear(8, 8)
    self$output <- nn_linear(8, 1)
    self$sigmoid <- nn_sigmoid()
    self$relu <- nn_relu()
  },
  forward = function(x) {
    x <- self$layer1(x)
    x <- self$relu(x)
    x <- self$layer2(x)
    x <- self$relu(x)
    x <- self$layer3(x)
    x <- self$relu(x)
    x <- self$output(x)
    x <- self$sigmoid(x)
    x
  }
)

outcome_net <- nn_module(
  "OutcomeNet",
  initialize = function() {
    self$layer1 <- nn_linear(1, 8)
    self$layer2 <- nn_linear(8, 8)
    self$layer3 <- nn_linear(8, 8)
    self$output <- nn_linear(8, 1)
    self$relu <- nn_relu()
  },
  forward = function(x) {
    x <- self$layer1(x)
    x <- self$relu(x)
    x <- self$layer2(x)
    x <- self$relu(x)
    x <- self$layer3(x)
    x <- self$relu(x)
    x <- self$output(x)
    x
  }
)

esti.INT <- function(x, y, A, n0) {
  # Check if we have both treated and control units
  if(length(unique(A)) != 2) {
    return(c(NA, NA))
  }
  
  K <- 4
  folds <- sample(rep(1:K, length.out=n0))
  
  psi_vals <- numeric(n0)  # store influence function contributions
  
  for (k in 1:K) {
    # define train/test
    train_idx <- which(folds != k)
    test_idx  <- which(folds == k)
    
    x_train <- x[train_idx,]
    x_test  <- x[test_idx,]
    A_train <- A[train_idx]
    A_test  <- A[test_idx]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # Convert to tensors
    x_tensor_train <- torch_tensor(matrix(as.numeric(x_train), ncol = 1), dtype = torch_float())
    y_tensor_train <- torch_tensor(as.numeric(y_train), dtype = torch_float())
    A_tensor_train <- torch_tensor(as.numeric(A_train), dtype = torch_float())
    
    # Train propensity score model
    ps_model <- propensity_net()
    optimizer_ps <- optim_adam(ps_model$parameters, lr = 0.002)
    
    for(epoch in 1:400) {
      optimizer_ps$zero_grad()
      ps_pred_train <- ps_model(x_tensor_train)$squeeze()
      loss <- nnf_binary_cross_entropy(ps_pred_train, A_tensor_train)
      loss$backward()
      optimizer_ps$step()
    }
    
    # Train outcome models
    treated_idx_train <- A_train == 1
    control_idx_train <- A_train == 0
    
    # Treated group model
    model_mu1 <- outcome_net()
    optimizer_mu1 <- optim_adam(model_mu1$parameters, lr = 0.002)
    
    if (sum(treated_idx_train) > 0) {
      x_treated <- x_tensor_train[treated_idx_train, ]
      y_treated <- y_tensor_train[treated_idx_train]
      
      for(epoch in 1:400) {
        optimizer_mu1$zero_grad()
        mu1_pred_train <- model_mu1(x_treated)$squeeze()
        loss <- nnf_mse_loss(mu1_pred_train, y_treated)
        loss$backward()
        optimizer_mu1$step()
      }
    }
    
    # Control group model
    model_mu0 <- outcome_net()
    optimizer_mu0 <- optim_adam(model_mu0$parameters, lr = 0.002)
    
    if (sum(control_idx_train) > 0) {
      x_control <- x_tensor_train[control_idx_train, ]
      y_control <- y_tensor_train[control_idx_train]
      
      for(epoch in 1:400) {
        optimizer_mu0$zero_grad()
        mu0_pred_train <- model_mu0(x_control)$squeeze()
        loss <- nnf_mse_loss(mu0_pred_train, y_control)
        loss$backward()
        optimizer_mu0$step()
      }
    }
    
    x_tensor_test <- torch_tensor(matrix(as.numeric(x_test), ncol = 1), dtype = torch_float())
    
    with_no_grad({
      ps_test <- as.numeric(ps_model(x_tensor_test)$squeeze())
      mu1_test <- as.numeric(model_mu1(x_tensor_test)$squeeze())
      mu0_test <- as.numeric(model_mu0(x_tensor_test)$squeeze())
    })
    
    ps_test <- pmin(pmax(ps_test, 0.001), 0.999)
    
    psi_vals[test_idx] <- A_test * y_test / ps_test - (1 - A_test) * y_test / (1 - ps_test) -
      (A_test / ps_test - 1) * mu1_test + ((1 - A_test) / (1 - ps_test) - 1) * mu0_test
  }
  
  # Average the K folds
  htau_INT <- mean(psi_vals, na.rm = TRUE)
  variance <- mean((psi_vals - htau_INT)^2, na.rm = TRUE) / n0
  stderror <- sqrt(variance)
  
  return(c(htau_INT, stderror))
}

esti.PRM <- function(x, y, A, n0, n1, hbeta1) {
  n <- n0 + n1
  
  K <- 4
  folds <- sample(rep(1:K, length.out=n0))
  
  psi_vals <- numeric(n0)  # store influence function contributions
  
  for (k in 1:K) {
    # define train/test
    train_idx <- which(folds != k)
    test_idx  <- which(folds == k)
    
    x_train <- x[train_idx,]
    x_test  <- x[test_idx,]
    A_train <- A[train_idx]
    A_test  <- A[test_idx]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # Convert to tensors
    x_tensor_train <- torch_tensor(matrix(as.numeric(x_train), ncol = 1), dtype = torch_float())
    y_tensor_train <- torch_tensor(as.numeric(y_train), dtype = torch_float())
    A_tensor_train <- torch_tensor(as.numeric(A_train), dtype = torch_float())
    
    # Train propensity score model
    ps_model <- propensity_net()
    optimizer_ps <- optim_adam(ps_model$parameters, lr = 0.002)
    
    for(epoch in 1:400) {
      optimizer_ps$zero_grad()
      ps_pred_train <- ps_model(x_tensor_train)$squeeze()
      loss <- nnf_binary_cross_entropy(ps_pred_train, A_tensor_train)
      loss$backward()
      optimizer_ps$step()
    }
    
    # Train outcome models
    treated_idx_train <- A_train == 1
    control_idx_train <- A_train == 0
    
    # Treated group model
    model_mu1 <- outcome_net()
    optimizer_mu1 <- optim_adam(model_mu1$parameters, lr = 0.002)
    
    if (sum(treated_idx_train) > 0) {
      x_treated <- x_tensor_train[treated_idx_train, ]
      y_treated <- y_tensor_train[treated_idx_train]
      
      for(epoch in 1:400) {
        optimizer_mu1$zero_grad()
        mu1_pred_train <- model_mu1(x_treated)$squeeze()
        loss <- nnf_mse_loss(mu1_pred_train, y_treated)
        loss$backward()
        optimizer_mu1$step()
      }
    }
    
    # Control group model
    model_mu0 <- outcome_net()
    optimizer_mu0 <- optim_adam(model_mu0$parameters, lr = 0.002)
    
    if (sum(control_idx_train) > 0) {
      x_control <- x_tensor_train[control_idx_train, ]
      y_control <- y_tensor_train[control_idx_train]
      
      for(epoch in 1:400) {
        optimizer_mu0$zero_grad()
        mu0_pred_train <- model_mu0(x_control)$squeeze()
        loss <- nnf_mse_loss(mu0_pred_train, y_control)
        loss$backward()
        optimizer_mu0$step()
      }
    }
    
    x_tensor_test <- torch_tensor(matrix(as.numeric(x_test), ncol = 1), dtype = torch_float())
    
    with_no_grad({
      ps_test <- as.numeric(ps_model(x_tensor_test)$squeeze())
      mu1_test <- as.numeric(model_mu1(x_tensor_test)$squeeze())
      mu0_test <- as.numeric(model_mu0(x_tensor_test)$squeeze())
    })
    
    ps_test <- pmin(pmax(ps_test, 0.001), 0.999)
    
    psi_vals[test_idx] <- A_test * y_test / ps_test - (1 - A_test) * y_test / (1 - ps_test) -
      (A_test / ps_test - 1) * mu1_test + ((1 - A_test) / (1 - ps_test) - 1) * mu0_test
  }
  
  
  # Average the K folds
  design <- cbind(rep(1,n0), x, A )
  residue <- c(y - design%*%hbeta1)
  htau_PRM <- mean(psi_vals - colMeans(c(psi_vals*residue)*design)%*%
                     solve(t(design)%*%diag(residue^2)%*%design/n0)%*%t(residue*design))
  variance <- (mean((psi_vals - htau_PRM)^2) + (n0/n1 - 1)*colMeans(c(psi_vals*residue)*design)%*%
                 solve(t(design)%*%diag(residue^2)%*%design/n0)%*%colMeans(c(psi_vals*residue)*design))/n0
  stderror <- sqrt(variance)
  
  return(c(htau_PRM, stderror))
}

esti.EFF <- function(x,y,A,hbeta1,n0,n1){
  n <- n0+n1
  
  K <- 4
  folds <- sample(rep(1:K, length.out=n0))
  
  psi_vals <- numeric(n0)  # store influence function contributions
  
  for (k in 1:K) {
    # define train/test
    train_idx <- which(folds != k)
    test_idx  <- which(folds == k)
    
    x_train <- x[train_idx,]
    x_test  <- x[test_idx,]
    A_train <- A[train_idx]
    A_test  <- A[test_idx]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # Convert to tensors
    x_tensor_train <- torch_tensor(matrix(as.numeric(x_train), ncol = 1), dtype = torch_float())
    y_tensor_train <- torch_tensor(as.numeric(y_train), dtype = torch_float())
    A_tensor_train <- torch_tensor(as.numeric(A_train), dtype = torch_float())
    
    # Train propensity score model
    ps_model <- propensity_net()
    optimizer_ps <- optim_adam(ps_model$parameters, lr = 0.002)
    
    for(epoch in 1:400) {
      optimizer_ps$zero_grad()
      ps_pred_train <- ps_model(x_tensor_train)$squeeze()
      loss <- nnf_binary_cross_entropy(ps_pred_train, A_tensor_train)
      loss$backward()
      optimizer_ps$step()
    }
    
    # Train outcome models
    treated_idx_train <- A_train == 1
    control_idx_train <- A_train == 0
    
    # Treated group model
    model_mu1 <- outcome_net()
    optimizer_mu1 <- optim_adam(model_mu1$parameters, lr = 0.002)
    
    if (sum(treated_idx_train) > 0) {
      x_treated <- x_tensor_train[treated_idx_train, ]
      y_treated <- y_tensor_train[treated_idx_train]
      
      for(epoch in 1:400) {
        optimizer_mu1$zero_grad()
        mu1_pred_train <- model_mu1(x_treated)$squeeze()
        loss <- nnf_mse_loss(mu1_pred_train, y_treated)
        loss$backward()
        optimizer_mu1$step()
      }
    }
    
    # Control group model
    model_mu0 <- outcome_net()
    optimizer_mu0 <- optim_adam(model_mu0$parameters, lr = 0.002)
    
    if (sum(control_idx_train) > 0) {
      x_control <- x_tensor_train[control_idx_train, ]
      y_control <- y_tensor_train[control_idx_train]
      
      for(epoch in 1:400) {
        optimizer_mu0$zero_grad()
        mu0_pred_train <- model_mu0(x_control)$squeeze()
        loss <- nnf_mse_loss(mu0_pred_train, y_control)
        loss$backward()
        optimizer_mu0$step()
      }
    }
    
    x_tensor_test <- torch_tensor(matrix(as.numeric(x_test), ncol = 1), dtype = torch_float())
    
    with_no_grad({
      ps_test <- as.numeric(ps_model(x_tensor_test)$squeeze())
      mu1_test <- as.numeric(model_mu1(x_tensor_test)$squeeze())
      mu0_test <- as.numeric(model_mu0(x_tensor_test)$squeeze())
    })
    
    ps_test <- pmin(pmax(ps_test, 0.001), 0.999)
    
    psi_vals[test_idx] <- A_test * y_test / ps_test - (1 - A_test) * y_test / (1 - ps_test) -
      (A_test / ps_test - 1) * mu1_test + ((1 - A_test) / (1 - ps_test) - 1) * mu0_test
  }
  
  # Average the K folds
  lm_inter <- lm(y~x+A) 
  hbeta2 <- lm_inter$coef
  
  design <- cbind(rep(1,n0), x, A )
  residue <- c(y - design%*%hbeta2)
  htau_EFF <- mean(psi_vals) - n1/n*colMeans(c(psi_vals*residue)*design)%*%
    solve(t(design)%*%diag(residue^2)%*%design/n0)%*%(t(design)%*%design/n0)%*%(hbeta2 - hbeta1)
  variance <- (mean((psi_vals - c(htau_EFF))^2) - n1/n*colMeans(c(psi_vals*residue)*design)%*%
                 solve(t(design)%*%diag(residue^2)%*%design/n0)%*%colMeans(c(psi_vals*residue)*design))/n0
  stderror <- sqrt(variance)
  
  return(c(htau_EFF, stderror))
}

esti.KNW <- function(x,y,A,beta,n0,n1){
  n <- n0+n1
  
  K <- 4
  folds <- sample(rep(1:K, length.out=n0))
  
  psi_vals <- numeric(n0)  # store influence function contributions
  
  for (k in 1:K) {
    # define train/test
    train_idx <- which(folds != k)
    test_idx  <- which(folds == k)
    
    x_train <- x[train_idx,]
    x_test  <- x[test_idx,]
    A_train <- A[train_idx]
    A_test  <- A[test_idx]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # Convert to tensors
    x_tensor_train <- torch_tensor(matrix(as.numeric(x_train), ncol = 1), dtype = torch_float())
    y_tensor_train <- torch_tensor(as.numeric(y_train), dtype = torch_float())
    A_tensor_train <- torch_tensor(as.numeric(A_train), dtype = torch_float())
    
    # Train propensity score model
    ps_model <- propensity_net()
    optimizer_ps <- optim_adam(ps_model$parameters, lr = 0.002)
    
    for(epoch in 1:400) {
      optimizer_ps$zero_grad()
      ps_pred_train <- ps_model(x_tensor_train)$squeeze()
      loss <- nnf_binary_cross_entropy(ps_pred_train, A_tensor_train)
      loss$backward()
      optimizer_ps$step()
    }
    
    # Train outcome models
    treated_idx_train <- A_train == 1
    control_idx_train <- A_train == 0
    
    # Treated group model
    model_mu1 <- outcome_net()
    optimizer_mu1 <- optim_adam(model_mu1$parameters, lr = 0.002)
    
    if (sum(treated_idx_train) > 0) {
      x_treated <- x_tensor_train[treated_idx_train, ]
      y_treated <- y_tensor_train[treated_idx_train]
      
      for(epoch in 1:400) {
        optimizer_mu1$zero_grad()
        mu1_pred_train <- model_mu1(x_treated)$squeeze()
        loss <- nnf_mse_loss(mu1_pred_train, y_treated)
        loss$backward()
        optimizer_mu1$step()
      }
    }
    
    # Control group model
    model_mu0 <- outcome_net()
    optimizer_mu0 <- optim_adam(model_mu0$parameters, lr = 0.002)
    
    if (sum(control_idx_train) > 0) {
      x_control <- x_tensor_train[control_idx_train, ]
      y_control <- y_tensor_train[control_idx_train]
      
      for(epoch in 1:400) {
        optimizer_mu0$zero_grad()
        mu0_pred_train <- model_mu0(x_control)$squeeze()
        loss <- nnf_mse_loss(mu0_pred_train, y_control)
        loss$backward()
        optimizer_mu0$step()
      }
    }
    
    x_tensor_test <- torch_tensor(matrix(as.numeric(x_test), ncol = 1), dtype = torch_float())
    
    with_no_grad({
      ps_test <- as.numeric(ps_model(x_tensor_test)$squeeze())
      mu1_test <- as.numeric(model_mu1(x_tensor_test)$squeeze())
      mu0_test <- as.numeric(model_mu0(x_tensor_test)$squeeze())
    })
    
    ps_test <- pmin(pmax(ps_test, 0.001), 0.999)
    
    psi_vals[test_idx] <- A_test * y_test / ps_test - (1 - A_test) * y_test / (1 - ps_test) -
      (A_test / ps_test - 1) * mu1_test + ((1 - A_test) / (1 - ps_test) - 1) * mu0_test
  }
  
  # Average the K folds
  design <- cbind(rep(1,n0), x, A )
  residue <- c(y - design%*%beta)
  htau_KNW <- mean(psi_vals - colMeans(c(psi_vals*residue)*design)%*%
                     solve(t(design)%*%diag(residue^2)%*%design/n0)%*%t(residue*design))
  variance <- (mean((psi_vals - htau_KNW)^2) - colMeans(c(psi_vals*residue)*design)%*%
                 solve(t(design)%*%diag(residue^2)%*%design/n0)%*%colMeans(c(psi_vals*residue)*design))/n0
  stderror <- sqrt(variance)
  
  c(htau_KNW, stderror)
}

# Simulation function with cross-fitting
simuOne <- function(n, n0, n1) {
  alpha <- c(1,1,1)
  tau <- 1
  beta <- c(1.2306, 0.6069, 0.5929)
  
  x <- rnorm(n, mean=0, sd=1)
  A <- rbinom(n,1, 1/(1+exp(-(1-x))))
  y1 <- cbind(1,x,x^2)%*%alpha + rnorm(n, mean=0, sd=2)
  y0 <- cbind(1,x,x^2)%*%(alpha - c(0,0,1)) + rnorm(n, mean=0, sd=1)
  y <- A*y1+(1-A)*y0
  
  x_inter <- matrix(x[1:n0], ncol = 1)
  A_inter <- A[1:n0]
  y_inter <- y[1:n0]
  
  x_exter <- x[(n0+1):(n0+n1)]
  A_exter <- A[(n0+1):(n0+n1)]
  y_exter <- y[(n0+1):(n0+n1)]
  
  lm1 <- lm(y_exter~x_exter+A_exter)
  hbeta1 <- lm1$coef
  
  ## efficient estimator using only internal data
  result_INT <- esti.INT(x_inter, y_inter, A_inter, n0)
  if(any(is.na(result_INT))) {
    cat("Warning: esti.INT returned NA values\n")
    return(list(htau_INT=NA, htau_PRM=NA, htau_EFF=NA, htau_KNW=NA,
                stderror_INT=NA, stderror_PRM=NA, stderror_EFF=NA, stderror_KNW=NA,
                cv_INT=NA, cv_PRM=NA, cv_EFF=NA, cv_KNW=NA))
  }
  htau_INT <- result_INT[1]
  stderror_INT <- result_INT[2]
  cv_INT <- (tau >= htau_INT - 1.96*stderror_INT) & 
    (tau <= htau_INT + 1.96*stderror_INT)
  
  ## the crude estimator ignoring uncertainty of $\tilde{\beta}$
  result_PRM <- esti.PRM(x_inter, y_inter, A_inter, n0, n1, hbeta1)
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
  
  # list(htau_INT=htau_INT,stderror_INT=stderror_INT, cv_INT=cv_INT)
  
  list(htau_INT=htau_INT, htau_PRM=htau_PRM, htau_EFF=htau_EFF, htau_KNW=htau_KNW,
       stderror_INT=stderror_INT, stderror_PRM=stderror_PRM, stderror_EFF=stderror_EFF, stderror_KNW=stderror_KNW,
       cv_INT=cv_INT, cv_PRM=cv_PRM, cv_EFF=cv_EFF, cv_KNW=cv_KNW)
}

# Simulation parameters
nsims <- 1000
ncpus <- 1
n <- 3000 # = n0 + n1
n0 <- 1000 #internal sample size in the maintext
n1 <- 2000 #external sample size in the main text

duration <- Sys.time()
duration

# set RNG seed for reproducibility 
RNGkind("L'Ecuyer-CMRG")
set.seed(100)

# # Test single simulation first
# test_result <- simuOne(n, n0, n1)
# print("Single simulation test:")
# print(test_result)

# Run simulations
simu <- t(parReplicate(nsims, expr=simuOne(n, n0, n1), 
                       simplify=TRUE, mc.cores=ncpus, mc.set.seed=TRUE))

# Process results
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

# Save results
path.save <- paste0('./Results/')
dir.create(path.save, recursive=TRUE)
path.image <- paste0(path.save,'/',n1,'Simu7.RData')  # Changed to Simu5
path.out <- paste0(path.save,'/simu7.out')  # Changed to simu5
save.image(file=path.image)

# Create plot
path.graph <- paste0(path.save,'/m=',n1,'_NN_cross-fitting.pdf')  # Updated title for Chernozhukov's approach
pdf(file=path.graph, width=9, height=9)
op <- par(mfcol=c(1, 1), mai=c(0.6,0.5,0.6,0.4), 
          oma=c(2,2,2,2), cex=1)
ylim.box=c(0.4, 1.6)
boxplot(htau, outline=FALSE, xaxt='n', lwd=2, cex.axis=2.5,
        ylim=ylim.box, boxwex = 0.5)
title(main=bquote("m ="~ .(n1)~"(Neural network)"),cex.main=3, line=2)  # Updated title
axis(1, at = 1:4, labels = c("INT", "PRM", "EFF", "KNW"), line=1.5, tick=FALSE, cex.axis=2.5)
abline(h=1)
par(op)
dev.off() 
