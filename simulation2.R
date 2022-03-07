rm(list = ls())
library(MASS) # Import function mvrnorm
library(rTensor) # Tensor library
library(pracma) # Import function sqrtm
library(PMA) # PMD method
source("GLAA_SVD.R") # GLAA algorithm
source("models.R") # model settings
source("utility.R") # auxiliary functions

# ---------------------- Scenario 2 --------------------#
# This code reproduces one replicate with p1 = p2 = 100, and p3 = 20. The ranks for the three modes are r1 = r2 = 2, r3 = 1. The sparsity levels are s1 = s2 = s3 = 5. And the sample size n = 500. Please refer to Section 5.1 for more detailed information. 

times <- 1 # The number of replicates

## Set random seed
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

model <- Model2(pp = 20) # Set p3 = 20.
p <- model$p # The dimension p1, p2, and p3.
r <- model$r # The ranks r1, r2, and r3.
s <- model$s # The sparsity level s1, s2, and s3.
n <- model$n # The sample size.
sparse_mode <- model$sparse_mode # Specify the sparse modes.
a1.list <- model$a1.list # Tuning parameter for the first two modes.
a2.list <- model$a2.list # Tuning parameter for the third mode.
Gamma <- model$Gamma # The basis matrices Gamma1, Gamma2, Gamma3.
data.gen <- model$data.gen # The function generating the simulated data in Scenario 2.

output <- sapply(seq_len(times), function(i){
  cat("Time", i, '\n')
  # ------------------------ Data generation ------------------------ #
  data <- data.gen(n) # Generate three training data sets
  x <- data$x
  y <- data$y
  z <- data$z
  # standardization
  x <- scale(x)
  y <- scale(y)
  z <- scale(z)

  data.test <- data.gen(n) # Generate three test data sets
  x.test <- data.test$x
  y.test <- data.test$y
  z.test <- data.test$z
  # standardization
  x.test <- scale(x.test)
  y.test <- scale(y.test)
  z.test <- scale(z.test)

  # Compute the sample estimator Delta tilde for the test data set
  Delta.test <- tensor_mean(x.test, y.test, z.test)
  Delta.test <- as.tensor(Delta.test)
  
  # ----------------------- GLAA ------------------------- #
  a1.list.cycle <- rep(a1.list, each = length(a2.list))
  a2.list.cycle <- rep(a2.list, length(a1.list))
  fit.results <- mapply(function(a1, a2){ # Implement Algorithm 1 for each tuning parameter pair (a1, a2).
    fit <- STATSVD(x, y, z, r, s, sparse_mode = sparse_mode, tmax = 50, a1 = a1, a2 = a2) # Sparse tensor decomposition algorithm (Algorithm 1)
    dist.true <- sapply(1:length(p), function(i){
      subspace(fit[[i]], Gamma[[i]])
    })
    true.dist <- mean(dist.true) # The average subspace distance for all three modes.
    s.list <- lapply(1:length(p), function(i){
      which(apply(fit[[i]], 1, function(x){any(x != 0)})) # The estimated active set for each mode.
    })
    
    # ----- validation ------ #
    proj <- list(fit[[1]] %*% t(fit[[1]]), fit[[2]] %*% t(fit[[2]]), fit[[3]] %*% t(fit[[3]]))
    Delta.test.2 <- ttl(Delta.test, proj, ms = c(1,2,3)) # The projection of Delta tilde from the test data set onto the subspaces spanned by the estimated basis matrices.
    error <- rTensor::fnorm(Delta.test - Delta.test.2) # The error defined in (5)
    # ----------------------- #
    list(true.dist = true.dist, s.list = s.list, error = error)
  }, a1.list.cycle, a2.list.cycle)
  
  true.dist <- matrix(do.call(cbind, fit.results[1,]), length(a1.list), length(a2.list), byrow = TRUE)
  s.list <- matrix(fit.results[2,], length(a1.list), length(a2.list), byrow = TRUE)
  error <- matrix(do.call(cbind, fit.results[3,]), length(a1.list), length(a2.list), byrow = TRUE)
  
  ind <- which(error == min(error), arr.ind = TRUE) # Select the optimal tuning parameter.
  ind.i <- ind[1,1]
  ind.j <- ind[1,2]
  dist <- true.dist[ind.i, ind.j] # The subspace estimation error corresponding to the optimal tuning parameter.
  s.list.final <- s.list[[ind.i, ind.j]] # The selected active sets corresponding to the optimal tuning parameter.
  
  cat(paste0("The ", "(", ind.i, ", ", ind.j, ")", "-th parameter is the optimal one. \n"))
  
  TFPR.list <- sapply(seq_len(3), function(i){
    if(p[i]==s[i]){c(1,0)}
    else{
      TPR <- sum(s.list.final[[i]] %in% 1:s[i])/s[i] # The True Positive Rate on each mode
      FPR <- sum(s.list.final[[i]] %in% (s[i]+1):p[i])/(p[i] - s[i]) # The False Positive Rate on each mode
      c(TPR, FPR)
    }
  })
  cat("Set X: ", paste(s.list.final[[1]], collapse = " "), " | TPR(X):", TFPR.list[1,1], " | FPR(X):", TFPR.list[2,1], "\n",
      "Set Y: ", paste(s.list.final[[2]], collapse = " "), " | TPR(Y):", TFPR.list[1,2], " | FPR(Y):", TFPR.list[2,2], "\n",
      "Set Z: ", paste(s.list.final[[3]], collapse = " "), " | TPR(Z):", TFPR.list[1,3], " | FPR(Z):", TFPR.list[2,3], "\n",
      "Mean dist: ", dist, "\n\n", sep = "")
  
  GLAA.result <- c(c(TFPR.list), dist)
  # ---------------------------------- #

  # ----------------  Univariate LA -------------------- #
  ula <- tensor_mean(x, y ,z)
  ula <- as.tensor(ula) # Compute the sample estimator Delta tilde

  var.ula.list <- ula.variable(ula, p, s) # Variable selection for ULA: for each mode-k matricization, select rows with the first s_k largest l_2 norm.
  TFPR.ula <- sapply(1:3, function(i){ # The TPR and FPR for ULA
    ind <- var.ula.list[[i]]
    if((p[i] == s[i])){c(1,0)}
    else{
      TPR <- sum(ind %in% 1:s[i])/s[i]
      FPR <- sum(ind %in% (s[i]+1):p[i])/(p[i] - s[i])
      c(TPR, FPR)
    }
  })
  Gamma.ula <- lapply(1:length(p), function(k){ # Estimate the basis matrices for ULA: the first r_k left singular vectors of each mode-k matricization.
    Gamma <- svd(k_unfold(ula, k)@data)$u[,1:r[k], drop = FALSE]
    Gamma
  })
  dist.ula.list <- sapply(1:length(p), function(i){
    subspace(Gamma.ula[[i]], Gamma[[i]])
  })
  dist.ula <- mean(dist.ula.list) # The average subspace distance for ULA.
  ula.result <- c(c(TFPR.ula), dist.ula)
  # ----------------------------------------------- #

  c(GLAA.result, ula.result)
})

output <- as.data.frame(t(output)) # Record the output
colnames(output) <- c( "GLAA_T1", "GLAA_F1", "GLAA_T2", "GLAA_F2", "GLAA_T3", "GLAA_F3", "GLAA_D", "ULA_T1", "ULA_F1", "ULA_T2", "ULA_F2", "ULA_T3", "ULA_F3", "ULA_D")
print(output)
# --------------------------------------------------------#
