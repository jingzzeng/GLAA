# Implements the sparse tensor decomposition algorithm (Algorithm 1).

######################################################
# STAT-SVD for real data
# ###########
# Input:
# x,y,z: input data.
# r: specified ranks. 
# tmax: maximal iteration.
# sparse_mode: specifies which modes to impose sparsity.
# vartol: the convergence tolerance.
# a1: the tuning parameter in the initialization step.
# a2: the tuning parameter in the iteration step.
# LA: the pre-specified Delta tilde tensor, defined in Section 3.
######################################################
STATSVD.real <- function(x, y, z, r, tmax = 10, sparse_mode, vartol = 1e-6, a1 = 1, a2 = 1, LA){
  if(missing(LA)){
    try(if(missing("x") | missing("y") | missing("z")) stop("missing argument: X, Y or Z."))
    try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
    if(is.vector(x)){x <- as.matrix(x)}
    if(is.vector(y)){y <- as.matrix(y)}
    if(is.vector(z)){z <- as.matrix(z)}
    if((dim(x)[1] != dim(y)[1]) | (dim(x)[1] != dim(z)[1]) | (dim(y)[1] != dim(z)[1])) {stop("unmatched dimension: sample sizes of X, Y and Z are supposed to be equal.")}
    
    if(dim(z)[2] == 1){
      p <- c(dim(x)[2], dim(y)[2])
      d <- 2
    }
    else{
      p <- c(dim(x)[2], dim(y)[2], dim(z)[2])
      d <- 3
    }
    x <- scale(x)
    y <- scale(y)
    z <- scale(z)
    Delta <- tensor_mean(x,y,z)
    Delta <- as.tensor(Delta)
  }
  else{
    Delta <- LA
    p <- dim(Delta)
    if(p[3] == 1){
      p <- p[1:2]
      d <- 2
    }
    else{
      p <- p
      d <- 3
    }
  }
  
  try(if(missing("sparse_mode")) sparse_mode<-rep(TRUE, d))
  
  p.prod <- prod(p)
  p.minus.k <- p.prod / p
  if(is.atomic(r) && length(r)==1){
    r <- rep(r, d)
  }
  r.minus.k <- prod(r) / r # Introduce r_{-k}
  
  try(if(d != length(r)) stop("invalid input: r and the order of Y is incompatible."))
  
  Gamma_list <- list(); I_0 <- list()
  for(k in 1:d){ # Initialization Step1: (1) find the active set I_k
    if(isTRUE(sparse_mode[k])){
      Delta_k <- k_unfold(Delta, k)@data
      Delta_k_row_max <- apply(abs(Delta_k), 1, max)
      thrd.k <- a1[k]
      I_0 <- c(I_0, list(Delta_k_row_max > thrd.k))
    }
    else{
      I_0 <- c(I_0, list(rep(TRUE,p[k])))
    }
  }
  
  Delta_sparse <- Delta
  for (k in 1:d){ # Initialization Step1: (2) construct Delta_sparse
    Delta_sparse <- ttm(Delta_sparse, diag(I_0[[k]]*1), k)
  }
  for (k in 1:d){ # Initialization Step1: (3) find the initial basis matrix
    if(isTRUE(sparse_mode[k])){
      datai <- k_unfold(Delta_sparse, k)@data # mode-k matricization of Delta_sparse
      selec <- (1:p[k])[I_0[[k]]]
      if(length(selec)<r[k]){  # In case that sk < rk
        newdata <- matrix(0, nrow=nrow(datai), ncol = ncol(datai))
        newdata[selec,] <- datai[selec,]
        Gamma <- (svd(newdata)$u[,1:r[k]])
      }
      else{ # Else, implement SVD on the selected rows first.
        Gamma <- matrix(0, nrow = p[k], ncol = r[k])
        Gamma[I_0[[k]],] <- svd(datai[I_0[[k]],,drop=FALSE])$u[,1:r[k]]
      }
      Gamma_list <- c(Gamma_list, list(t(Gamma)))
    }
    else{
      Gamma_list <- c(Gamma_list, list(t(svd(k_unfold(Delta_sparse, k)@data)$u[,1:r[k]])))
    }
  }
  
  t <- 1
  approx <- -1
  sv_k <- 0
  while(t<tmax){  # Step 2: iteratively update the basis matrices. Stop criterion: convergence or maximum number of iteration reached
    for(k in 1:d){ # For each mode k
      if(!isTRUE(sparse_mode[k])){ # For dense mode, skip Step 2a.
        A <- ttl(Delta, Gamma_list[-k], (1:d)[-k])
        A_k <- k_unfold(A, k)@data
        A_k.svd <- svd(A_k)
        Gamma_list[[k]] <- t(A_k.svd$u[,1:r[k]])
        sv_k <- A_k.svd$d[1:r[k]]
      }
      else{ # For sparse mode, implement Step 2a and 2b.
        A <- ttl(Delta, Gamma_list[-k], (1:d)[-k])
        A_k <- k_unfold(A, k)@data
        eta_k <- a2[k]
        I_k <- apply(A_k, 1, function(x){sum(x^2)}) > eta_k # Step 2a: Update the active set
        selec <- (1:p[k])[I_k]
        
        # If s_k < r_k, do not update
        if(length(selec) < r[k])
          next
        
        B_k <- A_k[I_k,,drop=FALSE] # Step 2b: Perform SVD
        ###
        Gamma <- matrix(0, nrow(A_k), r[k])
        B_k.svd <- svd(B_k)
        Gamma[I_k,] <- B_k.svd$u[,1:r[k]]
        Gamma_list[[k]] <- t(Gamma)
        sv_k <- B_k.svd$d[1:r[k]] # The first r_k singular values of B_k.
        ###
      }
    }
    if (abs(sum(sv_k^2) - approx) > vartol & t<tmax){
      t <- t+1
      approx <- sum(sv_k^2)
    }
    else {
      break
    }
  }

  for(k in 1:d){
    Gamma_list[[k]] <- t(Gamma_list[[k]])
  }
  return(Gamma_list)
}



######################################################
# STAT-SVD with two tuning parameters (for Scenario 2)
# ###########
# Input:
# Note: The input arguments are similar to those in the function STATSVD.real.
# s: Sparsity levels on each mode. Used in the theoretical forms for the tuning parameters.
# a1: Tuning parameter for the first two modes.
# a2: Tuning parameter for the third mode.
######################################################
STATSVD <- function(x, y, z, r, s, tmax = 10, sparse_mode, vartol = 1e-6, a1 = 1, a2 = 1){
  
  try(if(missing("x") | missing("y") | missing("z")) stop("missing argument: X, Y or Z."))
  try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
  if(is.vector(x)){x <- as.matrix(x)}
  if(is.vector(y)){y <- as.matrix(y)}
  if(is.vector(z)){z <- as.matrix(z)}
  if((dim(x)[1] != dim(y)[1]) | (dim(x)[1] != dim(z)[1]) | (dim(y)[1] != dim(z)[1])) {stop("unmatched dimension: sample sizes of X, Y and Z are supposed to be equal.")}
  
  if(dim(z)[2] == 1){
    p <- c(dim(x)[2], dim(y)[2])
    n <- dim(x)[1]
    d <- 2
  }
  else{
    p <- c(dim(x)[2], dim(y)[2], dim(z)[2])
    n <- dim(x)[1]
    d <- 3
  }
  
  # Sample estimator
  Delta <- tensor_mean(x,y,z)
  Delta <- as.tensor(Delta)
  
  try(if(missing("sparse_mode")) sparse_mode<-rep(TRUE, d))
  
  p.prod <- prod(p)
  p.minus.k <- p.prod / p
  s.prod <- prod(s)
  s.minus.k <- s.prod / s
  if(is.atomic(r) && length(r)==1){
    r <- rep(r, d)
  }
  r.minus.k <- prod(r) / r # Introduce r_{-k}
  
  try(if(d != length(r)) stop("invalid input: r and the order of Y is incompatible."))
  
  Gamma_list <- list(); I_0 <- list()
  for(k in 1:d){ # Initialization Step1: (1) find the active set I_k
    a <- ifelse(k == 3, a2, a1)
    if(isTRUE(sparse_mode[k])){
      Delta_k <- k_unfold(Delta, k)@data
      Delta_k_row_max <- apply(abs(Delta_k), 1, max)
      thrd.k <- sqrt(a*log(p.prod)/n)
      I_0 <- c(I_0, list(Delta_k_row_max > thrd.k))
    }
    else{
      I_0 <- c(I_0, list(rep(TRUE,p[k])))
    }
  }
  
  Delta_sparse <- Delta
  for (k in 1:d){ # Initialization Step1: (2) construct Delta_sparse
    Delta_sparse <- ttm(Delta_sparse, diag(I_0[[k]]*1), k)
  }
  for (k in 1:d){ # Initialization Step1: (3) find the initial basis matrix Gamma_0k
    if(isTRUE(sparse_mode[k])){
      datai <- k_unfold(Delta_sparse, k)@data # mode-k matricization of Delta_sparse
      selec <- (1:p[k])[I_0[[k]]]
      if(length(selec)<r[k]){  # in case that sk < rk
        newdata <- matrix(0, nrow=nrow(datai), ncol = ncol(datai))
        newdata[selec,] <- datai[selec,]
        Gamma <- (svd(newdata)$u[,1:r[k]])
      }
      else{ # Else, implement SVD on the selected rows first.
        Gamma <- matrix(0, nrow = p[k], ncol = r[k])
        Gamma[I_0[[k]],] <- svd(datai[I_0[[k]],,drop=FALSE])$u[,1:r[k]]
      }
      Gamma_list <- c(Gamma_list, list(t(Gamma)))
    }
    else{
      Gamma_list <- c(Gamma_list, list(t(svd(k_unfold(Delta_sparse, k)@data)$u[,1:r[k]])))
    }
  }
  
  t <- 1
  approx <- -1
  sv_k <- 0
  while(t<tmax){  # Step 2: iteratively update the basis matrices. Stop criterion: convergence or maximum number of iteration reached
    for(k in 1:d){ # For each mode k
      a <- ifelse(k == 3, a2, a1)
      
      if(!isTRUE(sparse_mode[k])){ # For dense mode, skip Step 2a.
        A <- ttl(Delta, Gamma_list[-k], (1:d)[-k])
        A_k <- k_unfold(A, k)@data
        A_k.svd <- svd(A_k)
        Gamma_list[[k]] <- t(A_k.svd$u[,1:r[k]])
        sv_k <- A_k.svd$d[1:r[k]]
      }
      else{ # For sparse mode, implement Step 2a and 2b.
        A <- ttl(Delta, Gamma_list[-k], (1:d)[-k])
        A_k <- k_unfold(A, k)@data
        eta_k <- (a*s.minus.k[k]*log(p.prod))/n
        I_k <- apply(A_k, 1, function(x){sum(x^2)}) > eta_k # Step 2a: Update the active set
        selec <- (1:p[k])[I_k]
        
        # If s_k < r_k, do not update
        if(length(selec) < r[k])
          next
        
        B_k <- A_k[I_k,,drop=FALSE] # Step 2b: Perform SVD
        ###
        Gamma <- matrix(0, nrow(A_k), r[k])
        B_k.svd <- svd(B_k)
        Gamma[I_k,] <- B_k.svd$u[,1:r[k]]
        Gamma_list[[k]] <- t(Gamma)
        sv_k <- B_k.svd$d[1:r[k]] # The first r_k singular values of B_k.
        ###
      }
    }
    if (abs(sum(sv_k^2) - approx) > vartol & t<tmax){
      t <- t+1
      approx <- sum(sv_k^2)
    }
    else {
      break
    }
  }
  
  for(k in 1:d){
    Gamma_list[[k]] <- t(Gamma_list[[k]])
  }
  return(Gamma_list)
}


######################################################
# STAT-SVD with one tuning parameters (for Scenarios 1,3)
# ###########
# Input:
# Note: The input arguments are similar to those in the function STATSVD.real.
# s: Sparsity levels on each mode. Used in the theoretical forms for the tuning parameters.
# a: Tuning parameter for the first two modes.
######################################################
STATSVD.one <- function(x, y, z, r, s, tmax = 10, sparse_mode, vartol = 1e-6, a = 1){
  
  try(if(missing("x") | missing("y") | missing("z")) stop("missing argument: X, Y or Z."))
  try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
  if(is.vector(x)){x <- as.matrix(x)}
  if(is.vector(y)){y <- as.matrix(y)}
  if(is.vector(z)){z <- as.matrix(z)}
  if((dim(x)[1] != dim(y)[1]) | (dim(x)[1] != dim(z)[1]) | (dim(y)[1] != dim(z)[1])) {stop("unmatched dimension: sample sizes of X, Y and Z are supposed to be equal.")}
  
  if(dim(z)[2] == 1){
    p <- c(dim(x)[2], dim(y)[2])
    n <- dim(x)[1]
    d <- 2
  }
  else{
    p <- c(dim(x)[2], dim(y)[2], dim(z)[2])
    n <- dim(x)[1]
    d <- 3
  }
  
  # Sample estimator
  Delta <- tensor_mean(x,y,z)
  Delta <- as.tensor(Delta)
  
  try(if(missing("sparse_mode")) sparse_mode<-rep(TRUE, d))
  
  p.prod = prod(p)
  p.minus.k <- p.prod / p
  s.prod <- prod(s)
  s.minus.k <- s.prod / s
  if(is.atomic(r) && length(r)==1){
    r <- rep(r, d)
  }
  r.minus.k <- prod(r) / r # Introduce r_{-k}
  
  try(if(d != length(r)) stop("invalid input: r and the order of Y is incompatible."))
  
  Gamma_list <- list(); I_0 <- list()
  for(k in 1:d){ # Initialization Step1: (1) find the active set I_k
    if(isTRUE(sparse_mode[k])){
      Delta_k <- k_unfold(Delta, k)@data
      Delta_k_row_max <- apply(abs(Delta_k), 1, max)
      thrd.k <- sqrt(a*log(p.prod)/n)
      I_0 <- c(I_0, list(Delta_k_row_max > thrd.k))
    }
    else{
      I_0 <- c(I_0, list(rep(TRUE,p[k])))
    }
  }
  
  Delta_sparse <- Delta
  for (k in 1:d){ # Initialization Step1: (2) construct Delta_sparse
    Delta_sparse <- ttm(Delta_sparse, diag(I_0[[k]]*1), k)
  }
  for (k in 1:d){ # Initialization Step1: (3) find the initial basis matrix Gamma_0k
    if(isTRUE(sparse_mode[k])){
      datai <- k_unfold(Delta_sparse, k)@data # mode-k matricization of Delta_sparse
      selec <- (1:p[k])[I_0[[k]]]
      if(length(selec)<r[k]){  # in case that sk < rk
        newdata <- matrix(0, nrow=nrow(datai), ncol = ncol(datai))
        newdata[selec,] <- datai[selec,]
        Gamma <- (svd(newdata)$u[,1:r[k]])
      }
      else{ # Else, implement SVD on the selected rows first.
        Gamma <- matrix(0, nrow = p[k], ncol = r[k])
        Gamma[I_0[[k]],] <- svd(datai[I_0[[k]],,drop=FALSE])$u[,1:r[k]]
      }
      Gamma_list <- c(Gamma_list, list(t(Gamma)))
    }
    else{
      Gamma_list <- c(Gamma_list, list(t(svd(k_unfold(Delta_sparse, k)@data)$u[,1:r[k]])))
    }
  }
  
  t <- 1
  approx <- -1
  sv_k <- 0
  while(t<tmax){  # Step 2: iteratively update the basis matrices. Stop criterion: convergence or maximum number of iteration reached
    for(k in 1:d){ # For each mode k
      if(!isTRUE(sparse_mode[k])){ # For dense mode, skip Step 2a.
        A <- ttl(Delta, Gamma_list[-k], (1:d)[-k])
        A_k <- k_unfold(A, k)@data
        A_k.svd <- svd(A_k)
        Gamma_list[[k]] <- t(A_k.svd$u[,1:r[k]])
        sv_k <- A_k.svd$d[1:r[k]]
      }
      else{ # For sparse mode, implement Step 2a and 2b.
        A <- ttl(Delta, Gamma_list[-k], (1:d)[-k])
        A_k <- k_unfold(A, k)@data
        eta_k <- (a*s.minus.k[k]*log(p.prod))/n
        I_k <- apply(A_k, 1, function(x){sum(x^2)}) > eta_k
        selec <- (1:p[k])[I_k]
        
        # If s_k < r_k, do not update
        if(length(selec) < r[k])
          next
        
        B_k <- A_k[I_k,,drop=FALSE] # Step 2b: Perform SVD
        ###
        Gamma <- matrix(0, nrow(A_k), r[k])
        B_k.svd <- svd(B_k)
        Gamma[I_k,] <- B_k.svd$u[,1:r[k]]
        Gamma_list[[k]] <- t(Gamma)
        sv_k <- B_k.svd$d[1:r[k]] # The first r_k singular values of B_k.
        ###
      }
    }
    if (abs(sum(sv_k^2) - approx) > vartol & t<tmax){
      t <- t+1
      approx <- sum(sv_k^2)
    }
    else {
      break
    }
  }
  
  for(k in 1:d){
    Gamma_list[[k]] <- t(Gamma_list[[k]])
  }
  return(Gamma_list)
}

######################################################
# Save algorithmic error over iterations. The function is based on function STATSVD, please refer to the comments in STATSVD for more details.
######################################################
STATSVD.algo.error <- function(x, y, z, r, s, tmax = 10, sparse_mode, vartol = 1e-6, a1 = 1, a2 = 1){
  
  try(if(missing("x") | missing("y") | missing("z")) stop("missing argument: X, Y or Z."))
  try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
  if(is.vector(x)){x <- as.matrix(x)}
  if(is.vector(y)){y <- as.matrix(y)}
  if(is.vector(z)){z <- as.matrix(z)}
  if((dim(x)[1] != dim(y)[1]) | (dim(x)[1] != dim(z)[1]) | (dim(y)[1] != dim(z)[1])) {stop("unmatched dimension: sample sizes of X, Y and Z are supposed to be equal.")}
  
  if(dim(z)[2] == 1){
    p <- c(dim(x)[2], dim(y)[2])
    n <- dim(x)[1]
    d <- 2
  }
  else{
    p <- c(dim(x)[2], dim(y)[2], dim(z)[2])
    n <- dim(x)[1]
    d <- 3
  }
  
  Delta <- tensor_mean(x,y,z)
  Delta <- as.tensor(Delta)
  
  try(if(missing("sparse_mode")) sparse_mode<-rep(TRUE, d))
  
  p.prod <- prod(p)
  p.minus.k <- p.prod / p
  s.prod <- prod(s)
  s.minus.k <- s.prod / s
  if(is.atomic(r) && length(r)==1){
    r <- rep(r, d)
  }
  r.minus.k <- prod(r) / r
  
  try(if(d != length(r)) stop("invalid input: r and the order of Y is incompatible."))
  
  Gamma_list <- list(); I_0 <- list()
  for(k in 1:d){
    a <- ifelse(k == 3, a2, a1)
    if(isTRUE(sparse_mode[k])){
      Delta_k <- k_unfold(Delta, k)@data
      Delta_k_row_max <- apply(abs(Delta_k), 1, max)
      thrd.k <- sqrt(a*log(p.prod)/n)
      I_0 <- c(I_0, list(Delta_k_row_max > thrd.k))
    }
    else{
      I_0 <- c(I_0, list(rep(TRUE,p[k])))
    }
  }
  
  Delta_sparse <- Delta
  for (k in 1:d){
    Delta_sparse <- ttm(Delta_sparse, diag(I_0[[k]]*1), k)
  }
  for (k in 1:d){
    if(isTRUE(sparse_mode[k])){
      datai <- k_unfold(Delta_sparse, k)@data
      selec <- (1:p[k])[I_0[[k]]]
      if(length(selec)<r[k]){
        newdata <- matrix(0, nrow=nrow(datai), ncol = ncol(datai))
        newdata[selec,] <- datai[selec,]
        Gamma <- (svd(newdata)$u[,1:r[k]])
      }
      else{
        Gamma <- matrix(0, nrow = p[k], ncol = r[k])
        Gamma[I_0[[k]],] <- svd(datai[I_0[[k]],,drop=FALSE])$u[,1:r[k]]
      }
      Gamma_list <- c(Gamma_list, list(t(Gamma)))
    }
    else{
      Gamma_list <- c(Gamma_list, list(t(svd(k_unfold(Delta_sparse, k)@data)$u[,1:r[k]])))
    }
  }
  
  t <- 1
  approx <- -1
  sv_k <- 0
  errors <- c()
  while(t<=tmax){
    error <- c()
    for(k in 1:d){
      a <- ifelse(k == 3, a2, a1)
      
      if(!isTRUE(sparse_mode[k])){
        A <- ttl(Delta, Gamma_list[-k], (1:d)[-k])
        A_k <- k_unfold(A, k)@data
        A_k.svd <- svd(A_k)
        Gamma_list[[k]] <- t(A_k.svd$u[,1:r[k]])
        sv_k <- A_k.svd$d[1:r[k]]
      }
      else{
        A <- ttl(Delta, Gamma_list[-k], (1:d)[-k])
        A_k <- k_unfold(A, k)@data
        eta_k <- (a*s.minus.k[k]*log(p.prod))/n
        I_k <- apply(A_k, 1, function(x){sum(x^2)}) > eta_k
        selec <- (1:p[k])[I_k]
        
        # If s_k < r_k, do not update
        if(length(selec) < r[k]){
          error <- c(error, 0)
          next
        }
        
        B_k <- A_k[I_k,,drop=FALSE]
        ###
        Gamma <- matrix(0, nrow(A_k), r[k])
        B_k.svd <- svd(B_k)
        Gamma[I_k,] <- B_k.svd$u[,1:r[k]]
        error <- c(error, subspace(t(Gamma_list[[k]]), Gamma))  # Save the algorithmic error for each basis matrix.
        Gamma_list[[k]] <- t(Gamma)
        sv_k <- B_k.svd$d[1:r[k]]
        ###
      }
    }
    errors <- rbind(errors, error)
    # Until the iteration reaches the t_max.
    t <- t+1
    
  }
  
  for(k in 1:d){
    Gamma_list[[k]] <- t(Gamma_list[[k]])
  }
  
  list(Gamma_hat = Gamma_list, errors = errors)
}



######################################################
# Save statistical error over iterations. The function is based on function STATSVD, please refer to the comments in STATSVD for more details.
# ###########
# Input:
# Note: The input arguments are similar to those in the function STATSVD.
# true_Gamma: the true basis matrices in the simulation model, used in computing the statistical error in each iteration.
######################################################
STATSVD.stat.error <- function(x, y, z, r, s, true_Gamma, tmax = 10, sparse_mode, vartol = 1e-6, a1 = 1, a2 = 1){
  
  try(if(missing("x") | missing("y") | missing("z")) stop("missing argument: X, Y or Z."))
  try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
  if(is.vector(x)){x <- as.matrix(x)}
  if(is.vector(y)){y <- as.matrix(y)}
  if(is.vector(z)){z <- as.matrix(z)}
  if((dim(x)[1] != dim(y)[1]) | (dim(x)[1] != dim(z)[1]) | (dim(y)[1] != dim(z)[1])) {stop("unmatched dimension: sample sizes of X, Y and Z are supposed to be equal.")}
  
  if(dim(z)[2] == 1){
    p <- c(dim(x)[2], dim(y)[2])
    n <- dim(x)[1]
    d <- 2
  }
  else{
    p <- c(dim(x)[2], dim(y)[2], dim(z)[2])
    n <- dim(x)[1]
    d <- 3
  }
  
  Delta <- tensor_mean(x,y,z)
  Delta <- as.tensor(Delta)
  
  try(if(missing("sparse_mode")) sparse_mode<-rep(TRUE, d))
  
  p.prod <- prod(p)
  p.minus.k <- p.prod / p
  s.prod <- prod(s)
  s.minus.k <- s.prod / s
  if(is.atomic(r) && length(r)==1){
    r <- rep(r, d)
  }
  r.minus.k <- prod(r) / r
  
  try(if(d != length(r)) stop("invalid input: r and the order of Y is incompatible."))
  
  errors <- c()
  Gamma_list <- list(); I_0 <- list()
  for(k in 1:d){
    a <- ifelse(k == 3, a2, a1)
    if(isTRUE(sparse_mode[k])){
      Delta_k <- k_unfold(Delta, k)@data
      Delta_k_row_max <- apply(abs(Delta_k), 1, max)
      thrd.k <- sqrt(a*log(p.prod)/n)
      I_0 <- c(I_0, list(Delta_k_row_max > thrd.k))
    }
    else{
      I_0 <- c(I_0, list(rep(TRUE,p[k])))
    }
  }
  
  Delta_sparse <- Delta
  for (k in 1:d){
    Delta_sparse <- ttm(Delta_sparse, diag(I_0[[k]]*1), k)
  }
  for (k in 1:d){
    if(isTRUE(sparse_mode[k])){
      datai <- k_unfold(Delta_sparse, k)@data
      selec <- (1:p[k])[I_0[[k]]]
      if(length(selec)<r[k]){
        newdata <- matrix(0, nrow=nrow(datai), ncol = ncol(datai))
        newdata[selec,] <- datai[selec,]
        Gamma <- (svd(newdata)$u[,1:r[k]])
      }
      else{
        Gamma <- matrix(0, nrow = p[k], ncol = r[k])
        Gamma[I_0[[k]],] <- svd(datai[I_0[[k]],,drop=FALSE])$u[,1:r[k]]
      }
      Gamma_list <- c(Gamma_list, list(t(Gamma)))
    }
    else{
      Gamma_list <- c(Gamma_list, list(t(svd(k_unfold(Delta_sparse, k)@data)$u[,1:r[k]])))
    }
  }
  # Save the statistical error at the initialization step for each basis matrix.
  error <- sapply(1:d, function(k){
    subspace(t(Gamma_list[[k]]), true_Gamma[[k]])
  })
  errors <- rbind(errors, error)
  ###
  
  t <- 1
  approx <- -1
  sv_k <- 0
  while(t<=tmax){
    error <- c()
    for(k in 1:d){
      a <- ifelse(k == 3, a2, a1)
      
      if(!isTRUE(sparse_mode[k])){
        A <- ttl(Delta, Gamma_list[-k], (1:d)[-k])
        A_k <- k_unfold(A, k)@data
        A_k.svd <- svd(A_k)
        Gamma_list[[k]] <- t(A_k.svd$u[,1:r[k]])
        sv_k <- A_k.svd$d[1:r[k]]
      }
      else{
        A <- ttl(Delta, Gamma_list[-k], (1:d)[-k])
        A_k <- k_unfold(A, k)@data
        eta_k <- (a*s.minus.k[k]*log(p.prod))/n
        I_k <- apply(A_k, 1, function(x){sum(x^2)}) > eta_k
        selec <- (1:p[k])[I_k]
        
        if(length(selec) < r[k]){
          next
        }
        
        B_k <- A_k[I_k,,drop=FALSE]
        Gamma <- matrix(0, nrow(A_k), r[k])
        B_k.svd <- svd(B_k)
        Gamma[I_k,] <- B_k.svd$u[,1:r[k]]
        Gamma_list[[k]] <- t(Gamma)
        sv_k <- B_k.svd$d[1:r[k]]
        ###
      }
    }
    ## Save the statistical error at each iteration for each basis matrix.
    error <- sapply(1:d, function(k){
      subspace(t(Gamma_list[[k]]), true_Gamma[[k]])
    })
    errors <- rbind(errors, error)
    ###
    
    # Until the iteration approches the t_max.
    t <- t+1
  }
  
  for(k in 1:d){
    Gamma_list[[k]] <- t(Gamma_list[[k]])
  }
  
  list(Gamma_hat = Gamma_list, errors = errors)
}
