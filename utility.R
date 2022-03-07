## Constructing AR covariance matrix
AR <- function(rho, p){
  mat <- diag(rep(1,p))
  for (i in seq_len(p)){
    for (j in seq_len(p)){
      mat[i,j] <- rho^(abs(i-j))
    }
  }
  mat
}

## Subspace distance defined in Section 5.1
subspace <- function(A,B){
  Pa <- qr.Q(qr(A))
  Pa <- Pa %*% t(Pa)
  Pb <- qr.Q(qr(B))
  Pb <- Pb %*% t(Pb)
  if(is.vector(A) || is.vector(B)){
    d <- 1
  }
  else{d <- dim(A)[2]}
  sqrt(sum((Pa-Pb)^2))/sqrt(2*d)
}

## Tensor mean of the three-way outer-product
tensor_mean <- function(x, y, z) {
  if(is.vector(x)){x <- as.matrix(x)}
  if(is.vector(y)){y <- as.matrix(y)}
  if(is.vector(z)){z <- as.matrix(z)}
  n <- dim(x)[1]
  result <- array(0, dim = c(dim(x)[2], dim(y)[2], dim(z)[2], n))
  for (i in seq_len(n)){
    result[,,,i] <- outer(outer(x[i,], y[i,]), z[i,])
  }
  drop(apply(result, c(1,2,3), mean))
}

## tanh function
tanh_func <- function(rho, xi, x){
  rho*(2/(1+exp(-2*xi*x)) - 1)
}

## Variable selection function for ULA
ula.variable <- function(ula, p, s){
  var.list <- lapply(1:length(p), function(i){
    tmp <- k_unfold(ula, i) # mode-i matricization of tensor ula.
    row_norm <- apply(tmp@data, 1, function(x){sum(x^2)}) # The l2-norm of each row
    order(row_norm, decreasing = TRUE)[1:s[i]] # The index for the rows with the largest s l2-norm.
  })
  var.list
}

## Generate index in different folders for cross-validation. (only used in `AD.R`)
folds.gen <- function(nfolds, n){
  seq <- sample(rep(1:nfolds, length=n), size = n, replace = FALSE)
  lab.list <- list()
  for (i in seq_len(nfolds)){
    lab.list[[i]] <- list(which(seq != i), which(seq == i))
  }
  lab.list
}
