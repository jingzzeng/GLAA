# ---------------------- Model 1 (model setting in Scenario 1) ---------------------- #
Model1 <- function(pp = 100){
  p <- c(pp,pp)
  s <- c(5,5)
  r <- c(2,2)
  sparse_mode <- c(TRUE, TRUE)
  n <- 500
  # The proper tuning parameter candidates for each dimension pp.
  if(pp==100){a.list <- seq(0.1, 2, length.out = 10)}
  else if(pp==200){a.list <- seq(0.1, 2, length.out = 10)}
  else if(pp==300){a.list <- seq(0.3, 2, length.out = 10)}
  else if(pp==400){a.list <- seq(0.35, 0.45, length.out = 10)}
  else if(pp==500){a.list <- seq(0.35, 0.45, length.out = 10)}
  
  sigma.x <- rbind(cbind(AR(0.3, s[1]), matrix(0, s[1], p[1] - s[1])), 
                   cbind(matrix(0, p[1]-s[1], s[1]), diag(1, p[1]-s[1], p[1]-s[1]))) # Sigma_X
  sigma.y <- rbind(cbind(AR(0.3, s[2]), matrix(0, s[2], p[2] - s[2])), 
                   cbind(matrix(0, p[2]-s[2], s[2]), diag(1, p[2]-s[2], p[2]-s[2]))) # Sigma_Y
  sigma.x.half <- sqrtm(sigma.x)$B # Simga_X^{1/2}
  sigma.y.half <- sqrtm(sigma.y)$B # Simga_Y^{1/2}
  
  alpha <- matrix(0, p[1], 2)
  alpha[1:s[1], 1] <- c(1,1,1,1,1)
  alpha[1:s[1], 2] <- c(0,0,0,-1,1)
  alpha[,1] <- alpha[,1]/sqrt(sum((alpha[,1])^2))
  alpha[,2] <- alpha[,2]/sqrt(sum((alpha[,2])^2))
  beta <- matrix(0, p[2], 2)
  beta[1:s[2], 1] <- c(1,1,1,1,1)
  beta[1:s[2], 2] <- c(0,0,0,-1,1)
  beta[,1] <- beta[,1]/sqrt(sum((beta[,1])^2))
  beta[,2] <- beta[,2]/sqrt(sum((beta[,2])^2))
  Gamma.x <- sigma.x.half %*% alpha # Gamma_1 for X
  Gamma.y <- sigma.y.half %*% beta # Gamma_2 for Y
  
  Gamma <- list(Gamma.x, Gamma.y)
  
  data.gen <- function(n){
    data <- lapply(seq_len(n), function(i){
      z <- rnorm(1)
      id.z <- ifelse(z > 0, 1, -1)
      rho.z <- id.z * diag(c(0.95,0.85), 2,2) # f matrix
      mean.xy.z <- Gamma.x%*%rho.z%*%t(Gamma.y) # Cov(X, Y)
      cov.xy.z <- rbind(cbind(sigma.x, mean.xy.z), cbind(t(mean.xy.z), sigma.y)) # Joint covariance matrix of (X, Y).
      xy <- mvrnorm(1, rep(0, p[1]+p[2]), cov.xy.z) # Generate (X, Y).
      x <- xy[1:p[1]]
      y <- xy[(p[1]+1):(p[1]+p[2])]
      list(x = x, y = y, z = z)
    })
    x <- do.call(rbind, lapply(data, '[[', 1))
    y <- do.call(rbind, lapply(data, '[[', 2))
    z <- do.call(rbind, lapply(data, '[[', 3))
    list(x = x, y = y, z = z)
  }
  
  list(data.gen = data.gen, Gamma = Gamma, a.list = a.list,  p = p, r = r, s = s, n = n, sparse_mode = sparse_mode)
}
# ------------------------------------------------ #


# ---------------------- Model 2 (model setting in Scenario 2) ---------------------- #
Model2 <- function(pp = 20){
  p <- c(100,100,pp)
  s <- c(5,5,5)
  r <- c(2,2,1)
  sparse_mode <- c(TRUE, TRUE, TRUE)
  n <- 500
  a1.list <- seq(0.05, 0.2, length.out = 5)
  a2.list <- seq(0.01,0.25, length.out = 10)

  sigma.x <- rbind(cbind(AR(0.3, s[1]), matrix(0, s[1], p[1] - s[1])),
                   cbind(matrix(0, p[1]-s[1], s[1]), diag(1, p[1]-s[1], p[1]-s[1])))
  sigma.y <- rbind(cbind(AR(0.3, s[2]), matrix(0, s[2], p[2] - s[2])),
                   cbind(matrix(0, p[2]-s[2], s[2]), diag(1, p[2]-s[2], p[2]-s[2])))
  sigma.z <- diag(1, p[3], p[3])
  sigma.x.half <- sqrtm(sigma.x)$B
  sigma.y.half <- sqrtm(sigma.y)$B

  alpha <- matrix(0, p[1], 2)
  alpha[1:s[1], 1] <- c(1,1,1,1,1)
  alpha[1:s[1], 2] <- c(0,0,0,-1,1)
  alpha[,1] <- alpha[,1]/sqrt(sum((alpha[,1])^2))
  alpha[,2] <- alpha[,2]/sqrt(sum((alpha[,2])^2))
  beta <- matrix(0, p[2], 2)
  beta[1:s[2], 1] <- c(1,1,1,1,1)
  beta[1:s[2], 2] <- c(0,0,0,-1,1)
  beta[,1] <- beta[,1]/sqrt(sum((beta[,1])^2))
  beta[,2] <- beta[,2]/sqrt(sum((beta[,2])^2))
  Gamma.x <- sigma.x.half %*% alpha
  Gamma.y <- sigma.y.half %*% beta

  Gamma.z <- matrix(0, p[3], 1)
  Gamma.z[1:s[3], 1] <- c(1,1,1,1,1)

  Gamma <- list(Gamma.x, Gamma.y, Gamma.z)

  data.gen <- function(n){
    data <- lapply(seq_len(n), function(i){
      z <- mvrnorm(1, rep(0,p[3]), sigma.z)
      id.z <- ifelse(c(t(Gamma.z) %*% z) > 0, 1, -1)
      rho.z <- id.z * diag(c(0.95,0.85), 2,2)
      mean.xy.z <- Gamma.x%*%rho.z%*%t(Gamma.y)
      cov.xy.z <- rbind(cbind(sigma.x, mean.xy.z), cbind(t(mean.xy.z), sigma.y))
      xy <- mvrnorm(1, rep(0, p[1]+p[2]), cov.xy.z)
      x <- xy[1:p[1]]
      y <- xy[(p[1]+1):(p[1]+p[2])]
      list(x=x, y=y, z=z)
    })
    x <- do.call(rbind, lapply(data, '[[', 1))
    y <- do.call(rbind, lapply(data, '[[', 2))
    z <- do.call(rbind, lapply(data, '[[', 3))
    list(x = x, y = y, z = z)
  }

  list(data.gen = data.gen, Gamma = Gamma, a1.list = a1.list, a2.list = a2.list, p = p, r = r, s = s, n = n, sparse_mode = sparse_mode)
}
# ------------------------------------------------ #


# # ---------------------- Model 3 (model setting in Scenario 3) ---------------------- #
Model3 <- function(n = 20){
  p <- c(100, 25)
  s <- c(5,5)
  r <- c(2,2)
  sparse_mode <- c(TRUE, TRUE)
  n <- n
  if(n ==60){a.list <- seq(0.1, 2, length.out = 50)
  }else{a.list <- seq(0.4, 0.6, length.out = 20)}
  
  sigma.x <- rbind(cbind(AR(0.3, s[1]), matrix(0, s[1], p[1] - s[1])), 
                   cbind(matrix(0, p[1]-s[1], s[1]), diag(1, p[1]-s[1], p[1]-s[1])))
  sigma.y <- rbind(cbind(AR(0.3, s[2]), matrix(0, s[2], p[2] - s[2])), 
                   cbind(matrix(0, p[2]-s[2], s[2]), diag(1, p[2]-s[2], p[2]-s[2])))
  sigma.x.half <- sqrtm(sigma.x)$B
  sigma.y.half <- sqrtm(sigma.y)$B
  
  alpha <- matrix(0, p[1], 2)
  alpha[1:s[1], 1] <- c(1,1,1,1,1)
  alpha[1:s[1], 2] <- c(0,0,0,-1,1)
  alpha[,1] <- alpha[,1]/sqrt(sum((alpha[,1])^2))
  alpha[,2] <- alpha[,2]/sqrt(sum((alpha[,2])^2))
  beta <- matrix(0, p[2], 2)
  beta[1:s[2], 1] <- c(1,1,1,1,1)
  beta[1:s[2], 2] <- c(0,0,0,-1,1)
  beta[,1] <- beta[,1]/sqrt(sum((beta[,1])^2))
  beta[,2] <- beta[,2]/sqrt(sum((beta[,2])^2))
  Gamma.x <- sigma.x.half %*% alpha
  Gamma.y <- sigma.y.half %*% beta
  
  Gamma <- list(Gamma.x, Gamma.y)
  
  data.gen <- function(n){
    data <- lapply(seq_len(n), function(i){
      z <- rnorm(1)
      id.z <- ifelse(z > 0, 1, -1)
      rho.z <- id.z * diag(c(0.95,0.85), 2,2)
      mean.xy.z <- Gamma.x%*%rho.z%*%t(Gamma.y)
      cov.xy.z <- rbind(cbind(sigma.x, mean.xy.z), cbind(t(mean.xy.z), sigma.y))
      xy <- mvrnorm(1, rep(0, p[1]+p[2]), cov.xy.z)
      x <- xy[1:p[1]]
      y <- xy[(p[1]+1):(p[1]+p[2])]
      list(x=x, y=y, z=z)
    })
    x <- do.call(rbind, lapply(data, '[[', 1))
    y <- do.call(rbind, lapply(data, '[[', 2))
    z <- do.call(rbind, lapply(data, '[[', 3))
    list(x = x, y = y, z = z)
  }
  list(data.gen = data.gen, Gamma = Gamma, a.list = a.list,  p = p, r = r, s = s, n = n, sparse_mode = sparse_mode)
}


# ---------------------- Model 4 (model setting in Scenario 1, with a tanh-type function f and xi = 1) ---------------------- #
Model4 <- function(pp = 100){
  p <- c(pp,pp)
  s <- c(5,5)
  r <- c(2,2)
  sparse_mode <- c(TRUE, TRUE)
  n <- 500
  
  if(pp==100){a.list <- seq(0.2, 2, length.out = 10)}
  else if(pp==200){a.list <- seq(0.2, 2, length.out = 10)}
  else if(pp==300){a.list <- seq(0.3, 0.5, length.out = 10)}
  else if(pp==400){a.list <- seq(0.3, 0.5, length.out = 10)}
  else if(pp==500){a.list <- seq(0.3, 0.5, length.out = 10)}
  else{a.list <- seq(0.2, 2, length.out = 10)}
  
  sigma.x <- rbind(cbind(AR(0.3, s[1]), matrix(0, s[1], p[1] - s[1])), 
                   cbind(matrix(0, p[1]-s[1], s[1]), diag(1, p[1]-s[1], p[1]-s[1])))
  sigma.y <- rbind(cbind(AR(0.3, s[2]), matrix(0, s[2], p[2] - s[2])), 
                   cbind(matrix(0, p[2]-s[2], s[2]), diag(1, p[2]-s[2], p[2]-s[2])))
  sigma.x.half <- sqrtm(sigma.x)$B
  sigma.y.half <- sqrtm(sigma.y)$B
  
  alpha <- matrix(0, p[1], 2)
  alpha[1:s[1], 1] <- c(1,1,1,1,1)
  alpha[1:s[1], 2] <- c(0,0,0,-1,1)
  alpha[,1] <- alpha[,1]/sqrt(sum((alpha[,1])^2))
  alpha[,2] <- alpha[,2]/sqrt(sum((alpha[,2])^2))
  beta <- matrix(0, p[2], 2)
  beta[1:s[2], 1] <- c(1,1,1,1,1)
  beta[1:s[2], 2] <- c(0,0,0,-1,1)
  beta[,1] <- beta[,1]/sqrt(sum((beta[,1])^2))
  beta[,2] <- beta[,2]/sqrt(sum((beta[,2])^2))
  Gamma.x <- sigma.x.half %*% alpha
  Gamma.y <- sigma.y.half %*% beta
  
  Gamma <- list(Gamma.x, Gamma.y)
  
  data.gen <- function(n){
    data <- lapply(seq_len(n), function(i){
      z <- rnorm(1)
      rho.z <- diag(c(tanh_func(0.95, 1, z), tanh_func(0.85, 1, z)), r[1], r[2]) # tanh
      mean.xy.z <- Gamma.x%*%rho.z%*%t(Gamma.y)
      cov.xy.z <- rbind(cbind(sigma.x, mean.xy.z), cbind(t(mean.xy.z), sigma.y))
      xy <- mvrnorm(1, rep(0, p[1]+p[2]), cov.xy.z)
      x <- xy[1:p[1]]
      y <- xy[(p[1]+1):(p[1]+p[2])]
      list(x=x, y=y, z=z)
    })
    x <- do.call(rbind, lapply(data, '[[', 1))
    y <- do.call(rbind, lapply(data, '[[', 2))
    z <- do.call(rbind, lapply(data, '[[', 3))
    list(x = x, y = y, z = z)
  }
  
  list(data.gen = data.gen, Gamma = Gamma, a.list = a.list,  p = p, r = r, s = s, n = n, sparse_mode = sparse_mode)
}
# ------------------------------------------------ #


# ---------------------- Model 5 (model setting in Scenario 1, with a tanh-type function f and xi = 5) ---------------------- #
Model5 <- function(pp = 100){
  p <- c(pp,pp)
  s <- c(5,5)
  r <- c(2,2)
  sparse_mode <- c(TRUE, TRUE)
  n <- 500
  
  if(pp==100){a.list <- seq(0.2, 2, length.out = 10)}
  else if(pp==200){a.list <- seq(0.2, 2, length.out = 10)}
  else if(pp==300){a.list <- seq(0.35, 0.5, length.out = 10)}
  else if(pp==400){a.list <- seq(0.35, 0.5, length.out = 10)}
  else if(pp==500){a.list <- seq(0.35, 0.5, length.out = 10)}
  else{a.list <- seq(0.2, 2, length.out = 10)}
  
  
  sigma.x <- rbind(cbind(AR(0.3, s[1]), matrix(0, s[1], p[1] - s[1])), 
                   cbind(matrix(0, p[1]-s[1], s[1]), diag(1, p[1]-s[1], p[1]-s[1])))
  sigma.y <- rbind(cbind(AR(0.3, s[2]), matrix(0, s[2], p[2] - s[2])), 
                   cbind(matrix(0, p[2]-s[2], s[2]), diag(1, p[2]-s[2], p[2]-s[2])))
  sigma.x.half <- sqrtm(sigma.x)$B
  sigma.y.half <- sqrtm(sigma.y)$B
  
  alpha <- matrix(0, p[1], 2)
  alpha[1:s[1], 1] <- c(1,1,1,1,1)
  alpha[1:s[1], 2] <- c(0,0,0,-1,1)
  alpha[,1] <- alpha[,1]/sqrt(sum((alpha[,1])^2))
  alpha[,2] <- alpha[,2]/sqrt(sum((alpha[,2])^2))
  beta <- matrix(0, p[2], 2)
  beta[1:s[2], 1] <- c(1,1,1,1,1)
  beta[1:s[2], 2] <- c(0,0,0,-1,1)
  beta[,1] <- beta[,1]/sqrt(sum((beta[,1])^2))
  beta[,2] <- beta[,2]/sqrt(sum((beta[,2])^2))
  Gamma.x <- sigma.x.half %*% alpha
  Gamma.y <- sigma.y.half %*% beta
  
  Gamma <- list(Gamma.x, Gamma.y)
  
  data.gen <- function(n){
    data <- lapply(seq_len(n), function(i){
      z <- rnorm(1)
      rho.z <- diag(c(tanh_func(0.95, 5, z), tanh_func(0.85, 5, z)), r[1], r[2]) # tanh
      mean.xy.z <- Gamma.x%*%rho.z%*%t(Gamma.y)
      cov.xy.z <- rbind(cbind(sigma.x, mean.xy.z), cbind(t(mean.xy.z), sigma.y))
      xy <- mvrnorm(1, rep(0, p[1]+p[2]), cov.xy.z)
      x <- xy[1:p[1]]
      y <- xy[(p[1]+1):(p[1]+p[2])]
      list(x=x, y=y, z=z)
    })
    x <- do.call(rbind, lapply(data, '[[', 1))
    y <- do.call(rbind, lapply(data, '[[', 2))
    z <- do.call(rbind, lapply(data, '[[', 3))
    list(x = x, y = y, z = z)
  }
  
  list(data.gen = data.gen, Gamma = Gamma, a.list = a.list,  p = p, r = r, s = s, n = n, sparse_mode = sparse_mode)
}
# ------------------------------------------------ #
