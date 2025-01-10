
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# correlation matrix ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# adjusted autoregressive correlation matrix
corar1 <- function(rho,p){
  times <- 1:p
  sigma <- 1
  H <- abs(outer(times, times, "-"))
  V <- sigma * rho^H
  p <- nrow(V)
  V[cbind(1:p, 1:p)] <- V[cbind(1:p, 1:p)] * sigma
  # first two variables not correlated
  V <- cbind(V[,1],V)
  V <- rbind(V[1,], V)
  V[1,2] <- V[2,1] <- 0
  
  return(V)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# dgp and missing from Little2004 - Robust Likelihood based analysis ---- 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dgp <- function(n=100,structure,p=0,rho=0, transform=TRUE,error = "gaussian"){
  Sigma <- corar1(rho,p+1)
  X <- mvtnorm::rmvnorm(n= n, sigma = Sigma)
  # transform
  if(transform){
    X[,1] <- 2*(pnorm(X[,1])-0.5)
    X[,2] <- 2*(pnorm(X[,2])-0.5)
  }
  X1 <- X[,1] 
  X2 <- X[,2]
  colnames(X) <- paste0("X", 1:ncol(X))
  if(structure == "constant"){
    y <- 10
  }else if(structure=="linear"){
    y <- 10*(1 + X1 + 3*X2)
  }else if(structure=="additive"){
    y <- 118 + (3*X1 - 3)^3 + (3*X2 - 3)^3 
  }else if(structure=="non-additive"){
    y <- 10*(1 + X1 + X2 + 4*X1*X2) 
  }else if(structure=="sine"){
    y <- 6 + 15*sin(pi*X1) 
  }
  y_mean <- mean(y)
  # add noise
  if(error == "binomial"){
    y <- as.factor(rbinom(n = n, size = 1, p = plogis(y)))
  }else{
    y <- y + rnorm(n = n, mean = 0, sd = 2)
  }

  return(list("data" = data.frame(y, X), "y_mean" = y_mean))
}

mis_little <- function(data,missing){
  if(missing == "constant"){
    miss_p <- plogis(0.5) 
  }else if(missing=="linear"){
    miss_p <- with(data, plogis(X1 + X2))
  }else if(missing=="additive"){
    miss_p <- with(data, plogis(X1^3 + X2^3))
  }else if(missing=="non-additive"){
    miss_p <- with(data, plogis(X1 + X2 + 3*X1*X2))
  }else if(missing=="sine"){
    miss_p <- plogis(1.5*sin(pi*data$X1)) 
  }
  mi_bin <- rbinom(n = nrow(data), size = 1, prob = miss_p)
  amp <- data
  amp$y[mi_bin==1] <- NA 
  return( list("data" = data, "amp" = amp))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# dgp and missing combined ---- 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dgp_and_miss <- function(n=100,structure,missing,p=0,rho=0, transform=TRUE,error = "gaussian"){
  data <- dgp(n,structure,p,rho, transform,error)
  y_mean <- data$y_mean
  data <- mis_little(data = data$data,missing)
  return(list(
    "data" = data$data
    ,"amp" = data$amp
    ,"y_mean" = y_mean
    ,"n_true" = nrow(data$data)
    ,"n_cc" = sum(!is.na(data$amp$y))
  ))
}

