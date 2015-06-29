# Libraries
require(MASS) # for normal bivariate

# ===================================================================
# Functions to build the matrix M
# ===================================================================
# Sample pairs for random communities and "Eyeball" (food web) matrices
SampleNormal <- function(NumPairs, mu, sigma, rho, Eyeball = FALSE){
  # Take pairs from a bivariate normal distribution
  if (Eyeball == FALSE) {
    # the marginals have the same mean
    mus <- c(mu, mu)
  } else {
    # the marginals have means with different signs, 
    # with mu_upper <- 3 * |mu_lower|
    mus <- c(-mu, 3 * mu)     
  }
  covariance.matrix <- matrix(sigma^2 * c(1, rho, rho, 1), 2, 2)
  Pairs <- mvrnorm(NumPairs, mus, covariance.matrix)
  return(Pairs)
}

Sample4C <- function(NumPairs, mu, sigma, rho, Eyeball = FALSE){
  # Take pairs from the four-corner distribution of Allesina et al Nature Communication 2015
  gamma <- (rho + 1) / 2
  PosNeg <- sign(rnorm(NumPairs))
  SwitchSign <- sign(runif(NumPairs, -(1 - gamma), gamma))
  if (Eyeball == FALSE){
    # The marginals have the same mean
    Pairs <- cbind(mu + PosNeg * sigma, mu + PosNeg * SwitchSign * sigma)
  } else {
    # the marginals have means with different signs, 
    # with mu_upper <- 3 * |mu_lower|
    mus <- c(-mu, 3 * mu)
    Pairs <- cbind(mus[1] + PosNeg * sigma, mus[2] + PosNeg * SwitchSign * sigma)  
  }
  return(Pairs)
}

# Build matrices
BuildMatrices <- function(S, connectance, a, mu, sigma, rho, Q, Distr = "Normal", Eyeball = FALSE){
  parameters <- data.frame()
  matrices <- list(M = matrix(0, 1, 1), A = matrix(0, 1, 1), B = matrix(0, 1, 1))
  
  # Matrix to be built
  M <- matrix(0, S, S)
  # Vector storing membership
  membership <- rep(0, S)
  # Build data frame for parameters
  parameters <- data.frame(S = S, connectance = connectance, a = a, mu = mu, sigma = sigma, 
                           rho = rho, Q = Q, Distr = Distr, Eyeball = Eyeball, Cw = NA, Cb = NA)
  # First, check whether this level of modularity is viable: Cw and Cb should be bounded by 0 and 1
  Cw <- round(connectance * (Q / (a^2 + (1-a)^2) + 1), 6)
  Cb <- round(connectance * (-Q / (2 * a * (1-a)) + 1), 6)
  if ((Cw > 1) | (Cb > 1) | (Cw < 0) | (Cb < 0)){
    print(paste("Unfeasible Cw", Cw, "Cb", Cb))
    return(list(parameters = parameters, membership = membership, matrices = matrices))
  } else {
    parameters$Cw <- Cw
    parameters$Cb <- Cb
  }
  # M = W * K 
  # where W is a matrix of weights, and K the adjacency matrix of an undirected graph
  ###################################################################################
  # Matrix W
  ###################################################################################
  W <- matrix(0, S, S)
  # sample the pairs
  if (Distr == "4C"){
    Pairs <- Sample4C(S * (S-1) / 2, mu, sigma, rho, Eyeball)
  } else {
    Pairs <- SampleNormal(S * (S-1) / 2, mu, sigma, rho, Eyeball)
  }
  # arrange the pairs
  W[upper.tri(W)] <- Pairs[,1]
  W <- t(W)
  W[upper.tri(W)] <- Pairs[,2]
  W <- t(W)
  ###################################################################################
  # Matrix K
  ###################################################################################
  K <- matrix(0, S, S)
  # vector or membership
  membership <- (0:(S-1)) / S
  membership[membership >= a] <- 2
  membership[membership < 1] <- 1
  size1 <- sum(membership == 1)
  size2 <- sum(membership == 2)
  # If this is an Eyeball matrix, scramble the pairs
  if (Eyeball == TRUE) {
    membership <- sample(membership)
  }
  K <- matrix(0, S, S)
  K11 <- K[membership == 1, membership == 1]
  K11[upper.tri(K11)] <- sample(c(
    rep(1, round(Cw * size1 * (size1 - 1) / 2)),
    rep(0, size1 * (size1 - 1) / 2))[1:(size1 * (size1 - 1) / 2)])
  K22 <- K[membership == 2, membership == 2]
  K22[upper.tri(K22)] <- sample(c(
    rep(1, round(Cw * size2 * (size2 - 1) / 2)),
    rep(0, size2 * (size2 - 1) / 2))[1:(size2 * (size2 - 1) / 2)])
  K12 <- K[membership == 1, membership == 2]
  K12 <- sample(c(rep(1, round(Cb * size1 * size2 )),
                  rep(0, size1 * size2 ))[1:(size1 * size2)])
  K[membership == 1, membership == 1] <- K11
  K[membership == 2, membership == 2] <- K22
  K[membership == 1, membership == 2] <- K12
  K <- K + t(K)
  
  ###################################################################################
  # Matrices M, A, and B
  ###################################################################################
  # Matrix M
  M <- W * K
  
  # Matrix A
  A <- M
  A[membership == 1, membership == 1] <- parameters$mu * parameters$Cw
  A[membership == 1, membership == 2] <- parameters$mu * parameters$Cb
  A[membership == 2, membership == 1] <- parameters$mu * parameters$Cb
  A[membership == 2, membership == 2] <- parameters$mu * parameters$Cw
  
  # Matrix B
  B <- M - A
  
  # populate the list of matrices
  matrices$M <- M
  matrices$A <- A
  matrices$B <- B
  
  return(list(parameters = parameters, membership = membership, matrices = matrices))
}