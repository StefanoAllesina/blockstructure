source("Plotting.R")

Predictions <- function(matobj, plot = TRUE){
  ReL1.observed <- 0
  ReL1.predicted <- NA
  Re.bulk <- NA
  Re.outlier <- NA
  
  # first, calculate eigenvalues of M
  M <- matobj$matrices$M
  eM <- eigen(M, only.values = TRUE)$values
  # Store the max real part
  ReL1.observed <- max(Re(eM))
  
  if (plot == TRUE){
    deig <- data.frame(Real = Re(eM), Imaginary = Im(eM))
    dpoints <- data.frame(Real = numeric(0), Imaginary = numeric(0))
    dpoly <- data.frame(Real = numeric(0), Imaginary = numeric(0))
  }
  ###############################################################
  # Check whether it is possible to predict the eigenvalues
  ###############################################################
  # Extract the parameters
  pr <- matobj$parameters
  # get the stats needed for the calculations
  
  mu <- pr$mu
  sigma <- pr$sigma
  rho <- pr$rho
  
  S <- pr$S
  connectance <- pr$connectance
  a <- pr$a
  Q <- pr$Q
  
  Cw <- pr$Cw
  Cb <- pr$Cb
  
  E.w <- Cw * mu
  V.w <- Cw * (sigma^2 + (1 - Cw) * mu^2)
  rho.w <- (sigma^2 * rho + (1 - Cw) * mu^2) / (sigma^2 + (1 - Cw) * mu^2)
  
  # Effective parameters
  mu.b <- Cb * mu
  sigma.b <- sqrt(Cb * (sigma^2 + (1 - Cb) * mu^2))
  rho.b <- (sigma^2 * rho + (1 - Cb) * mu^2) / (sigma^2 + (1 - Cb) * mu^2)
  
  mu.w <- Cw * mu
  sigma.w <- sqrt(Cw * (sigma^2 + (1 - Cw) * mu^2))
  rho.w <- (sigma^2 * rho + (1 - Cw) * mu^2) / (sigma^2 + (1 - Cw) * mu^2)
  
  ########################################################################
  # Eigenvalues of A (outliers)
  ########################################################################
  A1 <- (S / 2) * (mu.w + sqrt((1 - 4 * a * (1 - a)) * mu.w^2 + 4 * a * (1 - a) * mu.b^2))
  A2 <- (S / 2) * (mu.w - sqrt((1 - 4 * a * (1 - a)) * mu.w^2 + 4 * a * (1 - a) * mu.b^2))
  
  ########################################################################
  # Eigenvalues of B (bulk), and correction for outliers
  ########################################################################
  if (a == 0.5){
    # Equal sized-communities
    rho.tilde <- (rho * sigma^2 + (1 - connectance - 4 * connectance * Q) * mu^2) / (sigma^2 + (1 - connectance - 4 * connectance * Q) * mu^2)
    sigma.tilde <- sqrt(connectance * (sigma ^2 + (1 - connectance - 4 * connectance * Q) * mu^2))
    Re.bulk <- sqrt(S) * sigma.tilde * (1 + rho.tilde) - mu.w
    # correction for outliers
    outlier1 <- A1
    outlier2 <- A2
    if (abs(A1) > sqrt(S) * sigma.tilde) outlier1 <- A1 + S * sigma.tilde^2 * rho.tilde / A1
    if (abs(A2) > sqrt(S) * sigma.tilde) outlier2 <- A2 + S * sigma.tilde^2 * rho.tilde / A2
    # Predictions
    Re.outlier <- max(outlier1, outlier2)
    ReL1.predicted <- max(Re.outlier, Re.bulk)  
    # For plotting
    if (plot == TRUE){
      dpoints <- rbind(dpoints, data.frame(Real = c(outlier1, outlier2), Imaginary = c(0, 0)))
      dpoly <- GetEllipse(-mu.w, sqrt(S) * sigma.tilde * (1 + rho.tilde), sqrt(S) * sigma.tilde * (1 - rho.tilde))
    }
  }
  if ((a != 0.5) & round(Cw, 2) == 0){
    # Perfectly Bipartite case
    # using equation 24. First, check which case we should use
    if ((rho.b) < 0 & (a < 4 * rho.b^2 / (1 + 6 * rho.b^2 + rho.b^4))){
        Re.bulk <-  sqrt(S) * sigma.b * sqrt((-2 * a * rho.b * (1 - rho.b)^2 ) / (rho.b^2)) / (2 * sqrt(2))
    } else {
      Re.bulk <- sqrt(S) * sigma.b * sqrt(rho.b + sqrt(a * (1 - a)) * (1 + rho.b))  
    }
    # no correction for outliers
    outlier1 <- A1
    outlier2 <- A2
    Re.outlier <- max(outlier1, outlier2)
    ReL1.predicted <- max(Re.outlier, Re.bulk)  
    # For plotting
    if (plot == TRUE){
      # calculate the ellipse for xy
      ellXY <- GetEllipse(S * rho.b * sigma.b^2,
                          S * sigma.b^2 * sqrt(a * (1 - a)) * (1 + rho.b^2),
                          S * sigma.b^2 * sqrt(a * (1 - a)) * (1 - rho.b^2))
      ellXY.comp <- complex(real = ellXY$Real, imaginary = ellXY$Imaginary)
      # Now do a square-root transformation 
      ellB.support <- (sqrt(ellXY.comp))
      ellB <- data.frame(Real = Re(ellB.support), Imaginary = Im(ellB.support))
      ellB <- rbind(ellB, data.frame(Real = -Re(ellB.support), Imaginary = Im(ellB.support)))
      dpoly <- ellB
      dpoints <- rbind(dpoints, data.frame(Real = c(0, outlier1, outlier2), Imaginary = c(0, 0, 0)))
    }
  }
  if ((a != 0.5) & round(Cb, 2) == 0){
    # Perfectly Modular case
    # the smaller subsystem has size aS
    if (a > 0.5) {
      a <- 1-a
      tmp <- A1
      A1 <- A2
      A2 <- tmp
    }
    # calculate semi-axes for largest subsystem
    rx <- sqrt(S * (1 - a)) * sigma.w * (1 + rho.w)
    ry <- sqrt(S * (1 - a)) * sigma.w * (1 - rho.w)
    Re.bulk <- -mu.w + rx
    # correction for outliers
    outlier1 <- A1
    outlier2 <- A2
    if (abs(A1) > sqrt((1 - a) * S) * sigma.w) outlier1 <- A1 + S * (1 - a) * sigma.w^2 * rho.w / A1
    if (abs(A2) > sqrt(a * S) * sigma.w) outlier2 <- A2 + S * a * sigma.w^2 * rho.w / A2
    # Predictions
    Re.outlier <- max(outlier1, outlier2)
    ReL1.predicted <- max(Re.outlier, Re.bulk)  
    # For plotting
    if (plot == TRUE){
      dpoints <- rbind(dpoints, data.frame(Real = c(outlier1, outlier2), Imaginary = c(0, 0)))
      dpoly <- GetEllipse(-mu.w, rx, ry)
    }
  }
  if (plot == TRUE){
    pl <- ggplot(data = deig, aes(Real, Imaginary)) +
          geom_point(alpha = 0.75) + 
          geom_point(data = dpoints, aes(Real, Imaginary), alpha = 0.25, colour = "blue", size = 4) + 
          geom_polygon(data = dpoly, aes(Real, Imaginary), alpha = 0.25, colour = NA, fill = "blue") + 
          scale_x_continuous(expression(Re(lambda[i]))) + 
          scale_y_continuous(expression(Im(lambda[i]))) + 
          theme_bw()
    print(pl)
  }
  matobj$ReL1.observed <- ReL1.observed
  matobj$ReL1.predicted <- ReL1.predicted
  matobj$Re.bulk <- Re.bulk
  matobj$Re.outlier <- Re.outlier
  matobj$eigenvalues.M <- eM 
  return(matobj)
}