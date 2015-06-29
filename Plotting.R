require(ggplot2)

GetEllipse <- function(centerx, radiusx, radiusy){
  thetas <- seq(pi / 2.0, 0.0, length = 1000)
  xbase <- radiusx * cos(thetas)
  ybase <- radiusy * sin(thetas)
  x <- c(xbase, rev(xbase), -xbase, rev(-xbase))
  y <- c(ybase, rev(-ybase), -ybase, rev(ybase))
  return(data.frame(Real = x + centerx, Imaginary = y))
}

