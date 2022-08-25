## Function to create the density values for an RACD map

# Here, 'xlim' and 'ylim' give the range of x- and y-coordinates for the pixel centres in the map area. Also, 'results' refers to a set of MCMC samples generated using run.MCMC() and 'M' is the size of the superpopulation.
# We are assuming that each pixel has an area of 1, so that 'xg' and 'yg' below contain the pixel centres. Then, '(length(xg) - 1) * (length(yg) - 1)' below gives the number of pixels we are using in the map area.

no.movement.density.vector <- function(xlim, ylim, results, M) {

  ## Points at which local density will be estimated
  xg <- seq(xlim[1], xlim[2], by=1)
  yg <- seq(ylim[1], ylim[2], by=1)

  ## Extracting z-values
  # Names of variables that have been monitored
  names <- names(results[1,])
  # Extracting "z" values from MCMC results
  Z <- results[,grep("z", names)]

  ## Extracting activity centres
  # Extracting "s" values from MCMC results (i.e. extracting sampled activity centres)
  S <- results[,grep("s[^i]", names)]
  # x-coordinates of all activity centres
  Sx <- S[,1:M]
  # y-coordinates of all activity centres
  Sy <- S[,-(1:M)]

  ## For each MCMC iteration, storing the number of animals alive and with their activity centres in each cell
  # Number of pixel centres
  npix <-  (length(xg) - 1) * (length(yg) - 1)
  Dn.vals <-  matrix(0, nrow=nrow(results), ncol=npix)
  for (i in 1:nrow(results)) {
    if ((i %% 100) == 0) print(i) # Track progress
    Sxout <-  Sx[i,][Z[i,] == 1]
    Sxout <-  cut(Sxout, breaks=xg, include.lowest=TRUE)
    Syout <-  Sy[i,][Z[i,] == 1]
    Syout <-  cut(Syout, breaks=yg, include.lowest=TRUE)
    Dn.vals[i,] <-  as.vector(table(Sxout, Syout))
  }

  ## Posterior mean for number of activity centres in each pixel
  density.vector <- apply(Dn.vals, 2, mean)
}

## Function to generate pixel centres across map area

# 'xlim' and 'ylim' are described as above. In addition, x.pixels and y.pixels give the number of pixels being used in the x- and y-direction.
# This function won't necessarily produce pixel centres representing pixels with an area of 1 (depending on arguments provided).

library("spatstat")
centres <- function(xlim, ylim, x.pixels, y.pixels) {
  # Creating an object of class 'owin' representing our map area
  window.2 <- owin(xrange=xlim, yrange=ylim)
  # Generating our set of pixel centres
  points <- gridcentres(window.2, x.pixels, y.pixels)
  # Converting the result to a matrix
  centres <- as.matrix(cbind(points$x, points$y))
  # Printing the result
  centres
}

