## Function to create the density values for an RACD map (the resulting density values do not take the movement of animals into account, unlike RUD maps)

## Here, 'xlim' and 'ylim' give the range of x- and y-coordinates for the map area (see the comment preceding the 'centres()' function in 'RUDMaps_Functions.R' for further explanation on this). Also, 'results' refers to a set of MCMC samples generated using run.MCMC() and M is the size of the superpopulation.

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
  Dn.vals = matrix(0, nrow=10000, ncol=2500)
  for (i in 1:10000) {
    if ((i %% 100) == 0) print(i) # Track progress
    Sxout = Sx[i,][Z[i,] == 1]
    Sxout = cut(Sxout, breaks=xg, include.lowest=TRUE)
    Syout = Sy[i,][Z[i,] == 1]
    Syout = cut(Syout, breaks=yg, include.lowest=TRUE)
    Dn.vals[i,] = as.vector(table(Sxout, Syout))
  }

  ## Posterior mean for number of activity centres in each pixel
  density.vector <- apply(Dn.vals, 2, mean)
}

