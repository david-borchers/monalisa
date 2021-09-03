## Here, 'xlim' and 'ylim' give the range of x- and y-coordinates for the map area (see the comment preceding the 'centres()' function 
## in 'RUDMaps_Functions.R' for further explanation on this). Also, 'results' refers to a set of MCMC samples generated using run.MCMC() 
## and M is the size of the superpopulation. 

no.movement.density.vector <- function(xlim, ylim, results, M) {
  ## Generating density values for each pixel, storing as a vector -- note that the resulting density values do not take the movement of 
  ## animals into account! (RUD maps take this movement into account, RACD maps do not). 
  
  ## Points at which local density will be estimated
  xg <- seq(xlim[1], xlim[2], by=1)
  yg <- seq(ylim[1], ylim[2], by=1)
  
  ## Extracting z-values
  # Names of variables that have been monitored
  names <- names(results[1,])
  # Extracting "z" values from MCMC results
  Z <- results[,grep("z", names)]
  # Logical vector - TRUE if z=1, FALSE if z=0 - for EACH z-value, each iteration
  logical <- (Z == 1)
  
  ## Extracting activity centres
  # Extracting "s" values from MCMC results (i.e. extracting sampled activity centres)
  S <- results[,grep("s[^i]", names)]
  # x-coordinates of all activity centres
  Sx <- S[,1:M]
  # y-coordinates of all activity centres
  Sy <- S[,-(1:M)]
  # Extracting relevant activity centres - x co-ordinates and y co-ordinates where z=1
  Sxout <- Sx[logical == 1]
  Syout <- Sy[logical == 1]
  
  ## Associating each activity centre from above with the pixel it falls into
  Sxout <- cut(Sxout, breaks=xg, include.lowest=TRUE)
  Syout <- cut(Syout, breaks=yg, include.lowest=TRUE)
  
  ## Tallying up how many activity centres are in each pixel
  Dn <- table(Sxout, Syout)
  
  ## Dividing the number of times an activity centre fell in each pixel by the number of iterations 
  table <- Dn/nrow(Z)
  density.vector <- as.vector(table)
}

