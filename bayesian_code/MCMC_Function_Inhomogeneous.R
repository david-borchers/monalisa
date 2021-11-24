## Function to generate MCMC samples when fitting an inhomogeneous density model

## The arguments we need to provide are:
# * The data object, created so that it is a list containing the elements: 'encounter.data', 'trap.loc', and 'n.occasions'
# * A data frame ('pixel.info')  with three columns: the first column gives the x-coordinates for pixel centres in the region of interest, the second column gives the y-coordinates of the pixel centres, and the third gives the associated covariate value for each pixel centre. NOTE we assume that: (1) these pixel centres are evenly-spaced, (2) the region of interest is a square, so the number of pixels in the x- and y-directions is the square root of the number of rows in this data frame and (3) there are no regions in the region where animals cannot go (so the mask would just be a matrix of 1's)
# * The number of pixels in the x- and y-directions (so the region doesn't have to be a square)
# * The size of the super-population, M
# * A vector containing the starting values for lambda0, sigma and beta1 ('inits.vec'), where the order is: c(lambda0, sigma, beta1). Later, we calculate the starting values for 'log_coeff' and 'beta0' so that they 'make sense', so we won't provide them here.
# * The dmax value to use for the getLocalObjects() function
# * The number of iterations to run the MCMC for
# * The number of burn-in iterations we want to use

run.MCMC.inhom = function(data, pixel.info, x.pixels, y.pixels, M, inits.vec, dmax = 56, n.iter = 10000, n.burn = 1000) {
  ## Therefore, subsetting the data we'll use in our NIMBLE model:
  # Encounter data
  y = data$encounter.data
  # Checking if there are any rows of 0's -- if there are, returning an error because these capture histories are unobserveable
  if (any(apply(y,1,sum)==0)) stop("The data shouldn't include any all-0 capture histories")

  # Trap locations matrix
  traplocs = data$trap.loc
  # Number of traps
  n.trap = nrow(traplocs)
  # xlim, ylim (based on centres in 'pixel.info' being evenly-spaced)
  xlim = c(min(pixel.info[,1])-(0.5*(pixel.info[1,1] - pixel.info[2,1])), max(pixel.info[,1]+(0.5*(pixel.info[1,1] - pixel.info[2,1])))
  ylim = c(min(pixel.info[,2])-(0.5*(pixel.info[1,2] - pixel.info[2,2])), max(pixel.info[,2]+(0.5*(pixel.info[1,2] - pixel.info[2,2])))
  # Number of animals detected
  n.observed = nrow(y)

  # Number of pixels in the map region
  nPix = nrow(pixel.info)

  # Pixel centres we need
  pixel.centres = pixel.info[,1:2]
  # Matrix containing numbering for each pixel centre -- i.e. denotes number in which values in 'pixel.centres' occur in our survey area
  pixel.centres.order = matrix(1:2500, ncol=sqrt(nPix), nrow=sqrt(nPix), byrow=T)
  pixel.centres.order = pixel.centres.order[nrow(pixel.centres.order):1,]

  # Area of each pixel we are considering, calculated based on centres in 'pixel.info' being evenly spaced
  pixel.area = (pixel.info[1,1] - pixel.info[2,1])^2


  ## Data augmentation
  # Setting the size of the superpopulation
  M = M
  # Adding all-0 rows to our encounter data matrix (so that, in total, we have encounter data for M animals)
  y = rbind(y, matrix(0, nrow=M-n.observed, ncol=ncol(y)))
  # Vector of 0's and 1's corresponding to our encounter data matrix - 1 if a 'real' individual (a detected individual), 0 for an 'added' individual (i.e. an individual that we are considering, but was never detected at a trap)
  z = c(rep(1, n.observed), rep(0, M-n.observed))


  ## Setting the starting values for s (activity centres for each animal)
  # For observed animals, we want these initialised activity centres to be at or near the traps at which individuals were captured. So, we'll make the starting activity centre the 'mean' location for all of the traps at which they were caught.
  # First, generating 'random' activity centres for ALL individuals in our supopulation
  sst = cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
  # Now, for observed animals, making their initial activity centre the 'mean' of all of the traps at which they were detected
  for (j in 1:n.observed) {
    if (sum(y[j,])==0) {next}
    sst[j,1] = mean(traplocs[y[j,]>0,1])
    sst[j,2] = mean(traplocs[y[j,]>0,2])
  }
  # And now, rounding these starting values so that they correspond to pixel centres
  sst = round(sst)
  # Finding the corresponding pixel index for each row of 'sst' - i.e. for each activity centre we have generated, we are finding the index of the corresponding pixel (where pixel '1' will be in the first row of the 'pixel.centres' object, and so on)
  starting.pixel.indices = vector("numeric", M)
  for (j in 1:M) {
    index = which(sst[j,1]==pixel.centres[,1] & sst[j,2]==pixel.centres[,2])
    starting.pixel.indices[j] = index
  }
  # And now, making sst equal to the 'starting.pixel.indices' object we have just created - so what we end up providing to our NIMBLE model in the way of initial activity centres is the index of the pixel that each animals' initial activity centre falls into.
  sst = starting.pixel.indices
  # Finding sx and sy -- these are the row/column indices for each entry in 'pixel.centres.order' that is stored in 'sst'
  sx_sy_init = matrix(0, ncol=2, nrow=length(sst))
  for (i in (1:length(sst))) {
    sx_sy_init[i,] = which(pixel.centres.order==sst[i], arr.ind=T)
  }
  # Subsetting the initial values for sx
  sx_init = sx_sy_init[,1]
  # Initial values for sy
  sy_init = sx_sy_init[,2]

  # Covariate values we want to use in our model
  mona.densities = pixel.info[,3]

  # As we want to use dpoisLocal_normal() below, we  need to scale the trap coordinates using the scaleCoordsToHabitatGrid() function
  # Before this, need to label columns in pixel.centres as 'x' and 'y'
  colnames(pixel.centres) = c("x", "y")
  scaledtrapcoords = scaleCoordsToHabitatGrid(coordsData = traplocs, coordsHabitatGridCenter = pixel.centres)
  scaledtrapcoords = scaledtrapcoords$coordsDataScaled
  # Scaling the pixel centres, as well
  scaledpixelcentres = scaleCoordsToHabitatGrid(coordsData = pixel.centres, coordsHabitatGridCenter = pixel.centres)
  scaledpixelcentres = scaledpixelcentres$coordsDataScaled


  ## Defining the NIMBLE model
    code <- nimbleCode({
    # Priors -- same as for the homogeneous PP model (for the parameters that are common to both models)
    lambda0 ~ dgamma(0.001, 0.001)
    log_coeff ~ dunif(-10, 10)
    coeff  <- exp(log_coeff)
    sigma <- sqrt(1/(2*coeff))
    # (The parameters below are unique to the inhomogeneous PP model)
    beta0 ~ dunif(-25, 25)
    beta1 ~ dunif(-10, 10)

    # Specifying prior probabilities for each pixel
    DPix[1:nPix] <- exp(beta0 + beta1*(mona.densities[1:nPix]))
    mu[1:nPix] <- DPix[1:nPix] * pixel.area
    probs[1:nPix] <- mu[1:nPix]/EN

    EN <- sum(mu[1:nPix])  # Expected value of N, E(N)
    psi <- EN/M  # Data augmentation parameter

    for (k in 1:M) {
      z[k] ~ dbern(psi) # Whether or not the ith animal exists
      sx[k] ~ dunif(1, 51) # Prior for column of matrix that represents where in 'pixel.centres.order' the activity centre lies
      sy[k] ~ dunif(1, 51) # Prior for row of matrix that represents where in 'pixel.centres.order' the activity centre lies
      ind_x[k] <- trunc(sx[k]) # Finding the row of 'pixel.centres.order' that contains the sampled activity centre
      ind_y[k] <- trunc(sy[k]) # Finding the column of 'pixel.centres.order' that contains the sampled activity centre
      s[k] <- pixel.centres.order[ind_x[k], ind_y[k]] # Finding the corresponding entry in 'pixel.centres.order', which represents the index on the pixel that contains the sampled activity centre
      ones[k] ~ dbern(probs[s[k]]) # Adds to the likelihood the probability of our activity centre falling into the given pixel centre

      # Likelihood, following Wolverine NIMBLE example found at: https://nimble-dev.github.io/nimbleSCR/wolverine_example.html
      y[k, 1:nMaxDetectors] ~ dpoisLocal_normal(detNums = nbDetections[k],
                                                  detIndices = yDets[k, 1:nMaxDetectors],
                                                  lambda = lambda0,
                                                  s = scaledpixelcentres[s[k],1:2], # s[k] is the index of the pixel centre into which the given activity centre falls -- so should index the row in 'scaledpixelcentres' that contains the chosen activity centre for the sampled animal
                                                  sigma = sigma,
                                                  trapCoords = scaledtrapcoords[1:n.trap, 1:2],
                                                  localTrapsIndices = detectorIndex[1:n.cells,1:maxNBDets],
                                                  localTrapsNum = nDetectors[1:n.cells],
                                                  resizeFactor = ResizeFactor,
                                                  habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                  indicator = z[k])
    }
    N <- sum(z[1:M])  # Realised value of N
    D <- (N/((xlim[2] - xlim[1]) * (ylim[2] - ylim[1]))) * 10000
  }
  )

  ## Values that we want to provide to our NIMBLE model
  # Data
  data <- list(scaledtrapcoords = scaledtrapcoords,
               scaledpixelcentres = scaledpixelcentres,
               mona.densities = mona.densities,
               ones = rep(1, M),
               pixel.centres.order = pixel.centres.order)
  # Constants
  constants <- list(nPix=nPix, M=M, n.trap=n.trap, xlim=xlim, ylim=ylim,
                    pixel.area=pixel.area)
  # Initial values
  inits = list (lambda0=inits.vec[1], log_coeff=log(1/(2*(inits.vec[2])^2)), z=z, beta0=log(n.observed/(nPix*pixel.area*10000)), beta1=inits.vec[3], s=sst, sx=sx_init, sy=sy_init)

  ## More constants that we want to provide to the NIMBLE model
  # And since we assume we have no unsuitable habitat, we are defining the 'habitatMask' argument as a matrix full of 1's.
  DetectorIndex <- getLocalObjects(habitatMask = matrix(1, ncol=50, nrow=50), coords = scaledtrapcoords, dmax = 56,resizeFactor = 1, plot.check=FALSE)
  # Generating the values of more constants that we want to provide to our NIMBLE model
  constants$y.maxDet <- dim(DetectorIndex$habitatGrid)[1]
  constants$x.maxDet <- dim(DetectorIndex$habitatGrid)[2]
  constants$ResizeFactor <- DetectorIndex$resizeFactor
  constants$n.cells <- nPix
  constants$maxNBDets <- DetectorIndex$numLocalIndicesMax
  data$detectorIndex <- DetectorIndex$localIndices
  data$nDetectors <- DetectorIndex$numLocalIndices
  data$habitatIDDet <- DetectorIndex$habitatGrid

  ## More data that we want to provide to the NIMBLE model
  # Generating a sparse representation of encounter data matrix
  ySparse <- getSparseY(x = y)
  data$y <- ySparse$y[,,1]
  data$yDets <- ySparse$detIndices[,,1]
  data$nbDetections <- ySparse$detNums[,1]
  constants$nMaxDetectors <- ySparse$maxDetNums

  ## Doing the clever NIMBLE stuff
  Rmodel <- nimbleModel(code, constants, data, inits, dimensions = list(pixel.centres.order = c(50,50)))
  # AF slice sampling for beta0 and beta1
  conf <- configureMCMC(Rmodel, monitors = c("lambda0", "sigma", "N", "D", "beta0", "beta1", "DPix"), print = FALSE)
  conf$removeSampler(c("beta0","beta1"))
  conf$addSampler(target = c("beta0","beta1"),
                  type = 'AF_slice',
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  samples <- runMCMC(Cmcmc, niter = n.iter+n.burn, progressBar=TRUE)

  # Returning the MCMC samples
  samples
}
