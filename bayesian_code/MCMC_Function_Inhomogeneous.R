## Function to generate MCMC samples when fitting an inhomogeneous density model

## The arguments we need to provide are:
# * The data object, created so that it is a list containing the elements: 'encounter.data', 'trap.loc', 'xlim', 'ylim' (where xlim and ylim are in increasing order) and 'n.occasions'
# * The number of pixels in the map region
# * The pixel area
# * The size of the super-population, M
# * The column from 'mona_df' that contains the covariate values we want to use
# * Starting values for lambda0 and sigma
# * The dmax value to use for the getLocalObjects() function
# * The number of iterations to run the MCMC for
# * The number of burn-in iterations we want to use
# * The length of the adaptation interval that we want to use
# * The number of pixels in the x- and y-directions

run.MCMC.inhom = function(data, nPix = 2500, pixel.area = 1, M, covariate, lambda0.start = 0.5, sigma.start = 4, dmax = 56, n.iter = 10000, n.burn = 1000) {
  # Therefore, subsetting the values we'll use in our NIMBLE model:
  ## Encounter data
  y = data$encounter.data
  # Removing rows of 0's in encounter data
  y = y[apply(y, 1, sum)>0, ]

  ## Trap locations matrix
  traplocs = data$trap.loc
  ## Number of traps
  trap.no = nrow(traplocs)
  ## Number of occassions for which data was collected
  n.occ = data$n.occasions
  ## xlim
  xlim = data$xlim
  ## ylim
  ylim = data$ylim
  ## Number of animals detected
  pop.size = nrow(y)

  ## Number of pixels in the map region
  nPix = nPix

  ## Pixel centres in the map region
  library("spatstat")
  # Function
  centres = function(xrange, yrange, x.pixels, y.pixels) {
    window.2 = owin(xrange=xrange, yrange=yrange)
    points = gridcentres(window.2, x.pixels, y.pixels)
    centres = as.matrix(cbind(points$x, points$y))
    centres
  }
  # Pixel centres we need
  pixel.centres = centres(xrange=c(0.5,50.5), yrange=c(0.5,50.5), x.pixels=50, y.pixels=50)
  # Matrix containing numbering for each pixel centre -- i.e. denotes number in which values in 'pixel.centres' occur in our survey area
  pixel.centres.order = matrix(1:2500, ncol=50, nrow=50, byrow=T)
  pixel.centres.order = pixel.centres.order[nrow(pixel.centres.order):1,]

  ## Area of each pixel we are considering
  pixel.area = pixel.area


  ### Data augmentation
  # Setting the size of the superpopulation - the total number of possible animals that could possibly be living in the area (this consists of animals
  # that were detected)
  M = M # Note that the true N (true population size) is 7451
  # Adding all-0 rows to our encounter data matrix (so that, in total, we have encounter data for M animals, out of which (M-pop.size) were unobserved
  # at the traps, but could possibly have existed in the area)
  y = rbind(y, matrix(0, nrow=M-pop.size, ncol=ncol(y)))
  # Vector of 0's and 1's corresponding to our encounter data matrix - 1 if a 'real' individual (a detected individual), 0 for an 'added' individual (i.e.
  # an individual that we are considering, but was never detected at a trap)
  z = c(rep(1, pop.size), rep(0, M-pop.size))


  ### Setting the starting values for s (activity centres for each animal)
  ## We want these initialised activity centres to be at or near the traps at which individuals were captured - so for observed animals, we'll make the
  ## starting activity centre the 'mean' location for all of the traps at which they were caught.
  ## NOTE - that we are initialising the activity centers for all animals in our superpopulation (observed or unobserved)
  # First, generating 'random' activity centres for ALL individuals in our supopulation
  sst = cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
  # Now, for observed animals, making their initial activity centre the 'mean' of all of the traps at which they were detected
  # NOTE - need a 'next' condition now, because some rows of y DO contain all zeroes (as we are now considering the superpopoulation, which contains unobserved
  # animals)
  for (j in 1:pop.size) {
    if (sum(y[j,])==0) {next}
    sst[j,1] = mean(traplocs[y[j,]>0,1]) # Setting the x-coord of the initial activity centre to be the mean of the x-coord of all traps at which the animal was detected
    sst[j,2] = mean(traplocs[y[j,]>0,2]) # Doing the same thing for the y-coord of the initial activity centre (and leaving the initial activity centre to be the randomised point above if the animal was not detected at all)
  }
  # And now, rounding these starting values for the activity centres so that they correspond to pixel centres - this means that we now know which pixel
  # centres each initialised activity centre falls into -- this is what we'll provide to our NIMBLE model (i.e. we'll be dealing with pixel centres
  # when considering the animals' activity centres).
  sst = round(sst)
  ## Finding the corresponding pixel index for each row of 'sst' - i.e. for each activity centre we have generated, we are finding the number/index of the
  ## corresponding pixel
  # Vector that we can store these indices in
  starting.pixel.indices = vector("numeric", M)
  for (j in 1:M) {
    # Generating the number/index of the pixel the activity centre falls into, based on the order of the pixel centres in our 'pixel.centres' object above
    # (generated using the centres() function)
    index = which(sst[j,1]==pixel.centres[,1] & sst[j,2]==pixel.centres[,2])
    # Storing this index in the 'starting.pixel.indices' vector
    starting.pixel.indices[j] = index
  }
  ## AND NOW, making sst equal to the 'starting.pixel.indices' object we have just created - so what we end up providing to our NIMBLE model in the way of
  ## initial activity centres is the index of the pixel that each animals' initial activity centre falls into.
  sst = starting.pixel.indices
  ## Using the values of sst to get the starting values for 'sx' and 'sy' in our model below -- these are the row and column number (respectively) for the corresponding
  ## entry in 'pixel.centres.order' (rather than x-/y-coordinates). That is, 'sx' and 'sy' give the row/column indices for the entry in 'pixel.centres.order' that represents the pixels in sst (based on the manner in which our survey area has been numbered)
  # Matrix containing the starting values for 'sx' and 'sy'
  sx_sy_init = matrix(0, ncol=2, nrow=length(sst))
  for (i in (1:length(sst))) {
    sx_sy_init[i,] = which(pixel.centres.order==sst[i], arr.ind=T)
  }
  # Subsetting the initial values for sx
  sx_init = sx_sy_init[,1]
  # Initial values for sy
  sy_init = sx_sy_init[,2]

  # Covariate values we provide
  mona.densities = covariate

  ## As we want to use dpoisLocal_normal() below, we apparently need to scale the trap coordinates to the 'habitat gird' using the scaleCoordsToHabitatGrid() function
  # Before this, need to label columns in pixel.centres as 'x' and 'y'
  colnames(pixel.centres) = c("x", "y")
  scaledtrapcoords = scaleCoordsToHabitatGrid(coordsData = traplocs, coordsHabitatGridCenter = pixel.centres)
  scaledtrapcoords = scaledtrapcoords$coordsDataScaled

  ## Now, we want to use NIMBLE to run MCMC - we are expecting to see correlation between beta0 and beta1 here
  code <- nimbleCode({
    ## Priors -- same as for the homogeneous PP model (for the parameters that are common to both models)
    lambda0 ~ dgamma(0.001, 0.001)
    log_coeff ~ dunif(-10, 10)
    coeff  <- exp(log_coeff)
    sigma <- sqrt(1/(2*coeff))
    # (The parameters below are unique to the inhomogeneous PP model)
    beta0 ~ dunif(-25, 25)
    beta1 ~ dunif(-10, 10)

    ## Specifying prior probabilities for each pixel
    DPix[1:nPix] <- exp(beta0 + beta1*(mona.densities[1:nPix]))
    mu[1:nPix] <- DPix * pixel.area
    probs[1:nPix] <- mu[1:nPix]/EN

    EN <- sum(mu[1:nPix])  # Expected value of N, E(N)
    psi <- EN/M  # Data augmentation parameter

    for (k in 1:M) {
      z[k] ~ dbern(psi) # Whether or not the ith animal exists
      sx[k] ~ dunif(1, 51) # Prior for column of matrix that pixel centre for activity centre lies
      sy[k] ~ dunif(1, 51) # Prior for row of matrix that pixel centre for activity centre lies
      ind_x[k] <- trunc(sx[k]) # Finding value of the row in which the x-coordinate of pixel centre for sampled activity centre lies
      ind_y[k] <- trunc(sy[k]) # Finding value of the column in which the y-coordinate of pixel centre for sampled activity centre lies
      s[k] <- pixel.centres.order[ind_x[k], ind_y[k]] # Specifying which pixel centre the given activity centre falls into (by the number of the pixel centre)
      ones[k] ~ dbern(probs[s[k]]) # Adds to the likelihood the probability of our activity centre falling into the given pixel centre

      ## Likelihood, following Wolverine NIMBLE example found at: https://nimble-dev.github.io/nimbleSCR/wolverine_example.html
      y[k, 1:nMaxDetectors] ~ dpoisLocal_normal(detNums = nbDetections[k],
                                                  detIndices = yDets[k, 1:nMaxDetectors],
                                                  lambda = lambda0,
                                                  s = pixel.centres[s[k],1:2], # s[k] is number of pixel centre that the activity centre falls into -- so should index the row in 'pixel.centres' that contains the chosen activity centre for the sampled animal
                                                  sigma = sigma,
                                                  trapCoords = scaledtrapcoords[1:trap.no, 1:2],
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

  ### Values that we want to provide to our NIMBLE model
  ## Data
  data <- list(scaledtrapcoords = scaledtrapcoords,
               pixel.centres = pixel.centres,
               mona.densities = mona.densities,
               ones = rep(1, M),
               pixel.centres.order = pixel.centres.order)
  ## Constants
  constants <- list(nPix=nPix, n.occ=n.occ, M=M, trap.no=trap.no, xlim=xlim, ylim=ylim,
                    pixel.area=pixel.area)
  ## Initial values.
  inits = list (lambda0=lambda0.start, log_coeff=log(1/(2*sigma.start^2)), z=z, beta0=log(pop.size/(nPix*pixel.area*10000)), beta1=0, s=sst, sx=sx_init, sy=sy_init)
  # True value of p0 is given by: 1-exp(-0.69)=0.4984 (4 sf), as we know that the true value of lambda0 is 0.69 for our simulation. So, starting lambda0 at 0.5 means we are starting near the true value.
  # And the true value of sigma is 2, so once again starting at the true value just to see if things will work!

  ## More constants that we want to provide to the NIMBLE model
  # Based on the Wolverine example, we are setting the dmax value below to 56, which seems to be okay.
  # dmax restricts calculations of detection probabilities to detectors within a 56 unit radius where detections are possible (i.e. within a 56 unit radius from the activity centre)
  # The idea of using dmax is that it potentially reduces the detectors for which we calculate the detection function later
  # The aim is to have dmax as small as possible to speed up computation, while being large enough to include detectors at which animals were detected
  # Here, we choose dmax of 56 as means that for all animals, will be considering all detectors (so aren't really speeding up computation here) -- e.g. if activity centre at (0,0), we are saying could still be detected at the furthest away detector!
  # We are using simulated animals here, there isn't a real sense of what detectors would be 'realistic' for animals to be caught at, so seems sensible to set dmax like this
  # And since we have no unsuitable habitat, we are defining the 'habitatMask' argument as a matrix full of 1's.
  DetectorIndex <- getLocalObjects(habitatMask = matrix(1, ncol=50, nrow=50), coords = scaledtrapcoords, dmax = 56,resizeFactor = 1)
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
  # Reinitialising constants
  #Cmodel$sigma = sigma.start
  #Cmodel$beta1 = 1
  #Cmodel$beta0 = 0
  #Cmodel$lambda0 = lambda0.start

  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  samples <- runMCMC(Cmcmc, niter = n.iter+n.burn, progressBar=TRUE)#, nburnin = n.burn)#, samplesAsCodaMCMC=FALSE, nburnin = n.burn))

  # Returning the MCMC samples
  samples
}
