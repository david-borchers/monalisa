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

run.MCMC.inhom = function(data, nPix = 2500, pixel.area = 1, M, mona.column, lambda0.start = 0.5, sigma.start = 4, dmax = 56, n.iter = 10000, n.burn = 1000, x.pixels = 50, y.pixels = 50) {

  ## Subsetting the data values we will use in the NIMBLE model:
  # Encounter data (with all-0 capture histories removed)
  y = data$encounter.data
  y = y[apply(y, 1, sum)>0, ]
  # Trap locations
  traplocs = data$trap.loc
  # Number of traps
  trap.no = nrow(traplocs)
  # Number of sampling occasions
  n.occ = data$n.occasions
  # xlim
  xlim = data$xlim
  # ylim
  ylim = data$ylim
  # Number of animals detected
  pop.size = nrow(y)


  ## Generating objects relating to the pixels in our map region that we'll use in the NIMBLE model:
  # Number of pixels in the map region
  nPix = nPix
  # Area of each pixel
  pixel.area = pixel.area
  # Generating matrix of pixel centres in the map region
  # Sourcing 'RUDMaps_Functions.R' for the centres() function
  source("RUDMaps_Functions.R")
  pixel.centres = centres(xrange=c(xlim[1], xlim[2]), yrange=c(ylim[1], ylim[2]), x.pixels = x.pixels, y.pixels = y.pixels)
  # Matrix containing numbering for each pixel centre -- i.e. denotes number in which values in 'pixel.centres' occur in our survey area
  pixel.centres.order = matrix(1:nPix, ncol=x.pixels, nrow=y.pixels, byrow=T)
  pixel.centres.order = pixel.centres.order[nrow(pixel.centres.order):1,]

  ## Data augmentation
  # Setting the size of the superpopulation
  M = M
  # Adding M-(pop.size) all-0 rows to the encounter matrix
  y = rbind(y, matrix(0, nrow=M-pop.size, ncol=ncol(y)))
  # Creating the z-vector (where we have 1's for each detected animal, 0's otherwise)
  z = c(rep(1, pop.size), rep(0, M-pop.size))


  ## Starting values for s (denotes index of pixel containing activity centre for each animal)
  # Generating 'random' activity centres for eall of the individuals in the sueprpopulation
  sst = cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
  # For observed animals, making their initial activity centre the 'mean' of all the traps at which they were observed
  for (j in 1:pop.size) {
    if (sum(y[j,])==0) {next}
    sst[j,1] = mean(traplocs[y[j,]>0,1]) # Setting the x-coord of the initial activity centre to be the mean of the x-coord of all traps at which the animal was detected
    sst[j,2] = mean(traplocs[y[j,]>0,2]) # Doing the same thing for the y-coord of the initial activity centre (and leaving the initial activity centre to be the randomised point above if the animal was not detected at all)
  }

  # Rounding the values in 'sst', so each row corresponds to a pixel centre
  sst = round(sst)
  # Finding the corresponding index for the pixels in each row of 'sst'
  starting.pixel.indices = vector("numeric", M)
  for (j in 1:M) {
    # Generating the number/index of the pixel the activity centre falls into, based on the order of the pixel centres in our 'pixel.centres' object above
    index = which(sst[j,1]==pixel.centres[,1] & sst[j,2]==pixel.centres[,2])
    # Storing this index in the 'starting.pixel.indices' vector
    starting.pixel.indices[j] = index
  }

  # Making 'sst' equal to the 'starting.pixel.indices' vector. This is what we will provide to the NIMBLE model (in the way of initial values for the activity centres)
  sst = starting.pixel.indices

  # Using the values of sst to generate starting values for 'sx' and 'sy' in the NIMBLE model. 'sx' and 'sy' are the row and column numbers (respectively) in 'pixel.centres.order' for each entry of sst.
  sx_sy_init = matrix(0, ncol=2, nrow=length(sst))
  for (i in (1:length(sst))) {
    sx_sy_init[i,] = which(pixel.centres.order==sst[i], arr.ind=T)
  }
  # Subsetting the initial values for sx
  sx_init = sx_sy_init[,1]
  # Initial values for sy
  sy_init = sx_sy_init[,2]


  ## Finding the covariate values for each pixel
  load("../output/mona_inputs.RData")
  # We are working with the 'Dblur' covariate for ch7b, extracting the values for this covariate (and the corresponding pixel centres)
  mona.densities = mona_df[, c("x", "y", mona.column)]
  # Each row in this data frame is copied 3 times -- removing the duplicated rows
  sequence = seq(1, 7500, by=3)
  mona.densities = mona.densities[sequence,]
  # Reordering this data frame, so the order of the pixel centres here matches the order of the centres in 'pixel.centres' (both currently have reverse ordering to each other wrt the y-values)
  split = split(mona.densities, mona.densities$y)
  mona.densities = do.call("rbind", split)
  # Extracting the vector of covariate values only, as this is what we want to use with NIMBLE
  mona.densities = mona.densities[, mona.column]


  ## Specifying the NIMBLE model
  code <- nimbleCode({
    ## Priors -- same as for the homogeneous PP model (for the parameters that are common to both models)
    lambda0 ~ dgamma(0.001, 0.001)
    log_coeff ~ dunif(-10, 10)
    coeff  <- exp(log_coeff)
    sigma <- sqrt(1/(2*coeff))
    # (The parameters below are unique to the inhomogeneous PP model)
    beta0 ~ dunif(-10, 10)
    beta1 ~ dunif(-10, 10)

    ## Specifying prior probabilities for each pixel
    mu[1:nPix] <- exp(beta0 + beta1*(log(mona.densities[1:nPix]))) * pixel.area
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
      y[k, 1:nMaxDetectors] ~ dpoisLocal_normal_2(detNums = nbDetections[k],
                                                  detIndices = yDets[k, 1:nMaxDetectors],
                                                  lambda0 = lambda0,
                                                  s = pixel.centres[s[k],1:2], # s[k] is number of pixel centre that the activity centre falls into -- so should index the row in 'pixel.centres' that contains the chosen activity centre for the sampled animal
                                                  sigma = sigma,
                                                  trapCoords = traplocs[1:trap.no, 1:2],
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
  data <- list(traplocs = traplocs, pixel.centres = pixel.centres, mona.densities = mona.densities, ones = rep(1, M), pixel.centres.order = pixel.centres.order)
  # Constants
  constants <- list(nPix=nPix, n.occ=n.occ, M=M, trap.no=trap.no, xlim=xlim, ylim=ylim, pixel.area=pixel.area)
  # Initial values (note the true value of lambda0 is 0.69, and the true value of sigma is 2)
  inits = list (lambda0=lambda0.start, sigma=sigma.start, z=z, beta0=0, beta1=1, s=sst, sx=sx_init, sy=sy_init)

  # More constants to provide to the NIMBLE model, based on how the nimbleSCR package works
  DetectorIndex <- getLocalObjects(habitatMask = matrix(1, ncol=x.pixels, nrow=y.pixels), coords = traplocs, dmax = dmax,resizeFactor = 1)
  constants$y.maxDet <- dim(DetectorIndex$habitatGrid)[1]
  constants$x.maxDet <- dim(DetectorIndex$habitatGrid)[2]
  constants$ResizeFactor <- DetectorIndex$resizeFactor
  constants$n.cells <- nPix
  constants$maxNBDets <- DetectorIndex$numLocalIndicesMax
  data$detectorIndex <- DetectorIndex$localIndices
  data$nDetectors <- DetectorIndex$numLocalIndices
  data$habitatIDDet <- DetectorIndex$habitatGrid

  # Generating a sparse representation of the encounter matrix, and amending the data/constants based on the way nimbleSCR works
  ySparse <- getSparseY(x = y)
  data$y <- ySparse$y[,,1]
  data$yDets <- ySparse$detIndices[,,1]
  data$nbDetections <- ySparse$detNums[,1]
  constants$nMaxDetectors <- ySparse$maxDetNums


  ## Running the MCMC!
  Rmodel <- nimbleModel(code, constants, data, inits, dimensions = list(pixel.centres.order = c(50,50)))
  # AF slice sampling for beta0 and beta1
  conf <- configureMCMC(Rmodel, monitors = c("lambda0", "sigma", "N", "D", "z", "s", "beta0", "beta1", "sx", "sy"), print = FALSE)#, control = list(adaptInterval = n.adapt))
  conf$removeSampler(c("beta0","beta1"))
  conf$addSampler(target = c("beta0","beta1"), type = 'AF_slice', control = list(adaptScaleOnly = TRUE), silent = TRUE)
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  system.time(samples <- runMCMC(Cmcmc, niter = n.iter+n.burn, nburnin = n.burn, progressBar=TRUE, samplesAsCodaMCMC=FALSE))

  # Returning the MCMC samples
  samples
}
