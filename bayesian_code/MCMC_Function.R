## Function to generate MCMC samples when fitting a homogeneous density model

## The arguments we need to provide (in order) are:
# * 'data':  a data object created in a similar manner to the data objects from the beginning of 'Figure9.R'. This is a list that should include elements labelled: 'encounter.data' (capture history matrix), 'trap.loc' (trap coordinates), 'xlim' (x-range of pixel centres), 'ylim' (y-range of pixel centres) and 'n.occasions' (number of sampling occasions)
# * 'M': the size of the superpopulation
# * 'n.iter': the number of MCMC iterations you want to run
# * 'n.adapt': the number of adaptation iterations you want (different from burn-in)
# * 'n.burn': the number of burn-in iterations that we want to use (these will be discarded for us)
# * 'lambda0.start': the initial value for lambda0 for the MCMC sampling
# * 'log_coeff.start': the initial value for 'log_coeff' for the MCMC sampling (see NIMBLE model below for an idea of what parameter this is)
# * 'thin': the value of the thinning parameter for the MCMC
# * 'parameters': the vector of parameters we want to monitor


run.MCMC <- function(data, M, n.iter=1000, n.adapt = 1000, n.burn = 100, lambda0.start = runif(1, 0, 50), log_coeff.start=runif(1, 0, 1), thin = 1, parameters) {
  ## Subsetting the data we'll use in the NIMBLE model:
  # Encounter matrix, no rows of zeroes
  y <- data$encounter.data
  # Removing rows of zeroes in encounter data
  y <- y[apply(y, 1, sum)>0, ]
  # Trap locations matrix
  traplocs <- data$trap.loc
  # Number of traps
  trap.no <- nrow(traplocs)
  # Number of occasions over which data was collected
  n.occ <- data$n.occasions
  # xlim
  xlim <- data$xlim
  # ylim
  ylim <- data$ylim
  # Number of animals detected
  pop.size <- nrow(y)

  ## Data augmentation
  # Setting the size of the super-population
  M <- M
  # Adding all-0 rows (so that, in total, we have encounter data for M animals)
  y <- rbind(y, matrix(0, nrow=M-pop.size, ncol=ncol(y)))
  # Vector of 0's and 1's corresponding to our encounter data matrix - 1 if a 'real' individual (a detected individual), 0 for an 'added' individual (i.e. an individual that we are considering, but was never detected at a trap)
  z <- c(rep(1, pop.size), rep(0, M-pop.size))

  ## Setting the starting values for s (activity centres for each animal)
  # For observed animals, we want these initialised activity centres to be at or near the traps at which individuals were captured. So, we'll make the starting activity centre the 'mean' location for all of the traps at which they were caught.
  # First, generating 'random' activity centres for ALL individuals in our supopulation
  sst <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
 # Now, for observed animals, making their initial activity centre the 'mean' of all of the traps at which they were detected
  for (i in 1:pop.size) {
    if (sum(y[i,])==0) {next}
    sst[i,1] <- mean(traplocs[y[i,]>0,1])
    sst[i,2] <- mean(traplocs[y[i,]>0,2])
  }

  ## NIMBLE model
  # We could monitor multiple things here: lambda0, psi, log_coeff, coeff, sigma, z-values and s-values (activity centres). Generally, we enter the argument: parameters = c("lambda0", "coeff", "sigma", "N", "D", "z", "s") - these are the parameters we tend to monitor!
  x <- nimbleCode({
    lambda0~dgamma(0.001, 0.001)
    psi ~ dunif(0,1)
    log_coeff ~ dunif(-10,10)
    coeff <- exp(log_coeff)
    sigma <- sqrt(1/(2*coeff))

    for (i in 1:M) {
      z[i] ~ dbern(psi)
      s[i,1] ~ dunif(xlim[1], xlim[2])
      s[i,2] ~ dunif(ylim[1], ylim[2])
      for (j in 1:trap.no) {
        d[i,j] <- sqrt((s[i,1] - traplocs[j,1])^2 + (s[i,2] - traplocs[j,2])^2)
        lambda[i,j] <- z[i] * lambda0 * exp(-coeff * d[i,j] * d[i,j])
        y[i,j] ~ dpois(lambda[i,j]*n.occ)
      }
    }

    N <-  sum(z[1:M])
    D <-  (N/((xlim[2] - xlim[1]) * (ylim[2] - ylim[1]))) * 10000
  })

  ## Data to provide to NIMBLE
  nim.data <- list(y=y, traplocs=traplocs)

  ## Constants to provide to NIMBLE
  constants <- list(n.occ=n.occ, M=M, trap.no=trap.no, xlim=xlim, ylim=ylim)

  ## Initial values for the MCMC
  inits <- list(lambda0=lambda0.start, log_coeff=log_coeff.start, s=sst, z=z)

  ## Parameters to monitor
  parameters <- parameters

  ## Running NIMBLE
  Rmodel <- nimbleModel(code=x, constants=constants, data=nim.data, inits=inits)
  conf <- configureMCMC(Rmodel, monitors=parameters, control = list(adaptInterval = n.adapt))
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

  ## Running the MCMC, generating and saving the final results
  results <- runMCMC(Cmcmc, niter=n.iter+n.burn, nburnin=n.burn, progressBar=TRUE, samplesAsCodaMCMC=T)
}
