## Here, 'data' refers to a data object created following the method in 'Gathering_Simulated_Data.R' (can also see examples in 'Figure6_Figure7_Code.R'). Note that due to the way in which we subset xlim and ylim below, the lower value must occur first in these vectors.

## 'M' is the size of the superpopulation and 'lambda0.prior' is the prior chosen for the lambda0 parameter.

## 'n.iter' is the number of MCMC iterations you want to run, 'n.adapt' is the number of iterations for adaptation you want (different from burn-in), 'n.burn' is the number of burn-in iterations that we want to use.

## 'lambda0.start' gives the initial value for lambda0 for the MCMC sampling.

## 'log_coeff.start' gives the initial value for 'log_coeff' (as we are using a log-uniform prior for the 'coeff' parameter).

## 'thin' sets the thinning parameter for the MCMC.

run.MCMC <- function(data, M, coeff.prior, parameters, n.iter=1000, n.adapt = 1000, n.burn = 100, lambda0.start = runif(1, 0, 50), log_coeff.start=runif(1, 0, 1), thin = 1) {

  ## Encounter matrix, no rows of zeroes - need to remove rows of zeroes manually!
  y <- data$encounter.data
  ## Removing rows of zeroes in encounter data
  y <- y[apply(y, 1, sum)>0, ]
  ## Trap locations matrix
  traplocs <- data$trap.loc
  ## Number of traps
  trap.no <- nrow(traplocs)
  ## Number of occasions over which data was collected
  n.occ <- data$n.occasions
  ## xlim
  xlim <- data$xlim
  ## ylim
  ylim <- data$ylim
  ## Number of animals detected
  pop.size <- nrow(y)

  ## Data augmentation
  # Total number of possible animals that could be living in the area - this INCLUDES the number of animals that
  # WERE detected/the number of animals from the original population size that weren't detected!
  M <- M
  # Adding all-0 rows (so that, in total, we have encounter data for 500 animals, out of which
  # M-pop.size were unobserved at the traps, but could possibly have existed in the area)
  y <- rbind(y, matrix(0, nrow=M-pop.size, ncol=ncol(y)))
  # Vector of 0's and 1's - 1 if a 'real' individual (a detected individual), 0 for an 'added' individual - i.e.
  # an individual that we are considering, but was never detected at a trap
  z <- c(rep(1, pop.size), rep(0, M-pop.size))

  ## Starting values for s (activity centres) - want activity centres at or near the traps at which individuals were captured - make activity
  ## centre the 'mean' for all the traps at which the individual was detected
  # NOTE - activity center created for every individual, real or potential
  ## Generating 'random' activity centres for ALL individuals (real/potential), and later will make activity centre
  ## the 'mean' of all the traps at which individuals were detected, for detected individuals ONLY
  sst <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
  # NOTE - need a 'next' condition below because some rows DO contain all zeroes

  for (i in 1:pop.size) {
    if (sum(y[i,])==0) {next}
    sst[i,1] <- mean(traplocs[y[i,]>0,1])
    sst[i,2] <- mean(traplocs[y[i,]>0,2])
  }

  ## NIMBLE model
  # We could monitor multiple things here: lambda0, psi, log_coeff, coeff, sigma, z-values and s-values (activity centres) -
  # MCMC will find posterior distirbutions for these values, so that will be able to determine the likely values of these (i.e.
  # the values of these under which the data is most likely and which we expected to be the most likely!) NOTE that our previous
  # beliefs for the values will only come into play if we have uniform priors (which we do tend to have - so here, we are
  # mainly looking for parameter values that maximise the likelihood function).
  # Generally, we enter the argument: parameters = c("lambda0", "coeff", "sigma", "N", "D", "z", "s") - these are the parameters we tend to monitor!

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
