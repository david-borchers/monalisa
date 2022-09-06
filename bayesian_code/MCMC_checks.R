## Code to check the validity of the MCMC samples we produce when fitting SCR models that assume the state process can be represented by an inhomogeneous Poisson process

## Libraries we need
library(secr)
## Objects we need
load("../output/revision/mona-inputs.RData")
## Functions we need
source("Functions.R")

## ---------------------------------------------------------------------------------------
# Checking the results for Figure 5
## ---------------------------------------------------------------------------------------

##### Loading in the MCMC samples #####

load("MCMC_Results/Figure5/InhomPP_18occ.RData")
load("MCMC_Results/Figure5/InhomPP_52occ.RData")
load("MCMC_Results/Figure5/InhomPP_111occ.RData")

## Discarding buring in (are discarding 1000 iterations as burn-in, see bayesian_code/Figure5_code.R to look at trace plots)
inhom.results.18occ <- inhom.results.18occ[-c(1:1000),]
inhom.results.52occ <- inhom.results.52occ[-c(1:1000),]
inhom.results.111occ <- inhom.results.111occ[-c(1:1000),]

inhom.results.18occ <- inhom.results.18occ[-c(1:2500),]
inhom.results.52occ <- inhom.results.52occ[-c(1:2500),]
inhom.results.111occ <- inhom.results.111occ[-c(1:2500),]

##### Comparing MCMC results to results from secr.fit() #####

## Looking at the first plot in Row 2
check.inhom.mcmc(results=inhom.results.18occ, j=1, mask=mona_mask)

## Looking at the second plot in Row 2
check.inhom.mcmc(results=inhom.results.52occ, j=2, mask=mona_mask)

## Looking at the third plot in Row 2
check.inhom.mcmc(results=inhom.results.111occ, j=3, mask=mona_mask)

# There are some differences in both sets of results, in all instances. However, things generally look good. 

