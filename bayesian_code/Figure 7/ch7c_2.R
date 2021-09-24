## Creating ch7c RData object with the function have written

## First, sourcing in 'posthoc_extract_chs.R'
#load("../output/mona_raw_outputs_100sim.RData")
#source("../code/posthoc_extract_chs.R")
load("mona_raw_outputs_100sim.RData") # for NeSI
source("posthoc_extract_chs.R") # for NeSI

## ch7c
encounterdat.ch7c = matrix(0, nrow=nrow(ch7c[,1,]), ncol=ncol(ch7c[,1,]))
for (i in 1:20) {
  encounterdat.ch7c = encounterdat.ch7c + ch7c[,i,]
}
xlim = c(0.5, 50.5)
ylim = c(0.5, 50.5)
ch7c.traploc = attributes(ch7c)$traps
data.ch7c = list(encounter.data = encounterdat.ch7c, trap.loc = ch7c.traploc, xlim = xlim, ylim = ylim, n.occasions = 20)

## ---------------------------------------------------------------------------------------

# Running the MCMC

## ---------------------------------------------------------------------------------------

# Libraries we need
library(nimble)
library(coda)
library(nimbleSCR)

# Sourcing in the dpoisLocal_normal_2() function that we will need to run the MCMC
#source("dpoisLocal_normal.R")

# Sourcing in the function that we'll need to run the MCMC
source("MCMC_Function_Inhomogeneous.R")
library("spatstat")

ch7c.sample = run.MCMC.inhom(data.ch7c, M=9000, mona.column="Dgood", n.iter=100000)
save(ch7c.sample, file="ch7c.RData")
