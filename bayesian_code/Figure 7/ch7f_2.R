## Creating ch7f RData object using function

## First, sourcing in 'posthoc_extract_chs.R'
#load("../output/mona_raw_outputs_100sim.RData")
#source("../code/posthoc_extract_chs.R")
load("mona_raw_outputs_100sim.RData") # for NeSI
source("posthoc_extract_chs.R") # for NeSI

## ch7f
encounterdat.ch7f = matrix(0, nrow=nrow(ch7f[,1,]), ncol=ncol(ch7f[,1,]))
for (i in 1:20) {
  encounterdat.ch7f = encounterdat.ch7f + ch7f[,i,]
}
ch7f.traploc = attributes(ch7f)$traps
xlim = c(0.5, 50.5)
ylim = c(0.5, 50.5)
data.ch7f = list(encounter.data = encounterdat.ch7f, trap.loc = ch7f.traploc, xlim = xlim, ylim = ylim, n.occasions = 20)

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

ch7f.sample = run.MCMC.inhom(data.ch7f, M=9000, mona.column="Dblur", n.iter=100000)
save(ch7f.sample, file="ch7f.RData")
