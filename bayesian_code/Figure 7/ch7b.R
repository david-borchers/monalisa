## To send to NeSI to generate the MCMC samples for ch7b
load("mona_raw_outputs_100sim.RData") # for NeSI
source("posthoc_extract_chs.R") # for NeSI

# Encounter matrix
encounterdat.ch7b = matrix(0, nrow=nrow(ch7b[,1,]), ncol=ncol(ch7b[,1,]))
for (i in 1:20) {
  encounterdat.ch7b = encounterdat.ch7b + ch7b[,i,]
}
# Trap locations
ch7b.traploc = attributes(ch7b)$traps
# xlim, ylim
xlim = c(0.5, 50.5)
ylim = c(0.5, 50.5)
# Data object
data.ch7b = list(encounter.data = encounterdat.ch7b, trap.loc = ch7b.traploc, xlim = xlim, ylim = ylim, n.occasions = 20)

# Libraries we need
library(nimble)
library(coda)
library(nimbleSCR)

# Sourcing in the function that we'll need to run the MCMC
source("MCMC_Function_Inhomogeneous.R")
library("spatstat")

ch7b.sample = run.MCMC.inhom(data=data.ch7b, M=9000, mona.column="Dgood", n.iter=20000, n.burn=0)
save(ch7b.sample, file="ch7b.RData")
