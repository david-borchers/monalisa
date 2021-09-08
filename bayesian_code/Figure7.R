## Code to create Figure 7 of the paper

## ---------------------------------------------------------------------------------------

# Creating the necessary data objects

## ---------------------------------------------------------------------------------------

## First, sourcing in 'posthoc_extract_chs.R'
source("../code/posthoc_extract_chs.R")
# source("posthoc_extract_chs.R") for NeSI

## So, it appears that:
# ch7b is the plot from traps placed at the top right
# ch7c is the plot from traps placed at the bottom left
# ch7e is the plot from traps placed at the top right
# ch7f is the plot from traps placed at the bottom left

## Are going to create 4 EACD maps, two sets.
# One set uses the 'Dgood' covariate, and uses the trap array at the top right and the bottom left.
# The other uses the 'Dblur' covariate, and also uses the trap array at the top right and the bottom left.

## So, we need to use MCMC to fit an inhomogeneous density model -- i.e. to generate MCMC samples for beta0 and beta1, which we will then use to estimate the inhomogeneous PP, and therefore to create the maps.

## We begin by creating the necessary data objects. Here, we are working with one set of simulated data only, which has 20 sampling occasions. So, for each map, we need to create data objects that contain an encounter data matrix found by summing over these 20 sampling occasions

## ch7b
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

## ch7c
encounterdat.ch7c = matrix(0, nrow=nrow(ch7c[,1,]), ncol=ncol(ch7c[,1,]))
for (i in 1:20) {
  encounterdat.ch7c = encounterdat.ch7c + ch7c[,i,]
}
ch7c.traploc = attributes(ch7c)$traps
data.ch7c = list(encounter.data = encounterdat.ch7c, trap.loc = ch7c.traploc, xlim = xlim, ylim = ylim, n.occasions = 20)

## ch7e
encounterdat.ch7e = matrix(0, nrow=nrow(ch7e[,1,]), ncol=ncol(ch7e[,1,]))
for (i in 1:20) {
  encounterdat.ch7e = encounterdat.ch7e + ch7e[,i,]
}
ch7e.traploc = attributes(ch7e)$traps
data.ch7e = list(encounter.data = encounterdat.ch7e, trap.loc = ch7e.traploc, xlim = xlim, ylim = ylim, n.occasions = 20)

## ch7f
encounterdat.ch7f = matrix(0, nrow=nrow(ch7f[,1,]), ncol=ncol(ch7f[,1,]))
for (i in 1:20) {
  encounterdat.ch7f = encounterdat.ch7f + ch7f[,i,]
}
ch7f.traploc = attributes(ch7f)$traps
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

## ------------------------------

# ch7b

## ------------------------------


ch7b.sample = run.MCMC.inhom(data=data.ch7b, M=9000, mona.column="Dblur")
save(ch7b.sample, file="ch7b.RData")
ch7e.sample = run.MCMC.inhom(data.ch7e, M=9000, mona.column="Dblur")
save(ch7e.sample, file="ch7e.RData")
ch7c.sample = run.MCMC.inhom(data.ch7c, M=9000, mona.column="Dgood")
save(ch7c.sample, file="ch7c.RData")
ch7f.sample = run.MCMC.inhom(data.ch7f, M=9000, mona.column="Dgood")
save(ch7f.sample, file="ch7f.RData")

