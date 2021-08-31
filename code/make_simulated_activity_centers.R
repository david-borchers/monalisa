# generates the two sets of simulated activity centers used in the paper

library(tidyverse)
# library(devtools)
#install_github("rachelphillip/SCR-Book/scrmlebook")
library(scrmlebook)

source("code/make_simulated_points_plottable.R")

load("output/mona_inputs.RData")

mlmesh <- read.mask(data = mona_df)

# simulate points (activity centers) from the true density 

# to get a decent representation of the image, you need a lot of points
# each cell is 1m2 and D is density / 10000m2
# so expected number of points in each cell is D / 10000 and 
# expected total points generated is sum(D) / 10000. 
# If we want n_pts points then we need to multiply D by 10000 * (n_pts / sum(D))

# large number of activity centers (dataset 1 in paper)
n_pts <- 7500  # desired number of points to generate, in expectation
D_for_sim = covariates(mlmesh)$D / sum(covariates(mlmesh)$D) * 10000 * n_pts
simulated_points_lots <- sim.popn(D = D_for_sim, 
                                  core = mlmesh, 
                                  model2D = "IHP", 
                                  nDist = "poisson",
                                  seed = 123)

# small number of activity centers (dataset 2 in paper)
n_pts <- 80  # desired number of points to generate, in expectation
D_for_sim = covariates(mlmesh)$D / sum(covariates(mlmesh)$D) * 10000 * n_pts
simulated_points_few <- sim.popn(D = D_for_sim, 
                                 core = mlmesh, 
                                 model2D = "IHP", 
                                 nDist = "poisson",
                                 seed = 12345)

# if you want to use these later e.g. for plotting, then save here
simulated_densities_df <- make_simulated_points_plottable(simulated_points_lots)
simulated_densities_small_df <- make_simulated_points_plottable(simulated_points_few)

# add zeros where no observations in a cell
simulated_densities_small_df <- simulated_densities_small_df %>% 
  complete(x  = 1:50, y = 1:50, fill = list(value = 0))

# name changed from "output/simulated_densities_for_paper.Rdata" in tidy up (31/8/2021)
save(simulated_points_lots, simulated_points_few,
     simulated_densities_df, simulated_densities_small_df,
     file = "output/simulated_densities.Rdata")

