# master file

library(dplyr)
library(stringr)
library(purrr)
#library(devtools)
#install_github("rachelphillip/SCR-Book/scrmlebook")
#library(scrmlebook)
library(secr)
library(doParallel)

source("code/predicted_densities_for_D0.R")
source("code/run_secr.R")

load("output/mona_inputs.RData")

mona_df <- mona_df %>% select(-cc) %>% distinct()
mlmesh <- read.mask(data = mona_df)

# simulate activity centers from the true density 

#####
# to get a decent representation of the image, you need a lot of points
# each cell is 1m2 and D is density / 10000m2
# so expected number of points in each cell is D / 10000 and 
# expected total points generated is sum(D) / 10000. 
# If we want n_pts points then we need to multiply D by 10000 * (n_pts / sum(D))
#####

expn_pts <- 750  # desired number of points to generate, in expectation
D_for_sim = covariates(mlmesh)$D / sum(covariates(mlmesh)$D) * 10000 * expn_pts
my_simulated_points <- sim.popn(D = D_for_sim, 
                                core = mlmesh, 
                                model2D = "IHP", 
                                nDist = "poisson",
                                seed = 123)

# I load the activity centers generated for the paper, comment out as desired
load("output/simulated_densities.Rdata")

# example of single run (not used in paper) 
x <- run_secr(simulated_points = simulated_points_lots,
              secr.fitformula = "D~1", 
              dx = 4, dy = 4, nx = 3, ny = 4, 
              xorig = 15, yorig = 15, 
              sigma = 2, lambda0 = 0.69, noccasions = c(1))

# runs for figure 3 and 4: no covariate, after 1 and 20 occasions
parlist1 <- expand.grid(secr.fitformula = "D~1", 
                        dx = 4, dy = 4, nx = 3, ny = 4, 
                        xorig = c(15,27), yorig = c(15,31), 
                        sigma = 2, lambda0 = 0.69, 
                        noccasions = c(1,20),
                        stringsAsFactors = FALSE)

cl <- makeCluster(2)
registerDoParallel(cl)
getDoParWorkers()

set.seed(123)
res <- foreach(i = 1:100, .packages = c("secr", "stringr", "purrr")) %dopar% {
  print(i)
  fig34_results <- pmap(parlist1[1,], .f = run_secr, simulated_points = simulated_points_few)
}

stopCluster(cl)

# runs for figure 5: various covariate, after 1 and 20 occasions, bot-left and top-rt arrays only 
parlist2a <- expand.grid(secr.fitformula = c("D~Dgood", "D~Dblur", "D~Dldv"),
                         dx = 4, dy = 4, nx = 3, ny = 4, 
                         xorig = 15, yorig = 15, 
                         sigma = 2, lambda0 = 0.69, 
                         noccasions = c(1,20),
                         stringsAsFactors = FALSE)
parlist2b <- expand.grid(secr.fitformula = c("D~Dgood", "D~Dblur", "D~Dldv"),
                         dx = 4, dy = 4, nx = 3, ny = 4, 
                         xorig = 27, yorig = 31, 
                         sigma = 2, lambda0 = 0.69, 
                         noccasions = c(1,20),
                         stringsAsFactors = FALSE)
parlist2 <- as.list(rbind(parlist2a, parlist2b))

fig5_results <- pmap(parlist2, .f = run_secr, simulated_points = simulated_points_lots, my.seed = 123)
save(fig5_results, file = "output/res_fig5v2.RData")

# runs for figure 6: fewer activity centers and different array (more spaced, top-left),
# else as for figure 5
#parlist3 <- expand.grid(secr.fitformula = c("D~1", "D~Dgood", "D~Dblur")
parlist3 <- expand.grid(secr.fitformula = c("D~1", "D~Dgood", "D~Dblur"),
                        dx = 8, dy = 8, nx = 3, ny = 3, 
                        xorig = 10, yorig = 26, 
                        sigma = 4, lambda0 = 0.69, 
                        noccasions = c(1,3,10,20),
                        stringsAsFactors = FALSE)

set.seed(123)
res_fig6v2 <- foreach(i = 1:100, .packages = c("secr", "stringr", "purrr")) %dopar% {
  fig6_results <- pmap(parlist3, .f = run_secr, simulated_points = simulated_points_few)
}
save(res_fig6v2, file = "output/res_fig6v2.RData")

# runs for figure 7: for the space use plots we already have most of what we need,
# just need to get predicted ac densities for captured animals only (row 3 of plot)
parlist4 <- expand.grid(secr.fitformula = "D~1",
                        dx = 8, dy = 8, nx = 3, ny = 3, 
                        xorig = 10, yorig = 26, 
                        sigma = 4, lambda0 = 0.69, 
                        noccasions = c(1,3,10,20),
                        stringsAsFactors = FALSE)

fig7_results <- pmap(parlist4, .f = run_secr, simulated_points = simulated_points_few, my.seed = 123, captures.only = T)

save(fig34_results, fig5_results, fig6_results, fig7_results, file = "output/mona_raw_outputs.RData")
