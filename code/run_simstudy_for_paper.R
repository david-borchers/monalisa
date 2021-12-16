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
load("output/capthists.RData")

#mona_df <- mona_df %>% select(-cc) %>% distinct()
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
              secr.fitformula = "D~log(Dgood_bigD)", 
              dx = 4, dy = 4, nx = 3, ny = 4, 
              xorig = 15, yorig = 15, 
              sigma = 2, lambda0 = 0.69, noccasions = c(1), 
              capthist = parlist1$capthist[[1]])

##########
# runs for paper
##########

### runs for S4.2.1: Realised activity center densities with many activity centers

parlist1 <- expand.grid(i = 1:100,
                        secr.fitformula = "D~1", 
                        dx = 4, dy = 4, nx = 3, ny = 4, 
                        xorig = c(15,27), yorig = c(15,31), 
                        sigma = 2, lambda0 = 0.69, 
                        noccasions = c(20),
                        stringsAsFactors = FALSE)

# merge in capture histories if you have them
parlist1 <- left_join(parlist1, capthists_realised_acd_many,
                         by = c("i", "secr.fitformula", "noccasions", "xorig", "yorig"))

cl <- makeCluster(2)
registerDoParallel(cl)
getDoParWorkers()

set.seed(123)
res_realised_acd_many <- foreach(i = 1:400, .packages = c("secr", "stringr")) %dopar% {
  print(i)
  x <- run_secr(simulated_points = simulated_points_lots,
                secr.fitformula = parlist1$secr.fitformula[i], 
                dx = parlist1$dx[i], dy = parlist1$dy[i], 
                nx = parlist1$nx[i], ny = parlist1$ny[i], 
                xorig = parlist1$xorig[i], yorig = parlist1$yorig[i], 
                sigma = parlist1$sigma[i], lambda0 = parlist1$lambda0[i], 
                noccasions = parlist1$noccasions[i], 
                capthist = parlist1$capthist[[i]])
}

# set.seed(123)
# res_realised_acd_many <- foreach(i = 1:2, .packages = c("secr", "stringr", "purrr")) %dopar% {
#   print(i)
#   res <- pmap(parlist1, .f = run_secr, simulated_points = simulated_points_many)
# }

stopCluster(cl)
save(res_realised_acd_many, file = "output/res_realised_acd_many.RData")

### runs for S4.2.2: Expected activity center densities with many activity centers
parlist2a <- expand.grid(secr.fitformula = c("D~log(Dgood_bigD)", "D~log(Dblur_bigD)"),
                         dx = 4, dy = 4, nx = 3, ny = 4, 
                         xorig = 15, yorig = 15, 
                         sigma = 2, lambda0 = 0.69, 
                         noccasions = c(20),
                         stringsAsFactors = FALSE)
parlist2b <- expand.grid(secr.fitformula = c("D~log(Dgood_bigD)", "D~log(Dblur_bigD)"),
                         dx = 4, dy = 4, nx = 3, ny = 4, 
                         xorig = 27, yorig = 31, 
                         sigma = 2, lambda0 = 0.69, 
                         noccasions = c(20),
                         stringsAsFactors = FALSE)
parlist2 <- rbind(parlist2a, parlist2b)

# merge in capture histories if you have them
parlist2 <- left_join(parlist2, capthists_expected_acd_many %>% dplyr::select(-i),
                      by = c("secr.fitformula", "noccasions", "xorig", "yorig"))

set.seed(123)
res_expected_acd_many <- foreach(i = 1:4, .packages = c("secr", "stringr")) %dopar% {
  print(i)
  x <- run_secr(simulated_points = simulated_points_lots,
                secr.fitformula = parlist2$secr.fitformula[i], 
                dx = parlist2$dx[i], dy = parlist2$dy[i], 
                nx = parlist2$nx[i], ny = parlist2$ny[i], 
                xorig = parlist2$xorig[i], yorig = parlist2$yorig[i], 
                sigma = parlist2$sigma[i], lambda0 = parlist2$lambda0[i], 
                noccasions = parlist2$noccasions[i], 
                capthist = parlist2$capthist[[i]])
}

#res_expected_acd_many <- pmap(parlist2, .f = run_secr, simulated_points = simulated_points_lots, my.seed = 123)
save(res_expected_acd_many, file = "output/res_expected_acd_many.RData")

# fig5_results <- pmap(parlist2, .f = run_secr, simulated_points = simulated_points_lots, my.seed = 123)
# save(fig5_results, file = "output/res_fig5v2.RData")

### runs for S4.2.3: Realised and expected activity center densities with few activity centers

parlist3 <- expand.grid(i=1:100,
                        secr.fitformula = c("D~1", "D~log(Dgood_smallD)", "D~log(Dblur_smallD)"),
                        dx = 8, dy = 8, nx = 3, ny = 3, 
                        xorig = 10, yorig = 26, 
                        sigma = 4, lambda0 = 0.69, 
                        noccasions = c(3,10,20),
                        stringsAsFactors = FALSE)

# merge in capture histories if you have them
parlist3 <- left_join(parlist3, capthists_realised_and_expected_acd_few,
                      by = c("i", "secr.fitformula", "noccasions", "xorig", "yorig"))

sum(unlist(lapply(parlist3$capthist, is.null)))

set.seed(123)
res_realised_and_expected_acd_few <- foreach(i = 1:900, .packages = c("secr", "stringr")) %dopar% {
  print(i)
  x <- run_secr(simulated_points = simulated_points_few,
                secr.fitformula = parlist3$secr.fitformula[i], 
                dx = parlist3$dx[i], dy = parlist3$dy[i], 
                nx = parlist3$nx[i], ny = parlist3$ny[i], 
                xorig = parlist3$xorig[i], yorig = parlist3$yorig[i], 
                sigma = parlist3$sigma[i], lambda0 = parlist3$lambda0[i], 
                noccasions = parlist3$noccasions[i], 
                capthist = parlist3$capthist[[i]])
}
# set.seed(123)
# res_realised_and_expected_acd_few <- foreach(i = 1:100, .packages = c("secr", "stringr", "purrr")) %dopar% {
#   res <- pmap(parlist3, .f = run_secr, simulated_points = simulated_points_few)
# }
save(res_realised_and_expected_acd_few, file = "output/res_realised_and_expected_acd_few.RData")

# res_fig6v2 <- foreach(i = 1:100, .packages = c("secr", "stringr", "purrr")) %dopar% {
#   fig6_results <- pmap(parlist3, .f = run_secr, simulated_points = simulated_points_few)
# }
# save(res_fig6v2, file = "output/res_fig6v2.RData")

# # runs for figure 7: for the space use plots we already have most of what we need,
# # just need to get predicted ac densities for captured animals only (row 3 of plot)
# parlist4 <- expand.grid(secr.fitformula = "D~1",
#                         dx = 8, dy = 8, nx = 3, ny = 3, 
#                         xorig = 10, yorig = 26, 
#                         sigma = 4, lambda0 = 0.69, 
#                         noccasions = c(1,3,10,20),
#                         stringsAsFactors = FALSE)
# 
# fig7_results <- pmap(parlist4, .f = run_secr, simulated_points = simulated_points_few, my.seed = 123, captures.only = T)

save(res_realised_acd_many, res_expected_acd_many, res_realised_and_expected_acd_few, file = "output/mona_raw_outputs.RData")
