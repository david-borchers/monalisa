library(dplyr)
library(secr)

load("output/simulated_densities.Rdata")

simulated_points <- simulated_points_few

n_pts <- nrow(simulated_points) 

# make a grid of detectors
detectors <- make.grid(nx = 3, ny = 3, spacex = 8, spacey = 8,
                       originxy = c(10, 16), detector = "count")

# capture history for maximum number of occasions 
sigma = 4
lambda0 = 0.69
max_nocc <- 80
my.seed = 112233
capture_history_max_occasions <- list()
for(i in 1:100){
  capture_history_max_occasions[[i]] <- sim.capthist(detectors, popn = simulated_points, detectfn = "HHN",
                                                     detectpar = list(lambda0 = lambda0, sigma = sigma),
                                                     noccasions = max_nocc,
                                                     nsessions = 1,
                                                     seed = my.seed)
}


load("output/capthists-oldfewacs.RData")

ch3 <- lapply(capture_history_max_occasions, subset, occasions = 1:3, dropunused = FALSE)
ch10 <- lapply(capture_history_max_occasions, subset, occasions = 1:10, dropunused = FALSE)
ch20 <- lapply(capture_history_max_occasions, subset, occasions = 1:20, dropunused = FALSE)
ch40 <- lapply(capture_history_max_occasions, subset, occasions = 1:40, dropunused = FALSE)

cht <- tibble(expand.grid(i=1:100, noccasions = c(3,10,20,40,80)))
cht <- cht %>% mutate(capthist = c(ch3, ch10, ch20, ch40, capture_history_max_occasions))

capthists_realised_and_expected_acd_few <- tibble(expand.grid(i=1:100,
                                                              xorig = 10, yorig = 16,
                                                              secr.fitformula = c("D~1", "D~log(Dgood_smallD)", "D~log(Dblur_smallD)"),
                                                              noccasions = c(3,10,20,40,80)))
capthists_realised_and_expected_acd_few <- left_join(capthists_realised_and_expected_acd_few, cht, by = c("i", "noccasions"))

save(capthists_expected_acd_many, capthists_realised_acd_many, capthists_realised_and_expected_acd_few,
     file = "output/capthists.RData")
