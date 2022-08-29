library(dplyr)
library(secr)

load("output/simulated_densities.Rdata")

simulated_points <- simulated_points_few

n_pts <- nrow(simulated_points) 

# make a 3x3 grid of detectors
detectors3x3 <- make.grid(nx = 3, ny = 3, spacex = 8, spacey = 8,
                          originxy = c(10, 16), detector = "count")

# make a 7x7 grid of detectors
detectors7x7 <- make.grid(nx = 7, ny = 7, spacex = 8, spacey = 8,
                            originxy = c(1.5, 1.5), detector = "count")

## make capthist for 3x3 grid

# capture history for maximum number of occasions 
sigma = 4
lambda0 = 0.1
max_nocc <- 701
my.seed = 112233
capture_history_max_occasions <- sim.capthist(detectors3x3, popn = simulated_points, detectfn = "HHN",
                                              detectpar = list(lambda0 = lambda0, sigma = sigma),
                                              noccasions = max_nocc,
                                              nsessions = 1,
                                              seed = my.seed)

# manually check capthists
i <- 5
chc <- subset(capture_history_max_occasions, occasions = 1:i, dropunused = FALSE)
summary(chc, terse = TRUE)

occs3x3 <- c(18, 52, 111, seq(5, 544, length = 50))
cht <- list()
j <- 0
for(i in occs3x3){
  j <- j + 1
  cht[[j]] <- subset(capture_history_max_occasions, occasions = 1:i, dropunused = FALSE)
}

cht2 <- tibble(expand.grid(i=1:1, noccasions = occs3x3))
cht2 <- cht2 %>% mutate(capthist = cht)
rm(cht)

capthists_few_alloccs <- tibble(expand.grid(i=1,
                                            xorig = 10, yorig = 16,
                                            secr.fitformula = c("D~1"),
                                            noccasions = occs3x3))
capthists_few_alloccs <- left_join(capthists_few_alloccs, cht2, by = c("i", "noccasions"))
rm(cht2)

capthists_few_alloccs_3x3 <- capthists_few_alloccs

## make capthist for 7x7 grid

# capture history for maximum number of occasions 
sigma = 4
lambda0 = 0.1
max_nocc <- 301
my.seed = 112233
capture_history_max_occasions <- sim.capthist(detectors7x7, popn = simulated_points, detectfn = "HHN",
                                              detectpar = list(lambda0 = lambda0, sigma = sigma),
                                              noccasions = max_nocc,
                                              nsessions = 1,
                                              seed = my.seed)

# manually check capthists
i <- 54
chc <- subset(capture_history_max_occasions, occasions = 1:i, dropunused = FALSE)
summary(chc, terse = TRUE)

occs7x7 <- c(7, 25, 55, seq(2, 202, length.out = 51))
cht <- list()
j <- 0
for(i in occs7x7){
  j <- j + 1
  cht[[j]] <- subset(capture_history_max_occasions, occasions = 1:i, dropunused = FALSE)
}

cht2 <- tibble(expand.grid(i=1:1, noccasions = occs7x7))
cht2 <- cht2 %>% mutate(capthist = cht)
rm(cht)

capthists_few_alloccs <- tibble(expand.grid(i=1,
                                            xorig = 1.5, yorig = 1.5,
                                            secr.fitformula = c("D~1"),
                                            noccasions = occs7x7))
capthists_few_alloccs <- left_join(capthists_few_alloccs, cht2, by = c("i", "noccasions"))
rm(cht2)

capthists_few_alloccs_7x7 <- capthists_few_alloccs

rm(capthists_few_alloccs)

save(capthists_few_alloccs_7x7, capthists_few_alloccs_3x3, file = "output/capthists_manyocc.RData")
