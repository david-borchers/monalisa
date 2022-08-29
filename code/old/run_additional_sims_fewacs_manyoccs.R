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
load("output/capthists_manyocc.RData")
load("output/simulated_densities.Rdata")

simulated_points <- simulated_points_few
mlmesh <- read.mask(data = mona_df)

summary(capthists_few_alloccs_3x3$capthist[[22]], terse = TRUE)
summary(capthists_few_alloccs_7x7$capthist[[24]], terse = TRUE)

### runs for S4.2.3: Realised and expected activity center densities with few activity centers

occs3x3 <- capthists_few_alloccs_3x3$noccasions
parlist3a <- expand.grid(i=1:1,
                         secr.fitformula = c("D~1", "D~log(Dblur_bigD)"),
                         dx = 8, dy = 8, nx = 3, ny = 3, 
                         xorig = 10, yorig = 16, 
                         sigma = 4, lambda0 = 0.1, 
                         noccasions = occs3x3,
                         stringsAsFactors = FALSE)

# merge in capture histories if you have them
parlist3a <- left_join(parlist3a, capthists_few_alloccs_3x3,
                       by = c("i", "noccasions", "xorig", "yorig"))

occs7x7 <- capthists_few_alloccs_7x7$noccasions
parlist3b <- expand.grid(i=1:1,
                         secr.fitformula = c("D~1","D~log(Dblur_bigD)"),
                         dx = 8, dy = 8, nx = 7, ny = 7, 
                         xorig = 1.5, yorig = 1.5, 
                         sigma = 4, lambda0 = 0.1, 
                         noccasions = occs7x7,
                         stringsAsFactors = FALSE)

# merge in capture histories if you have them
parlist3b <- left_join(parlist3b, capthists_few_alloccs_7x7,
                       by = c("i", "noccasions", "xorig", "yorig"))

parlist3 <- rbind(parlist3a, parlist3b)
rm(parlist3a, parlist3b, capthists_few_alloccs_3x3, capthists_few_alloccs_7x7)

#parlist3 <- parlist3[-c(1,2), ]
sum(unlist(lapply(parlist3$capthist, is.null)))

# just running the 3 we need for 3x3 and 7x7 main plot
parlist3 <- parlist3 %>% filter( (nx == 3 & noccasions %in% c(occs3x3[c(1:3,33)]))|
                                   (nx == 7 & noccasions %in% c(occs7x7[c(1:3,32)])))

parlist3 <- parlist3 %>% filter( (nx == 3 & noccasions == occs3x3[22])|
                                   (nx == 7 & noccasions == occs7x7[24]))

set.seed(123)
#res_acd_few_alloccs <- foreach(i = 1:279, .packages = c("secr", "stringr")) %dopar% {
res_acd_few <- list()
for(i in 1:nrow(parlist3)){
  print(i)
  res_acd_few[[16+i]] <- run_secr(simulated_points = simulated_points_few,
                                  mask = mlmesh,
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
save(res_acd_few, file = "output/res_acd_few.RData")
