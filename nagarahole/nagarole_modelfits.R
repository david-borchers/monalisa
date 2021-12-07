library(secr)
library(dplyr)
library(ggplot2)
library(stringr)
library(rgdal)
library(viridis)
library(maptools)
library(colorspace)
library(pals)
library(gridExtra)
library(purrr)

source("code/predicted_densities_for_D0.R")

tigerch = read.capthist("nagarahole/data/NagaraholeCHtimes.csv","nagarahole/data/Nagaraholetraps.csv",detector="count",noccasions=1,covnames=c("day","hour"))
cams = traps(tigerch)
if (dim(tigerch)[2]>1) {
  par(mfrow=c(3,3))
  for(i in 1:dim(tigerch)[2]) plot(tigerch[[i]],border=0,tracks=TRUE)
} else {
  plot(tigerch,border=0,tracks=TRUE)
}

# Get boundary polygon and create mask
sregion = rgdal::readOGR("nagarahole/data/NHstatespace-utm.shp")
tigermask = make.mask(traps(tigerch),buffer=10000,spacing=500,type="trapbuffer",poly=sregion) # buffer and spacing consistent with Dorazio+Kranth (2017)

plot(tigermask)
plot(cams,add =T)

#### ALL TRAPS

# fit null model (all traps)
fit0 = secr.fit(tigerch,
                mask=tigermask)

# fit D ~ y model (all traps)
fity = secr.fit(tigerch,
                mask=tigermask,
                model=list(D~y),
                start=list(g0=detectpar(fit0)$g0,sigma=detectpar(fit0)$sigma))

# predicted density surface for D ~ 1 (all traps), note these are standardised to sum to one
# (see the predicted_densities_for_D0 function for details)
predicted_densities_0  <- predicted_densities_for_D0(fit0, cams, tigermask)
predicted_densities_0$traps <- "All traps, no cov."

# predicted intensity surface for D ~ y (all traps)
predy = predictDsurface(fity,mask=tigermask,se.D=TRUE,cl.D=TRUE)
covariates(predy)$cv = (covariates(predy)$SE.0/covariates(predy)$D.0)*100
covariates(predy)[c("D.0","SE.0","lcl.0","ucl.0")] = covariates(predy)[c("D.0","SE.0","lcl.0","ucl.0")]*100^2
predicted_densities_y <- data.frame(x = predy$x, y = predy$y, value = covariates(predy)$D.0, 
                                    traps = "All traps, northing")

#### REDUCED TRAPS

# reduced trap array #1: only keeps traps in the extreme S, N, and E

# get reduced traps 
cams$id <- 1:nrow(cams)
# find some quantiles of northing and easting, to identify where we want to keep traps
quantile(cams$y, probs = c(0.1,0.9))
quantile(cams$x, probs = 0.9)
cams_to_include <- cams %>% 
  filter((y > 1345475) | # keep northernmost cams
           (x < 625000 & y < 1316790) | # keep south-westernmost cams
           (y > 1320000 & y < 1340000 & x > 631126.7)) %>% # keep easternmost cams
  select(id) %>% unlist()
reduced_cams <- cams[cams_to_include,]

# fit null model (reduced traps)
fit0_rt = secr.fit(subset(tigerch, traps = reduced_cams$id),
                   mask=tigermask)

# fit D ~ y model (reduced traps)
fity_rt = secr.fit(subset(tigerch, traps = reduced_cams$id),
                   mask=tigermask,
                   model=list(D~y),
                   start=list(g0=detectpar(fit0_rt)$g0,sigma=detectpar(fit0_rt)$sigma))

# predicted density surface for D ~ 1 (reduced traps)
predicted_densities_0_rt  <- predicted_densities_for_D0(fit0_rt, reduced_cams, tigermask)
# add variables for later plotting
predicted_densities_0_rt$traps <- "Subset #1, no cov."

# predicted intensity surface for D ~ y (reduced traps)
predy = predictDsurface(fity_rt, mask=tigermask, se.D=TRUE, cl.D=TRUE)
covariates(predy)$cv = (covariates(predy)$SE.0/covariates(predy)$D.0)*100
covariates(predy)[c("D.0","SE.0","lcl.0","ucl.0")] = covariates(predy)[c("D.0","SE.0","lcl.0","ucl.0")]*100^2
predicted_densities_y_rt <- data.frame(x = predy$x, y = predy$y, 
                                       value = covariates(predy)$D.0, 
                                       traps = "Subset #1, northing")

# reduced trap array #2: removes a bunch of traps closest to two of the highest-density 
# points in the "all traps" density surface. Simulates the case where, by chance, 
# we did not have traps close to where there is an activity center.

# get reduced traps 
cams$id <- 1:nrow(cams)
# y-coordinates of the cams "at" the two activity centers (very close)
y_excl <- c(1323864, 1330531)
# picks out the two cams at the center of the area we want to clear of cams
highD_cams <- cams %>% filter(round(y, 0) %in% y_excl)
# remove the 8 cams closest to the activity center (can adjust the "8" and see what happens)
cams_to_include <- cams %>% 
  mutate(d1 = rank((x - highD_cams$x[1])^2 + (y - highD_cams$y[1])^2),
         d2 = rank((x - highD_cams$x[2])^2 + (y - highD_cams$y[2])^2)) %>% 
  filter(d1 > 8, d2 > 8) %>% select(id) %>% unlist()
reduced_cams2 <- cams[cams_to_include,]

# fit null model (reduced traps)
fit0_rt2 = secr.fit(subset(tigerch, traps = reduced_cams2$id),
                    mask=tigermask)

# fit D ~ y model (reduced traps)
fity_rt2 = secr.fit(subset(tigerch, traps = reduced_cams2$id),
                    mask=tigermask,
                    model=list(D~y),
                    start=list(g0=detectpar(fit0_rt2)$g0,sigma=detectpar(fit0_rt2)$sigma))

# predicted density surface for D ~ 1 (reduced traps)
predicted_densities_0_rt2  <- predicted_densities_for_D0(fit0_rt2, reduced_cams2, tigermask)
predicted_densities_0_rt2$traps <- "Subset #2, no cov."

# predicted intensity surface for D ~ y (reduced traps)
predy = predictDsurface(fity_rt2, mask=tigermask, se.D=TRUE, cl.D=TRUE)
covariates(predy)$cv = (covariates(predy)$SE.0/covariates(predy)$D.0)*100
covariates(predy)[c("D.0","SE.0","lcl.0","ucl.0")] = covariates(predy)[c("D.0","SE.0","lcl.0","ucl.0")]*100^2
predicted_densities_y_rt2 <- data.frame(x = predy$x, y = predy$y, value = covariates(predy)$D.0, 
                                        traps = "Subset #2, northing")

save.image("nagarahole/output/nagarole_modelruns.RData")