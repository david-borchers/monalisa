## Code to create uncertainty plots for the EACD plots shown in Figures 4 and 5

## Libraries we need
library(nimble)
library(coda)
library(nimbleSCR)
library(spatstat)
library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(secr)
library(patchwork)
library(ggpubr)
library(viridis)
## Objects we need
load("../output/revision/mona-inputs.RData")
load("../output/revision/mona-results.RData")
## Functions we need
source("Functions.R")

## For each EACD plot in Figures 4 and 5, we will create two uncertainty plots: one showing the lower 5% quantile for the posterior distribution of the density for each pixel, the other will show the upper 95% quantile.
## We will create two figures. One will show the uncertainty plots from Figure 4, along with the corresponding EACD plots. The other figure will show the same thing, but for Figure 5. We will refer to these figures as 'uncertainty figures' in the code below. 

## First, we will create the objects needed to create both figures. We will then put together both figures at the end. This is because we want the same maximum value to be used when colouring both figures.

## ---------------------------------------------------------------------------------------
#################### Objects needed for uncertainty figure for Figure 4 ##################
## ---------------------------------------------------------------------------------------

## Loading in the MCMC results from Row 2 of Figure 4
load("MCMC_Results/Figure4/InhomPP_18occ.RData")
load("MCMC_Results/Figure4/InhomPP_52occ.RData")
load("MCMC_Results/Figure4/InhomPP_111occ.RData")
# Discarding burn-in
inhom.results.18occ <- inhom.results.18occ[-c(1:1000),]
inhom.results.52occ <- inhom.results.52occ[-c(1:1000),]
inhom.results.111occ <- inhom.results.111occ[-c(1:1000),]

## Covariate that we use for these EACD plots (see 'Functions.R' for an explanation of the function)
log.dblur <- eacd.covariate()

## Creating the data frames that will contain all of the info required to create the uncertainty figure for Figure 4
# 18 sampling occasions
info.18occ <- eacd.quantile.info(results=inhom.results.18occ, covariate=log.dblur, nPix=2500, nocc=18)
# 52 sampling occasions
info.52occ <- eacd.quantile.info(results=inhom.results.52occ, covariate=log.dblur, nPix=2500, nocc=52)
# 111 sampling occasions
info.111occ <- eacd.quantile.info(results=inhom.results.111occ, covariate=log.dblur, nPix=2500, nocc=111)

## ---------------------------------------------------------------------------------------
#################### Objects needed for uncertainty figure for Figure 5 ##################
## ---------------------------------------------------------------------------------------

## Loading in the MCMC results from Row 2 of Figure 5
load("MCMC_Results/Figure5/InhomPP_7occ.RData")
load("MCMC_Results/Figure5/InhomPP_25occ.RData")
load("MCMC_Results/Figure5/InhomPP_55occ.RData")
# Discarding burn-in
inhom.results.7occ <- inhom.results.7occ[-c(1:1000),]
inhom.results.25occ <- inhom.results.25occ[-c(1:1000),]
inhom.results.55occ <- inhom.results.55occ[-c(1:1000),]

## Covariate we use is the same as above
log.dblur <- eacd.covariate()

## Creating the data frames that will contain all of the info require to create the uncertainty figure for Figure 5
# 7 sampling occasions
info.7occ <- eacd.quantile.info(results=inhom.results.7occ, covariate=log.dblur, nPix=2500, nocc=7)
# 25 sampling occasions
info.25occ <- eacd.quantile.info(results=inhom.results.25occ, covariate=log.dblur, nPix=2500, nocc=25)
# 55 sampling occasions
info.55occ <- eacd.quantile.info(results=inhom.results.55occ, covariate=log.dblur, nPix=2500, nocc=55)

## ---------------------------------------------------------------------------------------
####################### Objects needed for both uncertainty figures ######################
## ---------------------------------------------------------------------------------------

## Collating all information for figures into one data frame
values_all <- rbind(info.18occ, info.52occ, info.111occ, info.7occ, info.25occ, info.55occ)

## Detectors
detectors_df_all <- res_acd %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

## Maximum value for colour scale for *all* plots
maxval <- max(values_all$value)

## ---------------------------------------------------------------------------------------
######################## Creating uncertainty figure for Figure 4 ########################
## ---------------------------------------------------------------------------------------

nn <- 3 # Number of different simulated datasets used in Figure 4
occ <- capthists_few_alloccs_3x3$noccasions # Number of sampling occasions for each dataset
asz <- c("3x3")
chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_3x3$capthist, summary, terse = TRUE)))
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist() # Row labels we'll use

## Adding faceting info to 'values_all' for rows that we will use in this uncertainty figure
values_all$covtype2 <- factor(values_all$covtype, levels=unique(values_all$covtype),
                              labels=c("Expected AC"))
values_all$quantile <- factor(rep(c(rep("0.05 quantile", 2500), rep("Mean", 2500), rep("0.95 quantile", 2500)), 3),
                              levels=c("0.05 quantile", "Mean", "0.95 quantile"))
values_all$occasions2 <- factor(values_all$occasions,
                                levels = occ,
                                labels = capthist_labels)

## Creating and saving the uncertainty figure
uncertainty.fig4 <- values_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_distiller(limits=c(0, maxval)) + 
  facet_grid(occasions2 ~ quantile) +
  geom_point(data = detectors_df_all %>% filter(occasions %in% occ[1:nn], array_size %in% asz), inherit.aes = T,
             colour = "gray80", pch = 4, size = 2) +
  geom_point(data = simulated_points, inherit.aes = F, aes(x=x,y=y),
             colour = "darkorange", pch = 16, size = 1, alpha = 0.5) +
  coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.spacing=unit(-1, "lines"),
        strip.background = element_rect(fill=NA, colour = NA), 
        legend.position="right", legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1.3,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("mona_3x3_uncertainty.png", uncertainty.fig4, width=8, height=8, dpi=600, bg="white")

## ---------------------------------------------------------------------------------------
######################## Creating uncertainty figure for Figure 5 ########################
## ---------------------------------------------------------------------------------------

nn <- 3
occ <-capthists_few_alloccs_7x7$noccasions
asz <- c("7x7")
chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_7x7$capthist, summary, terse = TRUE)))
chs <- chs %>% dplyr::filter(Occasions %in% occ)
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist()

## Adding faceting info for rows that we will use to create this uncertainty figure  
values_all$occasions2 <- factor(values_all$occasions,
                                levels = occ,
                                labels = capthist_labels)

## Creating and saving the uncertainty figure
uncertainty.fig5 <- values_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_distiller(limits=c(0, maxval)) + 
  facet_grid(occasions2 ~ quantile) +
  geom_point(data = detectors_df_all %>% filter(occasions %in% occ[1:nn], array_size %in% asz), inherit.aes = T,
             colour = "gray80", pch = 4, size = 2) +
  geom_point(data = simulated_points, inherit.aes = F, aes(x=x,y=y),
             colour = "darkorange", pch = 16, size = 1, alpha = 0.5) +
  coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.spacing=unit(-1, "lines"),
        strip.background = element_rect(fill=NA, colour = NA), 
        legend.position="right", legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1.3,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("mona_7x7_uncertainty.png", uncertainty.fig5, width=8, height=8, dpi=600, bg="white")
