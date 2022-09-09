## Code to create uncertainty plots for the plots shown in Figures 4 and 5

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

## ---------------------------------------------------------------------------------------
####################### Uncertainty plots for Figure 4 (attempt 1) #######################
## ---------------------------------------------------------------------------------------

## Loading in the MCMC results for Row 1 of Figure 4 (note that as stated in 'bayesian_code/Plots_Code.R', we don't need to discard any burn-in, this has been done for us)
load("MCMC_Results/Figure4/HomPP_18occ.RData")
load("MCMC_Results/Figure4/HomPP_52occ.RData")
load("MCMC_Results/Figure4/HomPP_111occ.RData")

## When dealing with the RACD maps, the maps in the middle column will show posterior mean of number of activity centres in each pixel. This is what RACD maps show. Therefore, the middle column will just consist of RACD maps!
## The first column will contain maps where the pixels are coloured using the value of the 5% quantile for the credible interval of the number of activity centres in each pixel.
## The third column will contain maps where the pixels are coloured using the value of the 95% quantile for the credible interval of the number of activity centres in each pixel. 

## Function to create plots in the first and third columns
# Arguments are the same as 'no.movement.density.vector()', with the addition of 'quantile'. Use this argument to provide the quantile of the credible interval for the number of activity centres in each pixel that we want to use (when colouring the pixels). 
racd.uncertainty.plots <- function(xlim, ylim, results, M, quantile) {

  ## Points at which local density will be estimated
  xg <- seq(xlim[1], xlim[2], by=1)
  yg <- seq(ylim[1], ylim[2], by=1)

  ## Extracting z-values
  # Names of variables that have been monitored
  names <- names(results[1,])
  # Extracting "z" values from MCMC results
  Z <- results[,grep("z", names)]

  ## Extracting activity centres
  # Extracting "s" values from MCMC results (i.e. extracting sampled activity centres)
  S <- results[,grep("s[^i]", names)]
  # x-coordinates of all activity centres
  Sx <- S[,1:M]
  # y-coordinates of all activity centres
  Sy <- S[,-(1:M)]

  ## For each MCMC iteration, storing the number of animals alive and with their activity centres in each cell -- building up posterior distribution of number of activity centres in each pixel
  # Number of pixel centres
  npix <-  (length(xg) - 1) * (length(yg) - 1)
  Dn.vals <-  matrix(0, nrow=nrow(results), ncol=npix)
  for (i in 1:nrow(results)) {
    if ((i %% 100) == 0) print(i) # Track progress
    Sxout <-  Sx[i,][Z[i,] == 1]
    Sxout <-  cut(Sxout, breaks=xg, include.lowest=TRUE)
    Syout <-  Sy[i,][Z[i,] == 1]
    Syout <-  cut(Syout, breaks=yg, include.lowest=TRUE)
    Dn.vals[i,] <-  as.vector(table(Sxout, Syout))
  }
  #browser()
  ## Finding specified quantile of credible interval for number of activity centres in each cell
  quantile.vector <- apply(Dn.vals, 2, function(i) quantile(i, probs=quantile))
}

test1 <- racd.uncertainty.plots(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.18occ, M=300, quantile=0.05)
test2 <- racd.uncertainty.plots(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.18occ, M=300, quantile=0.95)
racd <- no.movement.density.vector(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.18occ, M=300)

pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)

test1.dat <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~1", 2500), occasions=rep(18, 2500), array_size=rep("3x3", 2500), value=test1)
test2.dat <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~1", 2500), occasions=rep(18, 2500), array_size=rep("3x3", 2500), value=test2)
racd.dat <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~1", 2500), occasions=rep(18, 2500), array_size=rep("3x3", 2500), value=racd)
pred_all <- rbind(test1.dat, test2.dat, racd.dat)

detectors_df_all <- res_acd %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

pred_all$covtype2 <- factor(pred_all$covtype, levels=unique(pred_all$covtype),
                            labels=c("Realised AC"))
pred_all$quantile <- factor(c(rep("0.05", 2500), rep("0.95", 2500), rep("mean", 2500)),
                            levels=c("0.05", "mean", "0.95"))


#info <- racd.summary(nocc=18, fig=4)

maxval <- max(pred_all$value)

install.packages('viridis')
library('viridis')

test.plot <- pred_all %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0, maxval)) + 
  #scale_fill_gradientn(colours = hcl.colors(16, palette = "Light grays"), limits = c(0,maxval)) +
  facet_grid(covtype2 ~ quantile) +
  geom_point(data = detectors_df_all %>% filter(occasions %in% c(18), array_size %in% c("3x3")), inherit.aes = T,
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
        legend.position="none", legend.key.width = unit(3, "cm"),
        legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

test.plot

## If we go with this quantile/credible interval idea as above -- even though we don't really have much uncertainty in our estimates, unless we want plots to the left and the right that are completely populated by -'s, we end up with the misleading plots we see above (where certain spots randomly have very high values, as we only really sample 0's and 1's in each pixel with the given number of sampling occasions).

# So, it might be better to instead create a map using the CV values where if we have a count of 0 activity centres of 0 99.9% of the time, we just colour the corresponding pixels with a CV of 0, as we are fairly certain of our estimate. CV is used to quantify uncertainty in our estimates, so this seems reasonable to do? 

## ---------------------------------------------------------------------------------------
############################# Uncertainty plots for Figure 4 #### ########################
## ---------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------
############# Creating the objects needed for uncertainty plots for Figure 4 #############
## ---------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------
# Creating the objects we specifically need to create uncertainty plots for row 1 of Fig 4
## ---------------------------------------------------------------------------------------

## Loading in the MCMC results for Row 1 of Figure 4 (note that as stated in 'bayesian_code/Plots_Code.R', we don't need to discard any burn-in, this has been done for us)
load("MCMC_Results/Figure4/HomPP_18occ.RData")
load("MCMC_Results/Figure4/HomPP_52occ.RData")
load("MCMC_Results/Figure4/HomPP_111occ.RData")

## Function to calculate CV for each pixel in an RACD map
# Where 99% or more of the sampled number of activity centres in a pixel is 0, we will set the CV for the pixel to 0. This seems reasonable -- the CV measures the uncertainty in our estimate of the number of activity centres in a pixel, and we are fairly certain that 0 activity centres lie in such a pixel (since we have 10,000 sampled values for the number of activity centres in each pixel)
# This function is very similar to no.movement.density.vector() (uses the same arguments)
cv.values.racd <- function(xlim, ylim, results, M) {

  ## Points at which local density will be estimated
  xg <- seq(xlim[1], xlim[2], by=1)
  yg <- seq(ylim[1], ylim[2], by=1)

  ## Extracting z-values
  # Names of variables that have been monitored
  names <- names(results[1,])
  # Extracting "z" values from MCMC results
  Z <- results[,grep("z", names)]

  ## Extracting activity centres
  # Extracting "s" values from MCMC results (i.e. extracting sampled activity centres)
  S <- results[,grep("s[^i]", names)]
  # x-coordinates of all activity centres
  Sx <- S[,1:M]
  # y-coordinates of all activity centres
  Sy <- S[,-(1:M)]

  ## For each MCMC iteration, storing the number of animals alive and with their activity centres in each cell -- building up posterior distribution of number of activity centres in each pixel
  # Number of pixel centres
  npix <-  (length(xg) - 1) * (length(yg) - 1)
  Dn.vals <-  matrix(0, nrow=nrow(results), ncol=npix)
  for (i in 1:nrow(results)) {
    if ((i %% 100) == 0) print(i) # Track progress
    Sxout <-  Sx[i,][Z[i,] == 1]
    Sxout <-  cut(Sxout, breaks=xg, include.lowest=TRUE)
    Syout <-  Sy[i,][Z[i,] == 1]
    Syout <-  cut(Syout, breaks=yg, include.lowest=TRUE)
    Dn.vals[i,] <-  as.vector(table(Sxout, Syout))
  }

  ## Posterior mean for number of activity centres in each pixel
  posterior.mean <- apply(Dn.vals, 2, mean)
  ## Posterior standard deviation
  standard.deviation <- apply(Dn.vals, 2, sd)

  #browser()
  ## Calculating CV. If 99% or more of a posterior sample for one cell is equal to 0, setting the CV to 0
  cv.values <- standard.deviation/posterior.mean
  for (i in 1:ncol(Dn.vals)) {
    if (mean(Dn.vals[,i]==0) >= 0.99)  {
      cv.values[i] <- 0
    }
  }
  # Returning the CV values
  cv.values
}

## Finding the CV values for the pixels in each plot in the first row of Figure 4
cv.18occ.racd <- cv.values.racd(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.18occ, M=300)
cv.52occ.racd <- cv.values.racd(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.52occ, M=300)
cv.111occ.racd <- cv.values.racd(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.111occ, M=300)

## ---------------------------------------------------------------------------------------
# Creating the objects we specifically need to create uncertainty plots for row 2 of Fig 4
## ---------------------------------------------------------------------------------------

## Loading in the MCMC results for Row 2 of Figure 4 (note that as stated in 'bayesian_code/Plots_Code.R', we need to discard 1000 iterations of burn-in from these objects)
load("MCMC_Results/Figure4/InhomPP_18occ.RData")
load("MCMC_Results/Figure4/InhomPP_52occ.RData")
load("MCMC_Results/Figure4/InhomPP_111occ.RData")
inhom.results.18occ <- inhom.results.18occ[-c(1:1000),]
inhom.results.52occ <- inhom.results.52occ[-c(1:1000),]
inhom.results.111occ <- inhom.results.111occ[-c(1:1000),]

## Covariate used in these plots (see 'bayesian_code/Plots_Code.R' for more details on the code below)
mona.densities <-  small_blurry_mona_df[,c("x", "y", "Dblur")]
split <-  split(mona.densities, mona.densities$y)
mona.densities <-  do.call("rbind", split)
rownames(mona.densities) = NULL
dblur <-  mona.densities[,"Dblur"]
log.dblur <-  log(dblur)

## Function to calculate CV for each pixel in an EACD map. 
# This function is very similar to eacd.density.vector() (uses the same arguments). 
cv.values.eacd <- function(results, covariate, nPix) {
  posterior.mean <- numeric() # Initialising object that will contain posterior mean values for density in each pixel
  posterior.sd <- numeric() # Initialising object that will contain posterior standard deviation values 
  # For loop to calculate this posterior mean and the posterior standard deviation for the animal density in each pixel
  for (i in 1:nPix) {
    if(i%%100==0) print(i) # Track progress
    # Posterior distribution for the density of the ith pixel
    density.posterior <- exp(results[,'beta0'] + results[,'beta1'] * (covariate[i]))
    # Posterior mean of this density
    posterior.mean[i] <- mean(density.posterior)
    # Posterior standard deviation of this density
    posterior.sd[i]  <- sd(density.posterior)
  }

  # Vector of CV values for the EACD map
  cv.values <- posterior.sd/posterior.mean
}

## Finding the CV values for the pixels in each plot in the second row of Figure 4
cv.18occ.eacd <- cv.values.eacd(results=inhom.results.18occ, covariate=log.dblur, nPix=2500)
cv.52occ.eacd <- cv.values.eacd(results=inhom.results.52occ, covariate=log.dblur, nPix=2500)
cv.111occ.eacd <- cv.values.eacd(results=inhom.results.111occ, covariate=log.dblur, nPix=2500)

## ---------------------------------------------------------------------------------------
############# Creating the objects needed for uncertainty plots for Figure 5 #############
## ---------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------
# Creating the objects we specifically need to create uncertainty plots for row 1 of Fig 5
## ---------------------------------------------------------------------------------------

## Loading in the MCMC results for Row 1 of Figure 5 (no burn-in needs to be discarded)
load("MCMC_Results/Figure5/HomPP_7occ.RData")
load("MCMC_Results/Figure5/HomPP_25occ.RData")
load("MCMC_Results/Figure5/HomPP_55occ.RData")

## Finding the CV values for each plot in the first row of Figure 5 (where once again, if at least 99% of the sampled number of activity centres for a pixel is 0, the CV for that pixel is set to 0)
cv.7occ.racd <- cv.values.racd(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.7occ, M=300)
cv.25occ.racd <- cv.values.racd(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.25occ, M=300)
cv.55occ.racd <- cv.values.racd(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.55occ, M=300)

## ---------------------------------------------------------------------------------------
# Creating the objects we specifically need to create uncertainty plots for row 2 of Fig 5
## ---------------------------------------------------------------------------------------

## Loading in the MCMC results for Row 2 of Figure 5 (need to discard 1000 iterations as burn-in)
load("MCMC_Results/Figure5/InhomPP_7occ.RData")
load("MCMC_Results/Figure5/InhomPP_25occ.RData")
load("MCMC_Results/Figure5/InhomPP_55occ.RData")
inhom.results.7occ <- inhom.results.7occ[-c(1:1000),]
inhom.results.25occ <- inhom.results.25occ[-c(1:1000),]
inhom.results.55occ <- inhom.results.55occ[-c(1:1000),]

## Covariate used in plots for row 2 of Figure 5 is same as the covariate used for row 2 of Figure 4, so is given by:
log.dblur # As seen above!

## Finding the CV values for each plot in the second row of Figure 5
cv.7occ.eacd <- cv.values.eacd(results=inhom.results.7occ, covariate=log.dblur, nPix=2500)
cv.25occ.eacd <- cv.values.eacd(results=inhom.results.25occ, covariate=log.dblur, nPix=2500)
cv.55occ.eacd <- cv.values.eacd(results=inhom.results.55occ, covariate=log.dblur, nPix=2500)

## ---------------------------------------------------------------------------------------
############# Creating the objects needed for both sets of uncertainty plots #############
## ---------------------------------------------------------------------------------------

##### Creating a data frame that contains all of the info required for all uncertainty plots #####

## Creating a data frame, labelled 'cv_values_all' that summarises the information we use to create all of the uncertainty plots (for Figures 4 and 5)

## Function to summarise the info for the uncertainty plots for the RACD maps
# The 'nocc' argument is the number of sampling occasions, and 'fig' is the number of the figure we are working with
cv.racd.summary <- function(nocc, fig) {
  # Pixel centres that we are working with
  pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)
  # Name of array we are working with -- if Figure 4, the name is '3x3' and if Figure 5, the name is '7x7'
  if (fig==4) {
    array <- "3x3"
  } else {
    if (fig==5) {
      array <- "7x7"
    }
  }
  # Obtaining the CV values for the uncertainty plot
  cv.vals <- get(paste0("cv.", nocc, "occ.racd"))
  # Data frame of information we want
  dat <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~1", 2500), occasions=rep(nocc, 2500), array_size=rep(array, 2500), value=cv.vals)
  # Returning this data frame
  dat
}

## Function to summarise the info for the uncertainty plots for the EACD maps
cv.eacd.summary <- function(nocc, fig) {
  pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)
  if (fig==4) {
    array <- "3x3"
  } else {
    if (fig==5) {
      array <- "7x7"
    }
  }
  cv.vals <- get(paste0("cv.", nocc, "occ.eacd"))   
  dat <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~log(Dblur)", 2500), occasions=rep(nocc, 2500), array_size=rep(array, 2500), value=cv.vals)
  dat
}

## Function to create the 'cv_values_all' data frame described above
# Here, 'nocc' is the vector of sampling occasions we are working with; 'fig' is the corresponding figure number for each uncertainty plot; 'type' is the corresponding type of map for which we want to create an uncertainty plot (enter as 'RACD' or 'EACD')
cv.overall.summary <- function(nocc, fig, type) {
   # Initialising data frame
  dat <- data.frame()
  for (i in 1:length(nocc)) {
    # If type="RACD", creating a data frame summarising the info for the RACD map corresponding to the given number of sampling occasions and given figure number
    if (type[i]=="RACD") {
      dat.add <- cv.racd.summary(nocc=nocc[i], fig=fig[i])
    } else {
      # If type="EACD", creating a data frame summarising the info for the EACD map
      if (type[i]=="EACD") {
        dat.add <- cv.eacd.summary(nocc=nocc[i], fig=fig[i])
      }
    }
    # Adding the data frame created by racd.summary() or eacd.summary() to our 'dat' data frame
    dat <- rbind(dat, dat.add)
  }
  # Returning the final data frame
  dat
}

## Creating the 'cv_values_all' data frame
nocc <- c(rep(c(18, 52, 111), 2), rep(c(7, 25, 55), 2))
fig <- c(rep(4, 6), rep(5, 6))
type <- c(rep(c(rep("RACD", 3), rep("EACD", 3)), 2))
cv_values_all <- cv.overall.summary(nocc=nocc, fig=fig, type=type)

##### Object that contains information about the detectors in both figures #####

detectors_df_all <- res_acd %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

##### Defining the max value of the colour scale #####

## We want this value to be the same for both sets of uncertainty plots. 

nn <- 3 # Number of simulated datasets we use in each figure  
xx <- cv_values_all %>% filter(array_size == "3x3", occasions %in% capthists_few_alloccs_3x3$noccasions[1:nn])
maxval1 <- max(xx$value) # Max density value found across all uncertainty plots for Figure 4
xx <- cv_values_all %>% filter(array_size == "7x7", occasions %in% capthists_few_alloccs_7x7$noccasions[1:nn])
maxval2 <- max(xx$value) # Max density value found across all uncertainty plots for Figure 5
maxval <- max(maxval1, maxval2) # Using the higher of these two max values as the top value of the colour scale for plots in both figures

## ---------------------------------------------------------------------------------------
####################### Creating the uncertainty plots for Figure 4 ######################
## ---------------------------------------------------------------------------------------

##### Adding the column and row labels used in Figure 4 to the 'cv_values_all' and 'detectors_df_all' objects #####

nn <- 3 # Number of different simulated datasets used in Figure 4
occ <- capthists_few_alloccs_3x3$noccasions # Number of sampling occasions for each dataset
asz <- c("3x3")
chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_3x3$capthist, summary, terse = TRUE)))
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist() # Column lables for Figure 4 -- we'll use these for the uncertainty plots, too!

## Adding the column labels from Figure 4 to the 'cv_values_all' and 'detectors_df_all' objects
cv_values_all$occasions2 <- factor(cv_values_all$occasions, 
                                   levels = occ,
                                   labels = capthist_labels)
detectors_df_all$occasions2 <- factor(detectors_df_all$occasions, 
                                      levels = occ,
                                      labels = capthist_labels)

## Adding the row labels from Figure 4 to the 'predicted_densities_all' and 'detectors_df_all' objects
cv_values_all$covtype2 <- factor(cv_values_all$covtype, 
                                 levels = unique(cv_values_all$covtype),
                                 labels = c("Realised AC", "Expected AC"))
detectors_df_all$covtype2 <- factor(detectors_df_all$covtype, 
                                    levels = unique(detectors_df_all$covtype),
                                    labels = c("Realised AC", "Expected AC"))

##### Creating and saving the uncertainty plots for Figure 4 #####

p2a <- cv_values_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0, maxval), breaks = c(0, 5, 10)) + 
  #scale_fill_gradientn(colours = hcl.colors(16, palette = "Light grays"), limits = c(0,maxval)) +
  facet_grid(covtype2 ~ occasions2) +
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

p2a

## Saving the uncertainty plots for Figure 4
ggsave("mona_3x3_uncertainty.png", p2a, width=8, height=6, dpi=600, bg="white")

## ---------------------------------------------------------------------------------------
####################### Creating the uncertainty plots for Figure 5 ######################
## ---------------------------------------------------------------------------------------

##### Adding the column and row labels used in Figure 5 to the 'cv_values_all' and 'detectors_df_all' objects #####

nn <- 3
occ2 <-capthists_few_alloccs_7x7$noccasions
asz2 <- c("7x7")
chs2 <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_7x7$capthist, summary, terse = TRUE)))
chs2 <- chs2 %>% dplyr::filter(Occasions %in% occ2)
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels2 <- map2(.x = chs2$Detections, .y = chs2$Animals, .f = paster) %>% unlist()

## Adding the column labels from Figure 5 to the 'predicted_densities_all' and 'detectors_df_all' objects
cv_values_all$occasions2 <- factor(cv_values_all$occasions, 
                                  levels = occ2,
                                  labels = capthist_labels2)
detectors_df_all$occasions2 <- factor(detectors_df_all$occasions, 
                                      levels = occ2,
                                      labels = capthist_labels2)

## Adding the row labels from Figure 5 to the 'predicted_densities_all' and 'detectors_df_all' objects
cv_values_all$covtype2 <- factor(cv_values_all$covtype, 
                                 levels = unique(cv_values_all$covtype),
                                 labels = c("Realised AC", "Expected AC"))
detectors_df_all$covtype2 <- factor(detectors_df_all$covtype, 
                                    levels = unique(detectors_df_all$covtype),
                                    labels = c("Realised AC", "Expected AC"))

##### Creating and saving the uncertainty plots for Figure 5 #####

p2b <- cv_values_all %>%
  filter(occasions %in% occ2[1:nn], array_size %in% asz2) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0, maxval), breaks = c(0, 5, 10)) + 
  #scale_fill_gradientn(colours = hcl.colors(16, palette = "Light grays"), limits = c(0,maxval)) +
  facet_grid(covtype2 ~ occasions2) +
  geom_point(data = detectors_df_all %>% filter(occasions %in% occ2[1:nn], array_size %in% asz2), inherit.aes = T,
             colour = "gray80", pch = 4, size = 2) +
  geom_point(data = simulated_points, inherit.aes = F, aes(x=x,y=y),
             colour = "darkorange", pch = 16, size = 1, alpha = 0.5) +
  theme_bw(base_size = 14) +
  coord_equal() +
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

p2b

## Saving the uncertainty plots for Figure 5
ggsave("mona_7x7_uncertainty.png", p2b, width=8, height=6, dpi=600, bg="white")

## Note that in some of our uncertainty plots, we see an inversion of colours -- we see bright colours in areas that correspond to low density in the RACD/EACD maps. This is likely because in these pixels, the posterior mean that we work with is quite low (close to 0). So, even though the associated posterior standard deviation is likely to be fairly small, the fact that it isn't as close to 0 as the posterior mean we are working with means that the CV becomes inflated. Therefore, the CV may be slightly misleading as this means that at least some of the 'hotspots' we see are likely to be due to the fact that the posterior mean in those cells is close to 0, rather than due to the fact that the uncertainty associated with the estimated number of activity centres in these cells is especially high compared to other cells in the plot (which is what we interpret the CV as showing -- the CV tries to help us understand where uncertainty in estimates is high). 
