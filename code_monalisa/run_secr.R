### function for running secr on simulated points

# to assess effect of number of occasions, specify a vector of occasions to simulate over
# (noccasions). run_secr simulates over the maximum of noccasions, and then fits separate
# models after the specified number of occasions e.g. noccasions = c(3,5,10) simulates
# 10 occasions and fits models with (a) the first 3, (b) the first 5, (c) all 10 occasions

run_secr <- function(simulated_points,
                     secr.fitformula, 
                     nx, ny, # number of traps in x and y dim
                     dx, dy, # distance between traps in x and y dim
                     xorig, yorig, # origin of trap array
                     sigma, lambda0, # detection function
                     noccasions, # a VECTOR of occasions to simulate
                     my.seed = sample(1:100000, 1), # in case you want to set own seed
                     captures.only = FALSE){ # set to T if you only want results for captures (D~1 model only)
  
  print("next run")
  
  n_pts <- nrow(simulated_points) 
  
  dx <- dx
  dy <- dy
  
  # make a grid of detectors
  detectors <- make.grid(nx = nx, ny = ny, spacex = dx, spacey = dy,
                         originxy = c(xorig, yorig), detector = "count")
  
  detectors_df <- data.frame(x = as.integer(), 
                             y = as.integer(), 
                             covtype = as.character(),
                             occasions = as.integer(),
                             array_size = as.character(),
                             array_spacing = as.character(),
                             array_origin = as.character(),
                             lambda0 = as.numeric(),
                             sigma = as.numeric(),
                             n_pts = as.integer())
  
  predicted_densities <- data.frame(x = as.integer(), 
                                    y = as.integer(), 
                                    value = as.numeric(),
                                    covtype = as.character(),
                                    occasions = as.integer(),
                                    array_size = as.character(),
                                    array_spacing = as.character(),
                                    array_origin = as.character(),
                                    lambda0 = as.numeric(),
                                    sigma = as.numeric(),
                                    n_pts = as.integer())
  
  estimated_sigma <- data.frame(est_sigma = as.numeric(),
                                covtype = as.character(),
                                occasions = as.integer(),
                                array_size = as.character(),
                                array_spacing = as.character(),
                                array_origin = as.character(),
                                lambda0 = as.numeric(),
                                sigma = as.numeric(),
                                n_pts = as.integer())
  
  ## Generate capture histories
  lambda0 <- lambda0
  sigma <- sigma
  D <- n_pts * 4 #(as 50x50 area is 1/4 of a hectare)
  
  capture_history_max_occasions <- sim.capthist(detectors, popn = simulated_points, detectfn = "HHN", 
                                                detectpar = list(lambda0 = lambda0, sigma = sigma), 
                                                noccasions = max(noccasions), 
                                                nsessions = 1,
                                                seed = my.seed)
  
  for(i in noccasions){
    
    filename <- paste(i, ifelse(i == 1, "occasion", "occasions"))
    
    capture_history <- subset(capture_history_max_occasions, occasions = 1:i)
    
    n <- dim(capture_history)[1]
    
    # fit model specified by secr.fitformula
    cfit <- secr.fit(capture_history, 
                     model = list(as.formula(secr.fitformula)),  
                     mask = mlmesh, detectfn = "HHN", 
                     start = list(D = D, lambda0 = lambda0, sigma = sigma),
                     ncores = 1, method = "Nelder-Mead")
    
    # with covariates predicted densities come straight from the fitted model, 
    # with no covariates need to compute density of activity centers in this realization
    # this is done by the predicted_densities_for_D0 function.
    
    # caution: if there is ANY "1" in secr.fitformula then assumes its a "no covariate" model!
    if (str_detect(secr.fitformula, "1")) { 
      preds <- predicted_densities_for_D0(cfit, detectors, mlmesh, captures.only)$value
    } else {
      preds <- covariates(predictDsurface(cfit))$D.0
    }
    
    # append results
    predicted_densities <- rbind(predicted_densities,
                                 data.frame(x = mlmesh$x, y = mlmesh$y, 
                                            value = preds,
                                            covtype = secr.fitformula,
                                            occasions = i,
                                            array_size = paste0(nx,"x",ny),
                                            array_spacing = paste0(dx,"_",dy),
                                            array_origin = paste0(xorig,"_",yorig),
                                            lambda0 = lambda0,
                                            sigma = sigma,
                                            n_pts = length(mlmesh$x)))
    
    detectors_df <- rbind(detectors_df,
                          data.frame(x = detectors$x, y = detectors$y, 
                                     covtype = secr.fitformula,
                                     occasions = i,
                                     array_size = paste0(nx,"x",ny),
                                     array_spacing = paste0(dx,"_",dy),
                                     array_origin = paste0(xorig,"_",yorig),
                                     lambda0 = lambda0,
                                     sigma = sigma,
                                     n_pts = length(mlmesh$x)))
    
    estimated_sigma <- rbind(estimated_sigma, 
                             data.frame(est_sigma = ifelse(str_detect(secr.fitformula,"1"),
                                                           exp(cfit$fit$par[3]),
                                                           exp(cfit$fit$par[4])),
                                        covtype = secr.fitformula,
                                        occasions = i,
                                        array_size = paste0(nx,"x",ny),
                                        array_spacing = paste0(dx,"_",dy),
                                        array_origin = paste0(xorig,"_",yorig),
                                        lambda0 = lambda0,
                                        sigma = sigma,
                                        n_pts = length(mlmesh$x)))
    
  }
  
  return(list(predicted_densities = predicted_densities, 
              detectors_df = detectors_df,
              estimated_sigma = estimated_sigma))
  
}
