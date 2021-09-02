predicted_densities_for_D0 <- function(cfit, detectors, mymesh, captures.only){
  
  n <-  dim(cfit$capthist)[1]
  noccasions <- dim(cfit$capthist)[2]
  
  # values of the density for the estimated location of the activity center of each point, 
  # evaluated at the center of each cell in the mask, for each capture history/individual
  prob_activity_center_per_seen_indiv <- fxi.secr(cfit, i=1:n)
  
  # add up the densities over captured individuals to get expected activity centers
  # per cell for captured individuals, divide by n to get density
  expnumber_activity_center_all_seen_indiv <- purrr::reduce(prob_activity_center_per_seen_indiv, `+`) 
  prob_activity_center_all_seen_indiv <- expnumber_activity_center_all_seen_indiv / n
  
  # values of the density for the estimated location of the activity center of each point, 
  # evaluated at the center of each cell in the mask, for UNSEEN individuals
  
  # we want p(ac is in (x,y)|not detected at any traps)
  # p(ac|nd) = p(nd|ac)p(ac) / p(nd) 
  #
  # p(nd|ac) is given, = 1 - pdot
  # p(ac) = 1 / number of cells
  # p(nd) = p(nd and ac here) + p(nd and ac not here)
  # = p(nd|ac)p(ac) + p(nd|nac)p(nac)
  
  p_nd.ac <- 1 - pdot(mymesh, detectors, detectfn="HHN", detectpar=detectpar(cfit), noccasions = noccasions)
  p_ac <- 1 / nrow(mymesh)
  p_nd.nac <- sum(p_nd.ac * 1/(nrow(mymesh) - 1)) - ((1/(nrow(mymesh)-1)) * p_nd.ac)
  p_nd <- p_nd.ac * p_ac + p_nd.nac * (1 - p_ac)
  
  prob_activity_center_all_unseen_indiv <- p_nd.ac * p_ac / p_nd
  
  # to get expected number of activity centers in each cell from UNSEEN animals we need
  # to know how many animals we missed
  estimated_N <- round(region.N(cfit)["R.N","estimate"])
  unseen_n <- estimated_N - n
  
  # density and expected number of activity centers from ALL animals
  expnumber_activity_center <- n * prob_activity_center_all_seen_indiv +
    unseen_n * prob_activity_center_all_unseen_indiv
  prob_activity_center <- expnumber_activity_center / estimated_N
  
  # # if you only want predicted densities for the CAPTURED animals (captures.only = T)
  # # ... could output both but the "captures only" scenario seems unlikely so won't bother 
  # if(captures.only == T){
  #   prob_activity_center <- prob_activity_center_all_seen_indiv
  #   expnumber_activity_center <- expnumber_activity_center_all_seen_indiv}
  
  # make data frame with actual densities
  predicted_densities <- data.frame(x = mymesh$x, y = mymesh$y, 
                                    prob_ac = prob_activity_center,
                                    expnumber_ac = expnumber_activity_center,
                                    prob_ac_seenonly = prob_activity_center_all_seen_indiv,
                                    expnumber_ac_seenonly = expnumber_activity_center_all_seen_indiv)
  
  return(predicted_densities)
  
}