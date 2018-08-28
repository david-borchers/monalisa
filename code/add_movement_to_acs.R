# function adding movement to the density of activity centers to get a "space use"
# (or some other name) density

# For each activity center, read the estimated activity center density at every 
# point (x, y), and then redistribute the mass at (x, y) across nearby locations 
# according to distance from (x, y) and the estimated HHN detection function. 
# Does this for all activity centers and sums up across centers to get the 
# "space use" density. 

add_movement_to_acs <- function(xt, yt, ac_densities, sigma){
  
  # reads the density at a point (x, y) on the density surface
  lambda0 <- as.numeric(ac_densities[(ac_densities$x == xt & ac_densities$y == yt), "value"])
  
  # smooths this out according to HHN function, basically redistributing the mass at (x, y) according
  # to distance from (x, y)
  dist <- sqrt((ac_densities$x - xt)^2 + (ac_densities$y - yt)^2)
  value_with_movement <- exp(-dist^2 / (2*sigma^2))
  
  # make sure the total mass allocated across cells is the same as the original mass at (x, y)
  value_with_movement <- value_with_movement / sum(value_with_movement) * lambda0
  
  return(value_with_movement)
}