# for some reason the simulated points from secr are shifted a bit, this just puts them 
# onto the same scale as the rest of the plots. Haven't looked into why this happens but
# pretty sure it doesn't affect anything else, just the plotting.

make_simulated_points_plottable <- function(simulated_points){
  
  simulated_points_shifted <- simulated_points %>% mutate(x = x - 0.5, y = y - 0.5)
  
  # extract simulated density surface 
  
  # extract the densities (first need the ggplot object, then extract what we need)
  simulated_densities <- ggplot(simulated_points_shifted, aes(x = x, y = y)) +
    stat_bin_2d(aes(fill = ..density..), binwidth = c(1,1), drop = FALSE) 
  simulated_densities <- ggplot_build(simulated_densities) %>% .[[1]] %>% .[[1]] 
  
  # put into a small df to convert into an image
  simulated_densities_df <- data.frame(x = simulated_densities$xbin, 
                                       y = simulated_densities$ybin, 
                                       value = simulated_densities$density)
  
  return(simulated_densities_df)
}
