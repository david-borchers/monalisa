# monalisa
This is not a species distribution model

### Instructions to recreate analysis and plots from paper

- in project folder create subfolders "output" and "mona_results"
- put image "monalist.jpg" in main project folder
- run the 3 data prep scripts, but probably just once (unless you want to simulate new points)
- run **run_simstudy_for_paper.R**
- run all the plotting scripts
- note: may not look identical, not 100% sure I've used the same seeds (which I can figure out in time) or standardisations (which I'm unclear about anyway). 

### Data preparation

**make_hires_mona_jpg.R**: takes the full size, high res version of the Mona Lisa (monalisa.jpg) and crops the image to the one used in subsequent analysis (output/hires_mona.jpg)

**make_mona_inputs.R**: takes hires_mona.jpg. Creates grayscale values from the hi-res image, a low-res version, and covariate surfaces. Saved as "output/mona_inputs.RData"

**make_simulated_activity_centers.R**: generates activity centers by simulating from the density surface given by the low res version of the Mona Lisa. User has to specify the desired number of generated points, in expectation. The current code generates the two sets of simulated activity centers used in the paper, which are saved as "output/simulated_densities_for_paper.Rdata". Note that this doesn't generate simulated capture histories, just the activity centers.

### Analysis 

**run_simstudy_for_paper**: main file for running SECR on simulated data under various trap configurations. Takes as input the covariate surfaces and the simulated activity centers. Calls "run_secr", which simulates capture histories and fits various SECR models. Outputs predicted density surfaces, detector positions, and estimated values of lambda from the fitted detection function. Note: a pattern match looks for a "1" in the SECR formula and, if one is found, it assumes a constant density model, so watch out for covariates with 1's in them! 

**run_secr**: used by **simstudy_for_paper** and explained above.

**predicted_densities_for_D0**: used by **run_secr**. Uses and SECR model object to produce a predicted density surface. If the model is a constant density, the surface is obtained by adding the density (or expected number of activity centers) of seen animals (obtained from the secr function fxi.secr) and unseen animals (see function for details, but essentially obtained from secr outputs and Bayes rule).

**add_movement_to_acs**: adds movement to the density of activity centers to get a "space use" (or some other name) density. For each activity center, read the estimated activity center density at every point (x, y), and then redistribute the mass at (x, y) across nearby locations according to distance from (x, y) and the estimated HHN detection function. Does this for all activity centers and sums up across centers to get the "space use" density. 

### Plotting

**make_simulated_points_plottable**: for some reason the simulated points from secr are shifted a bit up and right (by 0.5), this function just puts them onto the same scale as the rest of the plots. Haven't looked into why this happens but this function ONLY affects the plot of the simulated activity center surface. 

**plotfig2_inputs**: plots current figure 2 in the paper
**plotfig34_torch**: plots current figure 3, 4 in the paper
**plotfig5_covariates**: plots current figure 5 in the paper
**plotfig6_fewacs**: plots current figure 6 in the paper
**plotfig7_spaceuse**: plots current figure 7 in the paper

### Notes

- The *imager* package can cause conflicts with other packages used (both *secr* and *tidyverse* packages). If you get errors clear the workspace and just load the packages required by a particular script (will fix this up later)

- how to standardise plots for visualisation? Plots from more occasions, or with stronger covariates, have a much wider range than other plots, so facetting them with a single scale hides some of the patterns in the other plots. Most of the scripts here standardise each plot so min = 0 and max = 1, or so that the sum = 1 across all cells. 
