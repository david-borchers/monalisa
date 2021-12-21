# monalisa

Code, data, and model objects used to generate the results in *"That's not the Mona Lisa! How to interpret spatial capture-recapture density surface estimates"* by Ian Durbach, Rishika Chopara, David Borchers, Rachel Phillip, Koustubh Sharma, and Ben Stevenson (in review, Dec 2021).

There are three main sets of analyses/results in the paper: one in the introduction demonstrating the different density surface using a small made-up example; another a simulation study using the Mona Lisa as the true expected activity center; and finally a case study using data from a previously published SCR survey of tiger in the Nagarahole Tiger Sanctuary in India. The Mona Lisa simulations are done using SCR implemented using both maximum likelihood and Bayesian methods, so this repo has four main folders, one for each of the analyses:

- Code for the introductory example is in the **code_intro** folder.
- Code for the Mona Lisa simulations implemented using maximum likelihood methods is in the **code** folder.
- Code for the Mona Lisa simulations implemented using Bayesian methods is in the **bayesian_code** folder.
- Code for the Nagarahole analysis is in the **nagarahole** folder.

Other folders store model output (**output**) and input files for the paper (**paper**). 

**To reproduce the frequentist Mona Lisa analysis**

- to recreate the inputs to the simulation, run *make_mona_inputs.R* and *make_simulated_activity_centers.R*. These convert an image of the Mona Lisa (in *output/hires_mona.jpg*) into a density surface, generate two sets of activity centers, and create a set of potential covariates by blurring the image/density surface as described in the paper. Output is saved to *output/mona_inputs.RData* and *output/simulated_densities.RData* (these appear in the folder already, and running these scripts will overwrite them).
- run *run_sim_study_for_paper.R*. This calls other scripts in the folder, primarily *run_secr.R*, which contains most of the important parts of the simulation. Using multiple cores is recommended if rerunning all simulations. Note that this will overwrite results in *output/mona_raw_outputs.RData*.
- various plotting scripts recreate the plots in the paper: *plot_realised_acd_many.RData*, *plot_expected_acd_many.RData*, *plot_realised_and_expected_acd_few.RData*, *plot_realised_usage_few.RData*. 

**To reproduce the Bayesian Mona Lisa analysis**

- use the simulated capture histories from *output/capthists.RData* to run the MCMC models.
- for Figure 7, go to the **bayesian_code** folder, and run *Figure 7/ch7b.R*, *Figure 7/ch7c.R*, *Figure 7/ch7e.R* and *Figure 7/ch7f.R* to run the MCMC for plots (b), (c), (e) and (f), respectively. Using a machine with a lot of RAM is recommended, and a few hours will be needed to run the MCMC models. Then, run *Figure7.R* to create the plot. 
- for Figure 9, go to the **bayesian_code** folder and run *Figure9.R*. This file contains function calls that will run the MCMC, and then put together the final figure. The MCMC models will be fairly quick to run (even with a standard amount of RAM). 

**To reproduce the Nagarahole analysis**

- run *nagarahole_modelfits.R*. This fits SCR models under the assumption of the different detector arrays described in the paper.
- run *nagarahole_plots.R* to recreate the plots in the paper.
