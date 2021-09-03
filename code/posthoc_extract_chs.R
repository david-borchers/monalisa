load("../output/mona_raw_outputs_100sim.RData")

## Figure 6
# uses fig34_results_100sim
# results are averaged over models fitted to 100 simulated capture histories
# each capture history uses the big animal population (7451 ACs)
# fig34_results_100sim[[i]] is list length 8: 1-4 = 1 occasion (not used), 5-8 = 20 occasions

# capture histories for Figure 6
ch6a <- list(); ch6b <- list(); ch6c <- list(); ch6d <- list()
for(i in 1:100){ch6a[[i]] <- fig34_results_100sim[[i]][[5]]$capture_history} # 6b
for(i in 1:100){ch6b[[i]] <- fig34_results_100sim[[i]][[6]]$capture_history} # 6c
for(i in 1:100){ch6c[[i]] <- fig34_results_100sim[[i]][[7]]$capture_history} # 6d
for(i in 1:100){ch6d[[i]] <- fig34_results_100sim[[i]][[8]]$capture_history} # 6e

## Figure 7
# uses fig5_results_100sim
# results are for a single survey (didn't average over 100 sims cos takes long and estimated surfaces are based on covariate
# coeffients that will hardly change between surveys)
# each capture history uses the big animal population (7451 ACs)
# fig5_results is list length 12: 1-6 = bottom left array, 7-12 = top right array, 1-3, 7-9 = 1 occasion, 4-6,
# 10-12 = 20 occasions, 1,4,7,10 = good cov, 2,5,8,11 = moderate, 3,6,9,12 = Leonardo (not used)

# capture histories for Figure 7
ch7b <- fig5_results[[10]]$capture_history
ch7c <- fig5_results[[4]]$capture_history
ch7e <- fig5_results[[11]]$capture_history
ch7f <- fig5_results[[5]]$capture_history

# Figure 8
# uses fig67_results_100sim
# results are averaged over models fitted to 100 simulated capture histories
# each capture history uses the small animal population (84 ACs)
# fig67_results_100sim[[1]] is list length 12: 1-3 = 1 occasion, 4-6 = 3 occasions, 7-9 = 10 occasions, 10-12 = 20 occasions,
# 1,4,7,10 = D~1, 2,5,8,11 = D~good, 3,6,9,12 = D~moderate)

# capture histories for Figure 8
# top row
# capture histories for Figure 6
ch8a <- list(); ch8b <- list(); ch8c <- list()
for(i in 1:100){ch8a[[i]] <- fig67_results_100sim[[i]][[4]]$capture_history}
for(i in 1:100){ch8b[[i]] <- fig67_results_100sim[[i]][[7]]$capture_history}
for(i in 1:100){ch8c[[i]] <- fig67_results_100sim[[i]][[10]]$capture_history}
# 2nd row
ch8d <- list(); ch8e <- list(); ch8f <- list()
for(i in 1:100){ch8d[[i]] <- fig67_results_100sim[[i]][[5]]$capture_history}
for(i in 1:100){ch8e[[i]] <- fig67_results_100sim[[i]][[8]]$capture_history}
for(i in 1:100){ch8f[[i]] <- fig67_results_100sim[[i]][[11]]$capture_history}
# 3rd row
ch8g <- list(); ch8h <- list(); ch8i <- list()
for(i in 1:100){ch8g[[i]] <- fig67_results_100sim[[i]][[6]]$capture_history}
for(i in 1:100){ch8h[[i]] <- fig67_results_100sim[[i]][[9]]$capture_history}
for(i in 1:100){ch8i[[i]] <- fig67_results_100sim[[i]][[12]]$capture_history}

# Figure 9
# uses fig67_results_100sim
# results are for a single survey
# each capture history uses the small animal population (84 ACs)
# capture histories for Figure 9 use the FIRST simulated capture histories for each of those used in Figure 8
# i.e. ch8a[1], ch8b[1] etc

