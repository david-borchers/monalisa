library(secr)

load("output/mona_inputs.RData")
load("output/simulated_densities.Rdata")

mlmesh <- read.mask(data = mona_df)
simulated_points = simulated_points_lots

# make grid 
detectors <- make.grid(nx = 3, ny = 4, spacex = 4, spacey = 4,
                       originxy = c(15, 15), detector = "count")

## Generate capture histories
capture_history <- sim.capthist(detectors, popn = simulated_points, detectfn = "HHN", 
                                              detectpar = list(lambda0 = 0.69, sigma = 2), 
                                              noccasions = 1, 
                                              nsessions = 1,
                                              seed = 123)

# fit model to log of good covariate
cfit <- secr.fit(capture_history, 
                 model = list(D~log(Dgood_bigD)),  
                 mask = mlmesh, detectfn = "HHN", 
                 start = list(D = exp(0), lambda0 = 0.69, sigma = 2))

# fit model to log of bad covariate
cfit2 <- secr.fit(capture_history, 
                 model = list(D~log(Dblur_bigD)),  
                 mask = mlmesh, detectfn = "HHN", 
                 start = list(D = exp(0), lambda0 = 0.69, sigma = 2))

coef(cfit)
coef(cfit2)

# plotting log(trueD) against log(blurred image), with lines y=x and best fit overlaid
par(mfrow = c(2,2))
m0 <- lm(log(D_bigD) ~ log(Dgood_bigD), data = mona_df)
summary(m0)
plot(log(mona_df$Dgood_bigD), log(mona_df$D_bigD), xlim = c(5,12), ylim = c(5,12))
abline(0,1, col = "red")
abline(coef(m0)[1],coef(m0)[2], col = "green")

m1 <- lm(log(D_bigD) ~ log(Dblur_bigD), data = mona_df)
summary(m1)
plot(log(mona_df$Dblur_bigD), log(mona_df$D_bigD), xlim = c(5,12), ylim = c(5,12))
abline(0,1, col = "red")
abline(coef(m1)[1],coef(m1)[2], col = "green")

# plotting log(trueD) against log(blurred image), with lines y=x and best fit overlaid
m0 <- lm(log(D_bigD) ~ (Dgood_bigD), data = mona_df)
summary(m0)
plot((mona_df$Dgood_bigD), log(mona_df$D_bigD))
abline(0,1, col = "red")
abline(coef(m0)[1],coef(m0)[2], col = "green")

m1 <- lm(log(D_bigD) ~ (Dblur_bigD), data = mona_df)
summary(m1)
plot((mona_df$Dblur_bigD), log(mona_df$D_bigD))
abline(0,1, col = "red")
abline(coef(m1)[1],coef(m1)[2], col = "green")

par(mfrow = c(2,3))
hist(mona_df$D_bigD, breaks = seq(0,90000,5000), main = "Dorig")
hist(mona_df$Dgood_bigD, breaks = seq(0,90000,5000), main = "Dgood")
hist(mona_df$Dblur_bigD,breaks = seq(0,90000,5000), main = "Dblur")

hist(log(mona_df$D_bigD), breaks = seq(4,12,0.5), main = "log(Dorig)")
hist(log(mona_df$Dgood_bigD), breaks = seq(4,12,0.5), main = "log(Dgood)")
hist(log(mona_df$Dblur_bigD), breaks = seq(4,12,0.5), main = "log(Dblur)")

# comparing coefficients of lm and scr fits
coef(m0)
coef(cfit)

coef(m1)
coef(cfit2)
