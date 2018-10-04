## Some setup.
library(spatstat)
library(mvtnorm)
set.seed(1234)
## Parameter values for intensity surface.
b0 <- -7
b1 <- 0.025
## Standard deviation of Gaussian kernel for uncertainty.
kernel.sigma <- 1.5
## Standard deviation of Guassian utilisation distribution.
movement.sigma <- 10
## Number of rows/columns of pixels in outer region.
w <- 200
## Setting inner and outer limits of plot.
inner.lims <- rbind(c(0, 100), c(0, 100))
outer.lims <- rbind(c(-25, 125), c(-25, 125))
## Gettin unique x- and y-values for outer region.
unique.x <- seq(outer.lims[1, 1], outer.lims[1, 2], length.out = w)
unique.y <- seq(outer.lims[2, 1], outer.lims[2, 2], length.out = w)
## Pixel centre coordinates.
xs <- rep(unique.x, w)
ys <- rep(unique.y, each = w)
n.pixels <- w^2
## Intensities as a vector.
z1.vec <- exp(b0 + b1*xs)

## Simulating a binomial point pattern via rejection.
max.intensity <- max(z1.vec)
n.candidates <- round(max.intensity*w^2)
## Candidate points.
candidates <- cbind(runif(n.candidates, outer.lims[1, 1], outer.lims[1, 2]),
                    runif(n.candidates, outer.lims[2, 1], outer.lims[2, 2]))
## Finding nearest pixel.
candidates.dists <- crossdist(candidates[, 1], candidates[, 2], xs, ys)
nearest.pixel <- apply(candidates.dists, 1, function(x) which(x == min(x))[1])
## Finding acceptance probability.
candidate.prob <- z1.vec[nearest.pixel]/max.intensity
## Keeping the accepted points.
acs <- candidates[runif(n.candidates, 0, 1) < candidate.prob, ]
n.acs <- nrow(acs)

## Gaussian kernel densities for each individual to reflect measurement error.
kd <- matrix(NA, nrow = n.acs, ncol = n.pixels)
for (i in 1:n.acs){
    kd[i, ] <- dmvnorm(cbind(xs, ys), mean = acs[i, ], sigma = kernel.sigma^2*diag(2))
}
z2.vec <- apply(kd, 2, sum)

## Sum of utilisation distributions to reflect animal movement.
ud <- matrix(NA, nrow = n.acs, ncol = n.pixels)
for (i in 1:n.acs){
    ud[i, ] <- dmvnorm(cbind(xs, ys), mean = acs[i, ], sigma = movement.sigma^2*diag(2))
}
z3.vec <- apply(ud, 2, sum)

## Making things into matrices for image().
z1 <- matrix(NA, nrow = w, ncol = w)
z2 <- matrix(NA, nrow = w, ncol = w)
z3 <- matrix(NA, nrow = w, ncol = w)
## Putting everything into a matrix.
for (i in 1:n.pixels){
    ## For each pixel, extract x- and y-coordinates.
    x <- xs[i]
    y <- ys[i]
    index.x <- which(x == unique.x)
    index.y <- which(y == unique.y)
    ## Filling the matrices.
    z1[index.x, index.y] <- z1.vec[i]
    z2[index.x, index.y] <- z2.vec[i]
    z3[index.x, index.y] <- z3.vec[i]
}
## Trimming everything to meet the inner limits.
acs <- acs[acs[, 1] >= inner.lims[1, 1] &
           acs[, 1] <= inner.lims[1, 2] &
           acs[, 2] >= inner.lims[2, 1] &
           acs[, 2] <= inner.lims[2, 2], ]
x.keep <- unique.x >= inner.lims[1, 1] &
    unique.x <= inner.lims[1, 2]
y.keep <- unique.y >= inner.lims[2, 1] &
    unique.y <= inner.lims[2, 2]
unique.x <- unique.x[x.keep]
unique.y <- unique.y[y.keep]
z1 <- z1[x.keep, y.keep]
z2 <- z2[x.keep, y.keep]
z3 <- z3[x.keep, y.keep]

## Making the plot.
png("density-plots.png", width = 620, height = 220)
## Set up main plot with three panels.
par(mfrow = c(1, 3), xaxs = "i", yaxs = "i", mar = rep(0, 4), omi = rep(0.1, 4))
## First panel: intensity surface.
plot.new()
plot.window(xlim = inner.lims[1, ], ylim = inner.lims[2, ], asp = 1)
image(x = unique.x, y = unique.y, z = z1, add = TRUE)
points(acs, pch = 16)
box()
## Second panel: kernel density estimate.
plot.new()
plot.window(xlim = inner.lims[1, ], ylim = inner.lims[2, ], asp = 1)
image(x = unique.x, y = unique.y, z = z2, add = TRUE)
box()
## Third panel: sum of utilisation distributions.
plot.new()
plot.window(xlim = inner.lims[1, ], ylim = inner.lims[2, ], asp = 1)
image(x = unique.x, y = unique.y, z = z3, add = TRUE)
box()
dev.off()
