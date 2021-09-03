## Functions that are needed to produce the figures

## Loading the necessary packages
library("mvtnorm")
library("spatstat")


## Function to gather the sampled activity centres for each animal. This function will produce a list, where each element
## is a matrix containing the sampled activity centres for each animal in the superpopulation. The ith row of each
## matrix will contain 0's if the given animal didn't exist (i.e. had a z-value of 0) on the ith MCMC iteration.

# Here, 'results' is the MCMC samples and M is the size of the superpopulation.

activity.matrices <- function(results, M) {
  # Creating an empty list
  list.ans <- list()
  # Using a for loop to create a matrix for each animal in the superpopulation, containing all of the sampled activity centres
  # for each animal. Each matrix is stored in a separate element of the list created above.
  for (i in 1:M) {
    s1 <- paste("s[", i, ", 1]", sep="")
    s2 <- paste("s[", i, ", 2]", sep="")
    act.cent <- cbind(results[,s1], results[,s2])
    act.cent <- as.matrix(act.cent)
    list.ans[[i]] <- act.cent
  }

  ## Obtaining all of the z values from the MCMC results
  # Names of variables that have been monitored
  names <- names(results[1,])
  # Extracting "z" values from MCMC results -- each row will contain z-values for one MCMC iteration; columns will contain z-values for each animal
  z.values <- results[,grep("z", names)]

  ## Multiplying the activity centre matrices by the corresponding column of z-values for each animal - this means that
  ## any activity centres associated with a z-value of 0 will be set to (0, 0).
  # Creating an empty list
  new.list <- list()
  # Using a for loop to multiply the activity centre matrices by the corresponding z-values for each animal (as described above),
  # and storing the resulting matrices in a new list
  for (i in 1:M) {
    multiply <- list.ans[[i]] * z.values[,i]
    new.list[[i]] <- multiply
  }

  ## Printing the result
  new.list
}



## The function below generates a matrix containing all of the z-values from a set of MCMC samples. The result will be
## a matrix, where each row contains all of the z-values sampled in one MCMC iteration. Each column contains all of the z-values
## sampled for one animal in the superpopulation.

# Here, 'results' is a set of MCMC samples

extract.z.values <- function(results) {
  # Names of variables that have been monitored
  names <- names(results[1,])
  # Extracting "z" values from MCMC results
  z.values <- results[,grep("z", names)]
  # Printing the result
  z.values
}



## Function to generate pixel centres across map area

# Here, xrange and yrange give the range of x- and y-coordinates for the map area (i.e. the lower value of xrange
# gives the leftmost x-coordinate of the pixels that are furthest to the left in our map area, and the upper value gives the
# leftmost x-coordinates of pixels that are furthest to the right, and so on). In addition, x.pixels and y.pixels give the number
# of pixels being used in the x- and y-direction.

centres <- function(xrange, yrange, x.pixels, y.pixels) {
  # Creating an object of class 'owin' representing our map area
  window.2 <- owin(xrange=xrange, yrange=yrange)
  # Generating our set of pixel centres
  points <- gridcentres(window.2, x.pixels, y.pixels)
  # Converting the result to a matrix
  centres <- as.matrix(cbind(points$x, points$y))
  # Printing the result
  centres
}


## Function to generate density vectors for RUD maps (realised usage density maps). The vector will contain the final density
## values for each pixel

# Here, 'results' is a set of MCMC samples, 'activity.centres' is an object produced using the activity.matrices() function
# above, 'pixel.centres' is an object produced using the centres() function above; 'z.values' is an object produced using
# extract.z.values() above; 'n.iter' is the number of MCMC iterations run to produce 'results'; xlim and ylim give the range
# of x- and y-coordinates for the map area; ''points' indicates 'M' is the size of the superpopulation

density.vector <- function(results, activity.centres, pixel.centres, z.values, n.iter=10000,
                          xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), M=300) {

  # Pixel centres
  x <- pixel.centres

  ## Running a double-nested loop so we can work over each animal from each MCMC iteration
  # Empty list where each element will be a vector (generated using data from each MCMC iteration, so number of elements will
  # be equal to number of MCMC iterations we ran)
  iter.vectors <- list()
  # Generating a density vector for each MCMC iteration -- the inner loop generates the density values using sampled values for each
  # animal from the ith MCMC iteration, and the outer loop means that this is repeated for each MCMC iteration
  for (i in 1:n.iter) {
    print(i)
    # Sigma value for ith MCMC iteration
    s <- results[,"sigma"][i]
    # Resulting covariance matrix
    cov.matrix <- diag(2) * s^2
    # z vectors from all MCMC iterations
    z.values <- z.values

    # Empty list to store resulting vector in for each animal in superpopulation
    animal.vec <- list()

    # Loop to generate density values using sampled values for all animals in the ith iteration
    for (j in 1:M) {
      # Getting activity centre for jth animal for ith iteration
      mean <- activity.centres[[j]][i,]
      # Density vector (at each pixel centre) for jth animal for ith iteration
      density <- dmvnorm(x=x, mean=mean, sigma=cov.matrix)
      # Multiplying density vector by z value for the jth animal in the ith iteration - this means that if it was unobserved, the density
      # will be 0 for ALL pixels
      z.vector <- z.values[i,j]

      # Storing resulting density vector in list
      animal.vec[[j]] <- density * z.vector
    }

    # Creating a final density vector for the ith iteration where the density values for each pixel are found by
    # summing together all of the vectors produced above
    reduced.vec <- Reduce('+', animal.vec)
    # Putting this vector into the final list
    iter.vectors[[i]] <- reduced.vec
  }

  ## Now, summing across all of the elements of the 'final list'
  final.vector <- Reduce('+', iter.vectors)

  ## Dividing the result by the number of MCMC iterations that we ran
  final.vector <- final.vector/n.iter

  ## Printing result
  final.vector
}
