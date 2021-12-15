## Functions to create RUD maps using sampled values for captured animals only (i.e. not using sampled values associated with unobserved 
## animals in the superpopulation)


## Function to obtain activity centre matrices for observed animals only. Like in 'RUDMaps_Functions.R', this function will produce a list,
## where each element is a matrix of sampled activity centres (where the sampled activity centre on an MCMC iteration where an animal didn't
## 'exist'/had a z-value of 0 is set to (0,0)). However, this list will only contain one matrix for each animal that
## we observed in the given number of sampling occasions/in the given set of capture histories. 

# Below, 'results' is a set of MCMC samples generated using run.MCMC(); 'indices' is the row indices for the encounter data matrix provided to run.MCMC()
# once all-0 rows have been removed from the matrix. In other words, 'indices' should be a set of integers starting at 1 and ending at the
# number of observed animals in the given dataset. 

captured.centres <- function(results, indices) {
  ## Collecting all of the sampled MCMC activity centres for the observed animals only
  # Creating an empty list to store these activity centres in
  list.ans <- list()
  # Creating and storing a matrix in a separate element of the list created above; each matrix contains ALL of the sampled activity centres for the animals 
  # we have observed
  for (i in 1:length(indices)) {
    s1 <- paste("s[", indices[i], ",1]", sep="")
    s2 <- paste("s[", indices[i], ",2]", sep="")
    act.cent <- cbind(results[,s1], results[,s2])
    act.cent <- as.matrix(act.cent)
    list.ans[[i]] <- act.cent
  }
  
  ## Obtaining all of the relevent z-values from the MCMC results (i.e. z-values for the observed animals only)
  # Names of variables that have been monitored
  names <- names(results[1,])
  # Extracting all "z" values from MCMC results
  z.values <- results[,grep("z", names)]
  # Extracting the z values for the animals we have observed 
  z.values <- z.values[,indices]
  
  ## Multiplying each activity centre matrix stored in our list above by the corresponding set of 'z.values' for each animal.
  ## This means that any activity centres associated with a z-value of 0 will be set to (0,0). 
  # Creating an empty list to store the result
  new.list <- list()
  # Using a for loop to do the multiplication described above 
  for (i in 1:length(indices)) {
    multiply <- list.ans[[i]] * z.values[,i]
    new.list[[i]] <- multiply
  }
  
  ## Printing the result 
  new.list
}



## Function to obtain an array of the sampled z-values for the observed animals only. As with the 'extract.z.values()' function (from 'RUDMaps_Functions.R'),
## the resulting z-values for a given iteration can be found across the rows, and the sampled z-values for each (observed) animal can be found down the columns.
## So, the number of rows in the resulting array will be the same as the number of MCMC iterations we ran, and the number of columns will be equal to the number 
## of observed animals. 
captured.z.values <- function(results, indices) {
  ## Obtaining all of the z values from the MCMC results (for observed and unobserved animals)
  # Names of variables that have been monitored
  names <- names(results[1,])
  # Extracting all "z" values from MCMC results
  z.values <- results[,grep("z", names)]
  ## Extracting z values for the animals we know exist/that we have observed
  z.values <- z.values[,indices]
  ## Printing the result
  z.values
}


