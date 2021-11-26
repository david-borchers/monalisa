dpoisLocal_normal_2 = nimbleFunction(
  run = function(x = double(1), detNums = double(0), detIndices = double(1), lambda0 = double(0), sigma = double(0), s = double(1), trapCoords = double(2), localTrapsIndices = double(2), localTrapsNum = double(1), resizeFactor = double(0, default = 1), habitatGrid = double(2), indicator = double(0), log = double(0, default = 0)) {
    # Defining the type of value that will be returned (necessary for a nimbleFunction object)
    returnType(double(0)) # Will be returning a scalar value

    ## Now, finding the likelihood value:
    # We are deciding what to print if the individual isn't available for detection (i.e. if z-value is 0)
    # So, if the individual isn't available for detection (z=0) and there are no detections recorded for the individual
    # and we don't want the log-probability, the probability of the observed encounter history for that animal is 1.
    # If we want the log-probability in this case, the returned value is log(1)=0.
    # If the individual isn't available for detection and there ARE detections recorded for the detection, and if we don't want
    # the log-probability, the probability for the observed encounter history is 0 (we shouldn't have any detections).
    # If we want the log-probability in this case, the returned value is log(0)=-Inf
    if (indicator == 0) {
      if (detNums == 0) {
        if (log == 0)
          return(1)
        else return(0)
      }
      else {
        if (log == 0)
          return(0)
        else return(-Inf)
      }
    }

    # habitatGrid is a matrix of indices for each pixel that we are considering
    # And, resizeFactor is a value we specify that aggregates resizeFactor * resizeFactor pixel to represent one pixel (potentially speeding up computation)
    # So here, it looks like we are taking the sampled activity centre (given by s) and using it to find the index in 'habitatGrid' that corresponds to the chosen activity centre?
    sID <- habitatGrid[trunc(s[2]/resizeFactor) + 1, trunc(s[1]/resizeFactor) + 1]
    # And then are identifying the 'local' traps for which we calculate the detection function when finding the probability of a given ncounter history (aren't considering all traps to speed up computation)?
    # That is, we are getting a vector of the numbers/indices of the 'local' traps (so have '10' if the 10th trap is a local trap)
    theseLocalTraps <- localTrapsIndices[sID, 1:localTrapsNum[sID]]

     # If the given animal is detected at least once, going through each detection and if they weren't detected at any of the
    # 'local' traps detected above, are returning a probability of 0 (or a log-prob of -Inf, if we want the log-prob)
    # We aren't considering anything beyond the 'local' traps as we don't think animals would be detected at any of the other traps (for speed of computation),
    # so we set the overall probability of the encounter history to 0 -- we are saying it seems very, very unlikely for an animal to be detected at a non-local trap (probability 0)
    if (detNums > 0) {
      for (r in 1:detNums) {
        if (sum(detIndices[r] == theseLocalTraps) == 0) {
          if (log == 0)
            return(0)
          else return(-Inf)
        }
      }
    }

    # Coefficient for the encounter function
    alpha <- -1/(2 * sigma * sigma)
    # Initialising the log probability
    logProb <- 0
    # nimC() is the NIMBLE equivalent of c()
    # So, we are taking the encounter histories from the sparse matrix (sparse matrix look slightly different to the original encounter matrix, and also contains fewer columns)
    # and concatenating 0 to the end
    # Also note this is different from what we'd get from the 'raw' encounter matrix, as it contains the indices of the traps at which an animal was caught at, rather than the number of times an animal
    # was detected at each trap
    detIndices1 <- nimC(detIndices, 0)
    # Initialising the count
    count <- 1
    # For loop, running for each 'local' trap that we have identified for the given activity centre
    for (r in 1:localTrapsNum[sID]) {
      # Looks like if the index/number of the local trap is present in the sparse encounter history (i.e. if an animal was caught at a 'local' trap), then we:
      # * Find the distance between the activity centre and the given 'local' trap
      # * Then, use alpha and this distance to find the probability of detection at the given trap, using the encounter function
      # * And then adding this probability to the log probability (by considering the Poisson probability of being detected at the given trap the given number of times, where the mean of the Poisson distribution is the expected number of times we'd detect the given animal at the given device, found using the encounter function) -- this is an assumption we make in SCR)
      if (theseLocalTraps[r] == detIndices1[count]) {
        d2 <- pow(trapCoords[theseLocalTraps[r], 1] - s[1],
                  2) + pow(trapCoords[theseLocalTraps[r], 2] -
                             s[2], 2)
        p <- lambda0 * exp(alpha * d2)
        logProb <- logProb + dpois(x[count], lambda = p, log = TRUE)
        count <- count + 1
      }
      # And if an animal wasn't caught at the given local trap, then we find the probability of this occurring (i.e. the probability of the animal not being caught at this local trap),
      # by using the detection function in a similar manner
      else {
        d2 <- pow(trapCoords[theseLocalTraps[r], 1] - s[1],
                  2) + pow(trapCoords[theseLocalTraps[r], 2] -
                             s[2], 2)
        p <- lambda0 * exp(alpha * d2)
        logProb <- logProb + dpois(0, lambda = p, log = TRUE)
      }
    }

     # If we didn't specify wanting a log-probability, exponentiating the result (so we get just a probability)
    if (log) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
  })





