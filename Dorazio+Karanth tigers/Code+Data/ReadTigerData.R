library(raster)
library(fields)
library(maptools)
library(rgdal)

ReadCameraTrapData = function(resolution, shapefile, trapLocFile, trapOpFile, sunriseANDsunsetFile, detectionFile) {
    
    ## resolution =  resolution (meters) used to create discrete grid of spatial domain
    ## shapefile = name of shape file for spatial domain and spatial covariates
    ## trapLocFile = name of file containing trap locations and trap-specific covariates
    ## trapOpFile = name of file containing times of trap operation (start and stop times)
    ## sunriseANDsunsetFile = name of file containing sunrise and sunset for period of trapping
    ## detectionFile = name of file containing detections of individuals by date, time, and trap

    

    ### Utility function for computing daytime and nightime trapping effort during a single time interval
    ComputeTrapEffort = function(startDate, endDate, startTime, endTime, sunriseANDsunsetData) {

        dateChar = as.character(sunriseANDsunsetData[,'date'])
        date = as.Date(dateChar, format='%Y-%m-%d')
        sunriseChar = paste('0', as.character(sunriseANDsunsetData[,'sunrise']), sep='')
        sunsetChar  = as.character(sunriseANDsunsetData[,'sunset'])
        sunrise = as.numeric(strptime(paste(dateChar, sunriseChar, sep=' '), format='%Y-%m-%d %H%M'))
        sunset = as.numeric(strptime(paste(dateChar, sunsetChar, sep=' '), format='%Y-%m-%d %H%M'))
        ndates = length(date)
        daytime = sunset - sunrise
        nitetime = sunrise[2:ndates] - sunset[1:(ndates-1)]

     
        startTime = as.numeric(strptime(paste(startDate, startTime, sep=' '), format='%d-%b-%y %H:%M'))
        endTime = as.numeric(strptime(paste(endDate, endTime, sep=' '), format='%d-%b-%y %H:%M'))
        startDate = as.Date(startDate, format='%d-%b-%y')
        endDate = as.Date(endDate, format='%d-%b-%y')
        ind.start = match(startDate, date)
        ind.end = match(endDate, date)

        ## ... compute effort on day camera was started
        daytimeEffort = 0
        nitetimeEffort = 0
        if (startTime>=sunrise[ind.start] & startTime<sunset[ind.start]) {  # camera started during daytime
            daytimeEffort = sunset[ind.start] - startTime
            nitetimeEffort = min(nitetime[ind.start], (endTime-sunset[ind.start]))
        }
        else {  ## camera started during nitetime
         nitetimeEffort = min(sunrise[ind.start+1],endTime) - startTime
        }

        ind = ind.start + 1
        while(ind <= ind.end) {
            if (ind < ind.end) {
                daytimeEffort = daytimeEffort + daytime[ind]
                nitetimeEffort = nitetimeEffort + nitetime[ind]
            }
            else {
                if (endTime>=sunrise[ind] & endTime<sunset[ind]) {  # camera stopped during daytime
                    daytimeEffort = daytimeEffort + (endTime-sunrise[ind])
                }
                else {  ## camera stopped during nitetime
                    daytimeEffort = daytimeEffort + daytime[ind]
                    nitetimeEffort = nitetimeEffort + (endTime-sunset[ind])
                    }
            }
            ind = ind + 1
            }

        secondsPerDay = 24 * 60 * 60
        retVal = c(daytimeEffort, nitetimeEffort) / secondsPerDay
        retVal
        }


    ### Read data from files

    ## ... read shape file for spatial domain and create rasterized grid of domain
    spPolygonData=readOGR(dsn=getwd(), layer=shapefile, verbose=FALSE)

    spExtent = as.matrix(extent(spPolygonData))
    n.ColsAndRows = round((spExtent[,2]-spExtent[,1]) / resolution, digits=0)

    shape.grid=raster(ncol=n.ColsAndRows['x'], nrow=n.ColsAndRows['y'])
    extent(shape.grid) = extent(spPolygonData)

    spPolygonNames = names(spPolygonData)
    
    s=rasterize(spPolygonData, shape.grid, spPolygonNames[1])
    names(s) = spPolygonNames[1]
    for (i in 1:length(spPolygonNames)) {
        if (i>1) {
        temp = rasterize(spPolygonData, shape.grid, spPolygonNames[i])
        names(temp) = spPolygonNames[i]
        s = addLayer(s, temp)
        }
    }


    ## ... read sunrise and sunset data
    sunriseANDsunset.df = read.table(sunriseANDsunsetFile, header=TRUE)
    dateChar = paste(sunriseANDsunset.df[,'year'], sunriseANDsunset.df[,'month'], sunriseANDsunset.df[,'day'], sep='-')
    sunriseANDsunset.date = as.Date(dateChar, format='%Y-%m-%d')
    sunriseANDsunsetData = data.frame(date=dateChar, sunriseANDsunset.df[, c('sunrise','sunset')])

    ## ... read trap effort data
    trapData = read.csv(trapOpFile)
    trapDateVec = as.character(trapData[,'Date'])
    trapTimeVec = as.character(trapData[,'Time'])
    trapIDVec = as.character(trapData[, 'LOC_ID'])

    trapID = names(table(trapIDVec))
    ntraps = length(trapID)
    trapTimeMatrix = matrix(nrow=ntraps, ncol=max(table(trapIDVec)))
    for (i in 1:ntraps) {
        ind = (trapIDVec == trapID[i])
        trapTimes = strptime(paste(trapDateVec[ind], trapTimeVec[ind], sep=' '), format='%d-%b-%y %H:%M')
        trapTimeMatrix[i, 1:sum(ind)] = as.numeric(trapTimes)
        }


    trapDaytimeEffort = matrix(nrow=ntraps, ncol=max(table(trapIDVec))/2)
    trapNitetimeEffort = trapDaytimeEffort
    for (i in 1:ntraps) {
        ind = (trapIDVec == trapID[i])
        trapDate = trapDateVec[ind]
        trapTime = trapTimeVec[ind]
        ntimes = sum(ind)
        for (j in 1:(ntimes/2)) {
            effort = ComputeTrapEffort(startDate=trapDate[j*2-1], endDate=trapDate[j*2], startTime=trapTime[j*2-1], endTime=trapTime[j*2], sunriseANDsunsetData)
            trapDaytimeEffort[i,j] = effort[1]
            trapNitetimeEffort[i,j] = effort[2]
        }
    }


    ## ...  read trap location data
    trapLocData =  read.csv('trap_2015_coord.csv')
    trapIDVec = as.character(trapLocData[, 'LOC_ID'])
    trapCoord = as.matrix(trapLocData[, c('X_COORD', 'Y_COORD')])
    traploc = matrix(nrow=ntraps, ncol=2)
    for (i in 1:ntraps) {
        ind = (trapIDVec == trapID[i])
        traploc[i, ] = trapCoord[ind, ]
    }


    ## ... read capture data
    capData = read.csv(detectionFile)
    indivIDVec = as.character(capData[,'ANIMAL_ID'])
    capDateVec = as.character(capData[,'Date'])
    capTimeVec = as.character(capData[,'Time'])
    capTrapIDVec = as.character(capData[, 'LOC_ID'])

    indivID = names(table(indivIDVec))
    nindivs = length(indivID)
    capFreqMatrix = matrix(0, nrow=nindivs, ncol=ntraps)
    for (irow in 1:length(indivIDVec)) {
        i = match(indivIDVec[irow], table=indivID)
        j = match(capTrapIDVec[irow], table=trapID)
        capFreqMatrix[i,j] = capFreqMatrix[i,j] + 1
    }

    capTimeArray = array(dim=c(nrow(capFreqMatrix), ncol(capFreqMatrix), max(capFreqMatrix)))
    capTimeCovariateArray = capTimeArray
    capTimeCounter = matrix(0, nrow=nrow(capFreqMatrix), ncol=ncol(capFreqMatrix))
    for (irow in 1:length(indivIDVec)) {
        i = match(indivIDVec[irow], table=indivID)
        j = match(capTrapIDVec[irow], table=trapID)
        capTime = strptime(paste(capDateVec[irow], capTimeVec[irow], sep=' '), format='%d-%b-%y %I:%M:%S %p')
        capDate = as.Date(capDateVec[irow], format='%d-%b-%y')

        ind = match(capDate, sunriseANDsunset.date)
        sunriseChar = paste('0', as.character(sunriseANDsunsetData[ind,'sunrise']), sep='')
        sunsetChar  = as.character(sunriseANDsunsetData[ind,'sunset'])
        sunrise = as.numeric(strptime(paste(sunriseANDsunsetData[ind,'date'], sunriseChar, sep=' '), format='%Y-%m-%d %H%M'))
        sunset =  as.numeric(strptime(paste(sunriseANDsunsetData[ind,'date'], sunsetChar, sep=' '), format='%Y-%m-%d %H%M'))

        capTimeCovValue = ifelse(capTime>=sunrise & capTime<sunset, 1, 0)
        capTimeCounter[i,j] = capTimeCounter[i,j] + 1
        capTimeArray[i,j, capTimeCounter[i,j]] = as.numeric(capTime)
        capTimeCovariateArray[i,j, capTimeCounter[i,j]] = capTimeCovValue
    }

    ## ... check to ensure that each capture time resides within a trap's operational period
    for (i in 1:nrow(capFreqMatrix)) {
        for (k in 1:ncol(capFreqMatrix)) {
            if (capFreqMatrix[i,k]>0) {
                capTime = capTimeArray[i,k, 1:capFreqMatrix[i,k]]
                trapTime = trapTimeMatrix[k, ]
                J = sum(!is.na(trapTime))
                even = seq(2,J,  by=2)
                odd = seq(1,J-1, by=2)
                for (icap in 1:length(capTime)) {
                    if (!any(capTime[icap]>=trapTime[odd] & capTime[icap]<=trapTime[even])) {
                        stop("STOP: one or more capture times do not reside within a trap's operational period")
                    }
                }
            }
        }
    }

    ## ... check to ensure that capture times occur without multiplicity at each trap
    for (i in 1:nrow(capFreqMatrix)) {
        for (k in 1:ncol(capFreqMatrix)) {
            if (capFreqMatrix[i,k]>0) {
                capTime = capTimeArray[i,k, 1:capFreqMatrix[i,k]]
                if (length(unique(capTime)) < capFreqMatrix[i,k]) {
                    cat("capture times of individual ID ", indivID[i], " at trap ID", trapID[k], " are not unique", "\n")
                }
            }
        }
    }


    list(sgrid=s, traploc=traploc, trapcov=NA, trapDaytimeEffort=rowSums(trapDaytimeEffort, na.rm=TRUE), trapNitetimeEffort=rowSums(trapNitetimeEffort, na.rm=TRUE), captureFreqMatrix=capFreqMatrix, captureTimeCovariateArray=capTimeCovariateArray)
    
    }  # end of ReadCameraTrapData

    
