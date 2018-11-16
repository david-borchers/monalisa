library(raster)
library(fields)
library(mvtnorm)


### MODEL IS CONTINUOUS TIME BUT TIME-DEPENDENCE IS SPECIFIED USING A SINGLE, DISCRETE-VALUED COVARIATE (NITETIME VS DAYTIME)


### Simulate individuals distributed over rectangular region S

simData = function(beta=c(log(4), 0.8), alpha=c(log(0.5), 0.5), sigma=0.4, xi=-1) {
    ## beta = parameters for computing intensity of point process as a function of covariates
    ## alpha = parameters for computing baseline capture rate per trap as a function of trap-specific covariates
    ## sigma = std deviation of Gaussian kernel for individual movements
    ## xi = effect of time covariate on baseline capture rate per trap


## ... define the region S
s.xmin = -2.5
s.xmax =  2.5
s.ymin = -2.5
s.ymax =  2.5

s.area =  (s.xmax-s.xmin)*(s.ymax-s.ymin)

npixels.x = 1000
npixels.y = 1000

s = raster(ncol=npixels.x, nrow=npixels.y, xmn=s.xmin, xmx=s.xmax, ymn=s.ymin, ymx=s.ymax)
s.loc = xyFromCell(s, 1:ncell(s))


## ... simulate covariate values over discretization of region S
mu1 = c(0.75, -0.2)
Sigma1 = matrix(c(0.25, 0.25, 0.25, 1.00), nrow=2)

mu2 = c(-0.7, 0.6)
Sigma2 = matrix(c(1.0, -0.2, -0.2, 0.25), nrow=2)

mu3 = c(-1.5, -1.0)
Sigma3 = matrix(c(0.15, 0.1, 0.1, 0.25), nrow=2)


mu4 = c(1.75, 1.75)
Sigma4 = matrix(c(0.15, -0.1, -0.1, 0.25), nrow=2)


xcov = 0.25 * dmvnorm(s.loc, mean=mu1, sigma=Sigma1) + 0.45 * dmvnorm(s.loc, mean=mu2, sigma=Sigma2) + 0.15 * dmvnorm(s.loc, mean=mu3, sigma=Sigma3) +  0.15 * dmvnorm(s.loc, mean=mu4, sigma=Sigma4)
xcov = (xcov - mean(xcov))/sd(xcov)

values(s) = xcov
names(s) = 'x'


## ... compute expected density of individuals over discretization of region S
beta.param = beta
X = cbind(rep(1, length(xcov)), xcov)
temp = raster(s)
values(temp) = exp(X %*% beta.param)
names(temp) = 'lambda'
s = addLayer(s, temp)

## ... simulate point pattern of individuals over discretization of S
maxlambda = max(values(s)[,'lambda'])
N.hpp = rpois(1, maxlambda*s.area)

ind.hpp = sample(1:ncell(s), size=N.hpp, replace=FALSE)   #  sampling w/o replacement ensures only 1 indiv per pixel
loc.hpp = s.loc[ind.hpp, ]
lambda.hpp = values(s)[,'lambda'][ind.hpp]

ind.ipp = runif(N.hpp, 0,1) <= lambda.hpp/maxlambda
N.ipp = sum(ind.ipp)
loc.ipp = loc.hpp[ind.ipp, ]

### Simulate recaptures of N individuals in a rectangular array of K traps located within region S
    
ds = 1.5
trap.xmin = s.xmin + ds
trap.xmax = s.xmax - ds
trap.ymin = s.ymin + ds
trap.ymax = s.ymax - ds
ntraps.x = 10
ntraps.y = 10
trap.x = seq(trap.xmin, trap.xmax, length=ntraps.x)
trap.y = seq(trap.ymin, trap.ymax, length=ntraps.y)
traploc.x = rep(trap.x, length(trap.y))
traploc.y = rep(trap.y, each=length(trap.x))
traploc = cbind(traploc.x, traploc.y)

K = length(traploc.x)   # no. traps
T.day = rep(30/2, K)  # period of operation for each trap during day time (days)
T.nite = rep(30/2, K)  # period of operation for each trap during nite time (days)

sigma.param = sigma      # home-range size (std deviation of Gaussian kernel)
alpha.param = alpha
wcov = -0.1*traploc[,1] - 0.1*traploc[,2] + 0.01*traploc[,1]*traploc[,2]
wcov.mean = mean(wcov)
wcov = (wcov - wcov.mean)
W = cbind(rep(1,K), wcov)

psi.param = as.vector(exp(W %*% alpha.param))

xi.param = xi         # effect of time covariate on baseline capture rate per trap

    ## ## .... compute field of detection covariate values and baseline detection rates
    ## temp = raster(s)
    ## wcov.temp = -0.1*s.loc[,'x'] - 0.1*s.loc[,'y'] + 0.01*s.loc[,'x']*s.loc[,'y']
    ## values(temp) = (wcov.temp - wcov.mean)
    ## names(temp) = 'wcov'
    ## s = addLayer(s, temp)
    ## temp = raster(s)
    ## values(temp) = exp(alpha.param[1] + alpha.param[2]*values(s)[,'wcov'])
    ## names(temp) = 'psi'
    ## s = addLayer(s, temp)

ymat = matrix(nrow=N.ipp, ncol=K)
dmat = rdist(loc.ipp, traploc)
for (k in 1:K) {
    ymat[,k] = rpois(N.ipp, lambda=psi.param[k] * exp(-0.5*(dmat[,k]/sigma.param)^2) * (T.nite[k] + T.day[k]*exp(xi.param)) )
}

obs = apply(ymat,1,sum)>0
y = ymat[obs, ]
n = nrow(y)

Z = array(dim=c(nrow(y), ncol(y), max(y)))   # array  of time covariate measurement (binary values: daytime=1, nitetime=0)
for (i in 1:n) {
    for (k in 1:K) {
        if (y[i,k] > 0) {
            probNiteDetection = T.nite[k] / (T.nite[k] + T.day[k]*exp(xi.param))
            y.nite = rbinom(1, size=y[i,k], prob=probNiteDetection)
            if (y.nite==0)  Z[i,k, 1:y[i,k]] = 1
            if (y.nite==y[i,k])  Z[i,k, 1:y[i,k]] = 0
            if (y.nite!=0 & y.nite!=y[i,k]) {
                Z[i,k, 1:y.nite] = 0
                Z[i,k, (y.nite+1):y[i,k]] = 1
            }
        }
    }
}



### Establish discrete grid for integrating over region S
    
s.buffer = ds
sgrid.xmin = min(traploc.x) - s.buffer
sgrid.xmax = max(traploc.x) + s.buffer
sgrid.ymin = min(traploc.y) - s.buffer
sgrid.ymax = max(traploc.y) + s.buffer

sgrid.area =  (sgrid.xmax-sgrid.xmin)*(sgrid.ymax-sgrid.ymin)
sgrid.extent = extent(sgrid.xmin, sgrid.xmax, sgrid.ymin, sgrid.ymax)
s.cropped = crop(s, sgrid.extent)
s.cropped = dropLayer(s.cropped, 2:length(names(s.cropped)) )

gridfact = c(20,20)  #  form grid by aggregating every no. cols (and no. rows)
sgrid = aggregate(s.cropped, fact=gridfact, fun=mean)
sgrid.loc = xyFromCell(sgrid, 1:ncell(sgrid))
sgrid.ncell = ncell(sgrid)
sgrid.cellsize = prod(res(sgrid))


# .... compute matrix of distances between each location on discrete grid and each trap location
distmat = rdist(sgrid.loc, traploc)


# .... compute covariates at each location on discrete grid
X = cbind(rep(1, sgrid.ncell), values(sgrid))


    list(sgrid=sgrid, traploc=traploc, y=y, X=X, W=W, Z=Z, T.day=T.day, T.nite=T.nite, N.ipp=N.ipp)


}

