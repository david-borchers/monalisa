# This is code to create the various density plots from a Poisson point process 
# with linear trend, which are used in the Introduction.
# ---------------------------------------------------------------- dlb Oct 2022

library(secr)
library(scrmlebook)
library(sp)

usage.sigma = 5
error.sigma = 2.5
dmax = 60
errmult = 6 # ratio ofmax error to error.sigma
max.error.sigma = errmult*error.sigma
buffer = 30


traps = make.grid()
mask = make.mask(traps,buffer=buffer)
plot(mask)
plot(traps,add=TRUE)
summary(mask)
boundpoly = Polygons(list(Polygon(attributes(mask)$boundingbox)),"boundary")
boundary = SpatialPolygons(list(boundpoly))

tbbox = data.frame(x=c(min(traps$x),max(traps$x),max(traps$x),min(traps$x)),
                   y=c(min(traps$y),min(traps$y),max(traps$y),max(traps$y)))
trapboundpoly  = Polygons(list(Polygon(tbbox)),"trapboundary")
trapboundary = SpatialPolygons(list(trapboundpoly))
traplim = bbox(trapboundary)

D.linear.x = function(mask,Dleft,Dright) {
  xrange = range(mask$x)
  dx = diff(xrange)
  slope = (Dright-Dleft)/dx
  D = Dleft + slope*(mask$x-xrange[1])
}

distance = function(to,from) return(sqrt((to[1]-from[1])^2 + (to[2]-from[2])^2))

addnormal = function(ac,mask,sigma,lambda=1) {
  dists = unlist(apply(mask,1,distance,from=ac))
  dens = lambda*exp(-dists^2/(2*sigma^2))
  dens = dens/sum(dens)
  return(dens)
}


# Creat liner trend ac density surface and add to mask
covariates(mask)$D.0 = D.linear.x(mask,0,100)
# Make density surface a SpatialPixelsDataFrame, and create boundary polygon
spdf.D = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$D.0))
Dbbox = bbox(spdf.D)
Dbdf = data.frame(x=c(Dbbox[1,1],Dbbox[1,2],Dbbox[1,2],Dbbox[1,1]),
                  y=c(Dbbox[2,1],Dbbox[2,1],Dbbox[2,2],Dbbox[2,2]))
Dboundpoly  = Polygons(list(Polygon(Dbdf)),"bufferboundary")
Dboundary = SpatialPolygons(list(Dboundpoly))
# Make SpatialPixelsDataFrame for surface only inside trapboundary
inspdf.D = spdf.D[trapboundary,]

# Generate realisation of the point process
set.seed(1)
pop = sim.popn(covariates(mask)$D.0,mask,model2D="IHP",Ndist="fixed")
N = dim(pop)[1]
# make this population a SpatialPointsDataFrame
spop = SpatialPointsDataFrame(as.matrix(pop),data=data.frame(size=rep(1,N)))
# make the members inside trapboudary a SpatialPointsDataFrame
inpop=spop[trapboundary,]

## Create usage density surface (V slow!) and add do mask
#ausage = rep(0,dim(mask)[1])
#for(i in 1:dim(pop)[1]) ausage = ausage + addnormal(pop[i,],mask,sigma=usage.sigma) # This is slow!
#covariates(mask)$usage = ausage
## Make usage surface a SpatialPixelsDataFrame
#spdf.usage = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$usage))
#inspdf.usage = spdf.usage[trapboundary,]

# set up stuff:
traplim.x = bbox(trapboundary)[1,]
traplim.y = bbox(trapboundary)[2,]
errpts = data.frame(x=c(mean(traplim.x),mean(traplim.x),mean(traplim.x),traplim.x[1]   ,traplim.x[2]   ,
                        traplim.x[1],traplim.x[1],traplim.x[2],traplim.x[2]),
                    y=c(mean(traplim.y),traplim.y[1]   ,traplim.y[2]   ,mean(traplim.y),mean(traplim.y),
                        traplim.y[1],traplim.y[2],traplim.y[1],traplim.y[2]))
x0 = bbox(boundary)[1,1]
xmax = bbox(boundary)[1,2]
xmid = errpts[1,1]
x = c(x0,(xmid-dmax),xmid,(xmid+dmax),xmax)
sigma = c(max.error.sigma,max.error.sigma,error.sigma,max.error.sigma,max.error.sigma)

# Create unifrom error density surface (V slow!) and add do mask
acerr = rep(0,dim(mask)[1])
system.time(for(i in 1:dim(pop)[1]) acerr = acerr + addnormal(pop[i,],mask,sigma=error.sigma)) # This is slow!
covariates(mask)$acerr = acerr
# Make usage surface a SpatialPixelsDataFrame
spdf.acerr = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$acerr))
inspdf.acerr = spdf.acerr[trapboundary,]
plot(trapboundary)
plot(inspdf.acerr,col=terrain.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="darkgray",add=TRUE)

# Create unifrom error density surface with max.error.sigma (V slow!) and add do mask
maxacerr = rep(0,dim(mask)[1])
system.time(for(i in 1:dim(pop)[1]) maxacerr = maxacerr + addnormal(pop[i,],mask,sigma=max.error.sigma)) # This is slow!
covariates(mask)$maxacerr = maxacerr
# Make usage surface a SpatialPixelsDataFrame
spdf.maxacerr = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$maxacerr))
inspdf.maxacerr = spdf.maxacerr[trapboundary,]
plot(trapboundary)
plot(inspdf.maxacerr,col=terrain.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="darkgray",add=TRUE)

# Create increasing error density surface, up to max.error.sigma (V slow!) and add to mask
acerrtrend = rep(0,dim(mask)[1])
d = sigmad = rep(0,dim(pop)[1])
for(i in 1:dim(pop)[1]) d[i] = sqrt(sum((pop[i,]-errpts[1,])^2))
# make sigma increase up to sill of max.error.sigma at dmax:
for(i in 1:dim(pop)[1]) sigmad[i] = min(max.error.sigma,error.sigma + d[i]/dmax*(max.error.sigma-error.sigma))
system.time(for(i in 1:dim(pop)[1]) acerrtrend <- acerrtrend + addnormal(pop[i,],mask,sigma=sigmad[i])) # This is slow!
covariates(mask)$acerrtrend = acerrtrend
# Make usage surface a SpatialPixelsDataFrame
spdf.acerrtrend = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$acerrtrend))
inspdf.acerrtrend = spdf.acerrtrend[trapboundary,]
plot(trapboundary)
plot(inspdf.acerrtrend,col=terrain.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="darkgray",add=TRUE)




# Now plot with all plots in one jpg:
jpeg(file="./paper/example-densities.jpg",h=12,w=12,units="cm",res=720)
#par(mar=c(0.25,0.25,1.5,0.25))
layout(matrix(c(1,2,0,3,4,5,6,7,8),3,3,byrow=TRUE), 
       heights=c(lcm(4),lcm(4),lcm(2.5)), widths=rep(lcm(4),3))
#layout.show(8)
# density and AC plots
par(mar=c(0.25,0.25,1.5,0.25))
# Plot the ac density surface
plot(trapboundary,main="(a)")
plot(inspdf.D,col=heat.colors(40),what="image",add=TRUE)
# Plot the realised acs
#par(mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(b)")
plot(inpop,pch=19,cex=0.25,add=TRUE)
# realised density plots
# Plot with small error
#par(mfrow=c(1,3),mar=c(0.25,0.25,1.5,0.25))
#plot(trapboundary,main="(d)")
plot(trapboundary,main="(c)") # relabelled
plot(inspdf.acerr,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="darkgray",add=TRUE)
#par(mar=c(0.25,0.25,1.5,0.25))
#plot(trapboundary,main="(e)")
plot(trapboundary,main="(d)") # relabelled
plot(inspdf.maxacerr,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="darkgray",add=TRUE)
# Plot with small-to-large error
#par(mar=c(0.25,0.25,1.5,0.25))
#plot(trapboundary,main="(f)")
plot(trapboundary,main="(e)") # relabelled
plot(inspdf.acerrtrend,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="darkgray",add=TRUE)
# sigma plots
plot(c(x0,xmax),rep(error.sigma,2),type="l",xaxt="n",ylim=c(0,max.error.sigma+1),xlab="",ylab="",main="(f)")
segments(c(x0+buffer,xmax-buffer),rep(0,2),c(x0+buffer,xmax-buffer),rep(max.error.sigma+1,2),lty=2)
plot(c(x0,xmax),rep(max.error.sigma,2),type="l",xaxt="n",xlab="",ylab="",ylim=c(0,max.error.sigma+1),main="(g)")
segments(c(x0+buffer,xmax-buffer),rep(0,2),c(x0+buffer,xmax-buffer),rep(max.error.sigma+1,2),lty=2)
plot(x,sigma,type="l",xaxt="n",xlab="",ylab="",ylim=c(0,max.error.sigma+1),main="(h)")
segments(c(x0+buffer,xmax-buffer),rep(0,2),c(x0+buffer,xmax-buffer),rep(max.error.sigma+1,2),lty=2)
dev.off()



# Do an SCR survey to illustrate prediction error size change
# ===========================================================
scrtraps = make.grid(nx=4,ny=4,spacex=10,detector="count",origin=c(35,35))
plot(boundary)
plot(scrtraps,add=TRUE)
plot(trapboundary,add=TRUE)
set.seed(1)
simch = sim.capthist(scrtraps,pop,detectpar=list(g0=1,sigma=10),noccasions=1)
summary(simch)
plot(simch,tracks=TRUE,border=0)
simfit = secr.fit(simch,mask=mask)
fxtot = fx.total(simfit,mask=mask)
#plotcovariate(fxtot,covariate="D.sum",what="image")
#plot(scrtraps,add=TRUE)
simfit$capthist[1,,] = c(rep(0,15),1) # make top right only detection for this animal


pdf(file="./paper/screrr.pdf",h=4,w=4)
par(mar=c(1,1,1,1))
plot(trapboundary)
fxi.contour(simfit,i=1,nx=200,add=TRUE,drawlabels=FALSE)
fxi.contour(simfit,i=9,nx=200,add=TRUE,drawlabels=FALSE) # this is an animal in the centre of the grid
plot(scrtraps,add=TRUE)
ch9 = simfit$capthist[9,1,]
detected9 = which(ch9>0)
points(scrtraps$x[detected9],scrtraps$y[detected9],pch=15)
ch1 = simfit$capthist[1,1,]
detected1 = which(ch1>0)
points(scrtraps$x[detected1],scrtraps$y[detected1],pch=17)
dev.off()

