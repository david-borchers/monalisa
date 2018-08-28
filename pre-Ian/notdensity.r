library(jpeg)
library(fields)
library(scrmlebook)

#setwd("/Users/dlb/Research/SCR book/bookV1/")

ml=readJPEG("ml.jpg")
mlt=t(ml)
nx=dim(ml)[1]
ny=dim(ml)[2]
x=1:nx
y=1:ny
yt=x
xt=y
mltt=mlt[,order(yt,decreasing = TRUE)]
image(xt,yt,mltt,col=gray(seq(0,1,length=100)),useRaster=TRUE,asp=1,
      xlab="",ylab="",xaxt="n",yaxt="n",bty="n")

# get a mesh
load("monalisamesh.RData")

mlim = plotcovariate(mlmesh,covariate="D",contour=FALSE,col=gray.colors(50),asp=1,bty="n")

# "habitat" covariate
nmesh = length(mlmesh[,1]);nmesh
habitat = 10-log(10 + covariates(mlmesh)$D-mean(covariates(mlmesh)$D)) + rnorm(nmesh,0,0.02)
covariates(mlmesh)$habitat = habitat
plot(covariates(mlmesh)$D,habitat)
covariates(mlmesh)$habitat = habitat
mlveg = plotcovariate(mlmesh,covariate="habitat",contour=FALSE,col=tim.colors(50),xlim=c(30,90))
mlveg = image(mlveg)
#keep=mlpop$x<90 & mlpop$x>30
#mlcov = data.frame(x=mlmesh$x[keep],y=mlmesh$y[keep],veg=veg$z[keep])
save(mlveg,file="mlveg.RData")


dx = dy = 5
# make a grid of detectors
dets1=make.grid(nx=7,ny=7,spacex=dx,spacey=dy,originxy=c(25,20),detector="count")
plotcovariate(mlmesh,covariate="D",contour=FALSE,col=gray.colors(50),asp=1)
plot(dets1,add=TRUE)

dets2=make.grid(nx=7,ny=7,spacex=dx,spacey=dy,originxy=c(42.5,28.5),detector="count")
plotcovariate(mlmesh,covariate="D",contour=FALSE,col=gray.colors(50),asp=1)
plot(dets2,add=TRUE)

dets3=make.grid(nx=7,ny=7,spacex=dx,spacey=dy,originxy=c(60,38),detector="count")
plotcovariate(mlmesh,covariate="D",contour=FALSE,col=gray.colors(50),asp=1)
plot(dets3,add=TRUE)

pdf("arrays.pdf",h=4,w=4)
plotcovariate(mlmesh,covariate="D",contour=FALSE,col=gray.colors(50),asp=1,bty="n",key=FALSE,main="")
plot(dets1,add=TRUE,detpar=list(col="red"))
plot(dets2,add=TRUE,detpar=list(col="orange"))
plot(dets3,add=TRUE,detpar=list(col="green"))
dev.off()



#mlDbase = exp(max(covariates(mlmesh)$D)-covariates(mlmesh)$D)*20 - 20
#mlD = mlDbase * scaleup
covariates(mlmesh)$mlD = max(covariates(mlmesh)$D) - covariates(mlmesh)$D
scaleup = 6000
covariates(mlmesh)$mlD = covariates(mlmesh)$mlD * scaleup
mlD = covariates(mlmesh)$mlD

covariates(mlmesh)$mlD = max(covariates(mlmesh)$D) - covariates(mlmesh)$D
scaleup = 3000
covariates(mlmesh)$mlD = covariates(mlmesh)$mlD * scaleup
mlD = covariates(mlmesh)$mlD
mlpop=sim.popn(D=mlD, core=mlmesh, model2D="IHP", seed=1234565)
N=dim(mlpop)[1];N # check simulated population size
keep=mlpop$x<90 & mlpop$x>30
mlpop=data.frame(x=mlpop$x[keep],y=mlpop$y[keep])
plot(mlpop,asp=1,pch=19,cex=0.33)
save(mlpop,file="mlpop.RData")

covariates(mlmesh)$mlD = max(covariates(mlmesh)$D) - covariates(mlmesh)$D
scaleup = 2000
covariates(mlmesh)$mlD = covariates(mlmesh)$mlD * scaleup
mlD = covariates(mlmesh)$mlD
smlpop=sim.popn(D=mlD, core=mlmesh, model2D="IHP", seed=1234565)
N=dim(smlpop)[1];N # check simulated population size
keep=smlpop$x<90 & smlpop$x>30
smlpop=data.frame(x=smlpop$x[keep],y=smlpop$y[keep])
plot(smlpop,pch=19,asp=1,cex=0.5)
save(smlpop,file="smlpop.RData")

covariates(mlmesh)$mlD = max(covariates(mlmesh)$D) - covariates(mlmesh)$D
scaleup = 1000
covariates(mlmesh)$mlD = covariates(mlmesh)$mlD * scaleup
mlD = covariates(mlmesh)$mlD
vsmlpop=sim.popn(D=mlD, core=mlmesh, model2D="IHP", seed=54321)
N=dim(vsmlpop)[1];N # check simulated population size
keep=vsmlpop$x<90 & vsmlpop$x>30
vsmlpop=data.frame(x=vsmlpop$x[keep],y=vsmlpop$y[keep])
plot(vsmlpop,pch=19,asp=1,cex=0.5)
save(vsmlpop,file="vsmlpop.RData")


pop=sim.popn(D=mlD, core=mlmesh, model2D="IHP", seed=54321)
N=dim(pop)[1];N # check simulated population size

pdf("monalisaDpts.pdf",h=4,w=10)
par(mfrow=c(1,2))
mlDim = plotcovariate(mlmesh,covariate="D",contour=FALSE,col=gray.colors(50),asp=1,bty="n",key=FALSE,main="Density surface")
plot(pop$x,pop$y,pch=19,cex=0.2,asp=1,bty="n",main="Activity centres",xlab="x",ylab="y")
dev.off()

# plot mesh with individuals' locations and detectors overlaid
#plot(mlmesh, border=1,dots=FALSE, col="white", meshcol="gray")
#points(pop$x,pop$y,pch=19,cex=0.05)

## Generate capture histories
## first set number of occasions
g0=0.5;sigma=3
nt <- 1
capthist1=sim.capthist(dets1,popn=pop, detectfn="HN",detectpar=list(g0=g0,sigma=sigma), noccasions=nt, nsessions=1,seed=12345)
summary(capthist1)
n1=dim(capthist1)[1];n1
plot(mlmesh, border=5,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.25)
plot(dets1,add=TRUE)
plot(capthist1, border=sigma, tracks=TRUE, varycol=FALSE,gridlines=FALSE,rad=3,add=TRUE)

capthist2=sim.capthist(dets2,popn=pop, detectfn="HN",detectpar=list(g0=g0,sigma=sigma), noccasions=nt, nsessions=1,seed=12345)
summary(capthist2)
n2=dim(capthist2)[1];n2
plot(mlmesh, border=5,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.25)
plot(dets2,add=TRUE)
plot(capthist2, border=sigma, tracks=TRUE, varycol=FALSE,gridlines=FALSE,rad=3,add=TRUE)

capthist3=sim.capthist(dets3,popn=pop, detectfn="HN",detectpar=list(g0=g0,sigma=sigma), noccasions=nt, nsessions=1,seed=12345)
summary(capthist3)
n3=dim(capthist3)[1];n3
plot(mlmesh, border=5,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.25)
plot(dets3,add=TRUE)
plot(capthist3, border=sigma, tracks=TRUE, varycol=FALSE,gridlines=FALSE,rad=3,add=TRUE)

Nhat=rep(NA,3)
n=c(n1,n2,n3)
# estimate:
cfit1 = secr.fit(capthist1,mask=mlmesh)
fxi1 = fxi.secr(cfit1,i=1:n1)
N1 = region.N(cfit1)
Nhat[1] = round(N1["R.N","estimate"])

cfit2=secr.fit(capthist2,mask=mlmesh)
fxi2 = fxi.secr(cfit2,i=1:n2)
N2 = region.N(cfit2)
Nhat[2] = round(N2["R.N","estimate"])

cfit3=secr.fit(capthist3,mask=mlmesh)
fxi3 = fxi.secr(cfit3,i=1:n3)
N3 = region.N(cfit3)
Nhat[3] = round(N3["R.N","estimate"])


pd1 = pdot(mlmesh,dets1,detectfn="HN",detectpar=detectpar(cfit1),noccasions=nt)
pd2 = pdot(mlmesh,dets2,detectfn="HN",detectpar=detectpar(cfit2),noccasions=nt)
pd3 = pdot(mlmesh,dets3,detectfn="HN",detectpar=detectpar(cfit3),noccasions=nt)
pd = list(pd1/sum(pd1),pd2/sum(pd2),pd3/sum(pd3))
pdmiss = list((1-pd1)/sum(1-pd1),(1-pd2)/sum(1-pd2),(1-pd3)/sum(1-pd3))
covariates(mlmesh)$pd1 = pd[[1]]
covariates(mlmesh)$pd2 = pd[[2]]
covariates(mlmesh)$pd3 = pd[[3]]
covariates(mlmesh)$pd1miss = pdmiss[[1]]
covariates(mlmesh)$pd2miss = pdmiss[[2]]
covariates(mlmesh)$pd3miss = pdmiss[[3]]

plotcovariate(mlmesh,covariate="pd1miss",contour=FALSE,col=tim.colors(50))
plot(dets1,add=TRUE,detpar=list(cex=0.5,col="white"))




fxi. = matrix(rep(0,npts*3),nrow=3)
fxii = list(fxi1,fxi2,fxi3)
for(case in 1:3) {
  for(i in 1:n[case]) fxi.[case,] = fxi.[case,] + fxii[[case]][[i]]
  fxi.[case,] = (n[case]* fxi.[case,] + (Nhat[case]-n[case])*pdmiss[[case]])/Nhat[case]
}

covariates(mlmesh)$fxi1. = fxi.[1,] 
covariates(mlmesh)$fxi2. = fxi.[2,] 
covariates(mlmesh)$fxi3. = fxi.[3,] 
  
im1 = plotcovariate(mlmesh,covariate="fxi1.",contour=TRUE,col=gray.colors(50))
plot(dets1,add=TRUE,detpar=list(cex=0.5,col="white"))

im2 = plotcovariate(mlmesh,covariate="fxi2.",contour=TRUE,col=gray.colors(50))
plot(dets2,add=TRUE,detpar=list(cex=0.5,col="white"))

im3 = plotcovariate(mlmesh,covariate="fxi3.",contour=TRUE,col=gray.colors(50))
plot(dets3,add=TRUE,detpar=list(cex=0.5,col="white"))


xlim=matrix(rep(NA,2*3),nrow=3)
ylim=matrix(rep(NA,2*3),nrow=3)
dets = list(dets1,dets2,dets3)
for(case in 1:3) {
  xlim[case,] = range(dets[[case]]$x)
  ylim[case,] = range(dets[[case]]$y)
}
xlim[,1] = xlim[,1]-4*sigma
xlim[,2] = xlim[,2]+4*sigma
ylim[,1] = ylim[,1]-4*sigma
ylim[,2] = ylim[,2]+4*sigma

contour(mlim,bty="n",nlevels=3,xlim=xlim[1,],ylim=ylim[1,])
contour(im1,bty="n",xlim=xlim[1,],ylim=ylim[1,],add=TRUE)
#plot(dets1,add=TRUE,detpar=list(cex=0.5,col="red"))


pdf("globalests.pdf",h=8,w=12)
par(mfrow=c(2,3))
#image(mlim,bty="n",col=gray.colors(20))
contour(im1,bty="n")
#image(mlim,bty="n",col=gray.colors(20))
contour(im2,bty="n")
#image(mlim,bty="n",col=gray.colors(20))
contour(im3,bty="n")
contour(mlim,bty="n",nlevels=7)
contour(mlim,bty="n",nlevels=7)
contour(mlim,bty="n",nlevels=7)
dev.off()

pdf("localests.pdf",h=4,w=12)
par(mfrow=c(1,3))
image(mlim,bty="n",xlim=xlim[1,],ylim=ylim[1,],col=gray.colors(20))
contour(im1,bty="n",xlim=xlim[1,],ylim=ylim[1,],add=TRUE)
image(mlim,bty="n",xlim=xlim[2,],ylim=ylim[2,],col=gray.colors(20))
contour(im2,bty="n",xlim=xlim[2,],ylim=ylim[2,],add=TRUE)
image(mlim,bty="n",xlim=xlim[3,],ylim=ylim[3,],col=gray.colors(20))
contour(im3,bty="n",xlim=xlim[3,],ylim=ylim[3,],add=TRUE)
dev.off()

plot(pop$x,pop$y,pch=19,cex=0.1,bty="n",xlim=xlim[1,],ylim=ylim[1,])
contour(im1,bty="n",xlim=xlim[1,],ylim=ylim[1,],add=TRUE,col="gray")

plot(pop$x,pop$y,pch=19,cex=0.1,bty="n",xlim=xlim[2,],ylim=ylim[2,])
contour(im2,bty="n",xlim=xlim[2,],ylim=ylim[2,],add=TRUE,col="gray")

plot(pop$x,pop$y,pch=19,cex=0.1,bty="n",xlim=xlim[3,],ylim=ylim[3,])
contour(im3,bty="n",xlim=xlim[3,],ylim=ylim[3,],add=TRUE,col="gray")

pdf("overlap.pdf",h=4,w=8)
par(mfrow=c(1,2))
contour(im1,xlim=xlim[2,],ylim=ylim[2,],col="red")
contour(im2,xlim=xlim[2,],ylim=ylim[2,],add=TRUE,col="orange")
contour(im3,xlim=xlim[2,],ylim=ylim[2,],add=TRUE,col="blue")
contour(mlim,xlim=xlim[2,],ylim=ylim[2,],col=gray.colors(20))
dev.off()


contour(mlim,asp=1,bty="n")
contour(im2,asp=1,bty="n")
contour(mlim,asp=1,bty="n")
contour(im3,asp=1,bty="n")
