# --------------------- New functions ----------------------------
addcircle=function(x0,y0,r,n=720,...){
  da=2*pi/n
  a=0
  xi=x0+r*sin(a)
  yi=y0+r*cos(a)
  for(i in 1:n){
    lines(c(xi,x0+r*sin(a+da)),c(yi,y0+r*cos(a+da)),...)
    a=a+da
    xi=x0+r*sin(a)
    yi=y0+r*cos(a)
  }
}
# ---------------------------------------------------------------


library(jpeg)
library(fields)
#monalisaRGB=readJPEG("/Users/dlb/Research/SCR book/bookV1/monalisa.jpg")
monalisaRGB=readJPEG("monalisa.jpg")
z=0.3*monalisaRGB[,,1] + 0.59*monalisaRGB[,,2] + 0.11*monalisaRGB[,,3]
monalisa=tz=t(z)
nj=dim(tz)[2]
for(j in 1:nj) monalisa[,j]=tz[,nj-j+1]
x=1:2835
y=1:4289
x=x-ceiling(mean(x))
y=y-ceiling(mean(y))
image(x,y,monalisa)
persp(x,y,monalisa)

save(monalisa,file="/Users/dlb/Research/SECR Book/SECRbookR/monalisa.rda")
load(file="/Users/dlb/Research/SECR Book/SECRbookR/monalisa.rda")

# just look at top half:
x=1:2835
y=1:4289
#ml=monalisa[,(4289-2000):4289]
#y=y[(4289-2000):4289]
ml=monalisa[,1989:3989]
y=y[1989:3989]
y=y-ceiling(mean(y))
x=x-ceiling(mean(x))
ml=ml/max(ml) # scale between 0 and 1
writeJPEG(ml,target="/Users/dlb/Research/SECR Book/SECRbookR/ml.jpg",quality=0.1)



ml=readJPEG("/Users/dlb/Research/SECR Book/SCR book/bookV1/ml.jpg")
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
addcircle(0,0,100,col="red")


greymonalisa=mltt
ml.x=xt
ml.y=yt
image(ml.x,ml.y,greymonalisa,col=gray(seq(0,1,length=100)),useRaster=TRUE,asp=1,
      xlab="",ylab="",xaxt="n",yaxt="n",bty="n")

# create mesh
nx=length(ml.x) # unique xs
ny=length(ml.y) # unique ys
mesh.x=rep(ml.x,nx) # mesh xs (faster changing index)
mesh.y=rep(ml.y,rep(nx,ny)) # mesh ys (slower changing index)
meshn=nx*ny;meshn

mesh.D=rep(NA,meshn)
for(i in 1:nx) for(j in 1:ny) mesh.D[(i-1)*ny+j]=greymonalisa[i,j]

# That is way too big: reduce by d^2 (with a bit of cropping)
d=25
nis=round(nx/d)
njs=round(ny/d)
nis;njs;nis*njs
mlsml=data.frame(x=rep(1:nis,njs),y=rep(1:njs,rep(nis,njs)),D=rep(NA,nis*njs))
for(j in 1:njs) for(i in 1:nis) {
  is = ((i-1)*d) : (i*d)
  js = ((j-1)*d) : (j*d)
  ind=(j-1)*nis + i
  mlsml$D[ind] = mean(greymonalisa[is,js])
}
mlmesh=read.mask(data=mlsml,spacing=1,columns="D")
plotcovariate(mlmesh,covariate="D",col=gray.colors(50))
plotcovariate(mlmesh,covariate="D",contour=FALSE,col=gray.colors(50))
covariates(mlmesh)$logD=log(1+covariates(mlmesh)$D)
plotcovariate(mlmesh,covariate="logD",contour=FALSE,col=gray.colors(50))
save(mlmesh,file="monalisamesh.RData")



library(fields)
library(secr)
source("/Users/dlb/git/SCR-Book/scrmlebook/R/scrplotting.r")
load("/Users/dlb/Research/SCR book/bookV1/monalisamesh.RData")

plotcovariate(mlmesh,covariate="D",contour=FALSE,col=gray.colors(50))

# "habitat" covariate
habitat = covariates(mlmesh)$D-mean(covariates(mlmesh)$D)
covariates(mlmesh)$habitat = habitat
habitat = (log((cD+1)))
plot(cD,habitat)
covariates(mlmesh)$habitat = habitat
plotcovariate(mlmesh,covariate="habitat",contour=FALSE,col=gray.colors(50))

mlDbase = exp(max(covariates(mlmesh)$D)-covariates(mlmesh)$D)*20 - 20
scaleup = 100
mlD = mlDbase * scaleup
pop=sim.popn(D=mlD, core=mlmesh, model2D="IHP", seed=12345)
N=dim(pop)[1];N # check simulated population size
# plot mesh with individuals' locations and detectors overlaid
plot(mlmesh, border=1,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.05)

# make a grid of detectors
dets=make.grid(nx=10,ny=10,spacex=8,spacey=6,originxy=c(22,13),detector="count")
plot(mlmesh, border=1,dots=FALSE, col="white", meshcol="gray")
plot(dets,add=TRUE)


## Generate capture histories
## first set number of occasions
lambda0=0.5;sigma=5
nt <- 1
capthist=sim.capthist(dets,popn=pop, detectfn="HHN",detectpar=list(lambda0=lambda0,sigma=sigma), noccasions=nt, nsessions=1,seed=12345)
summary(capthist)
n=dim(capthist)[1];n
plot(mlmesh, border=5,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.25)
plot(dets,add=TRUE)
plot(capthist, border=sigma, tracks=TRUE, varycol=FALSE,gridlines=FALSE,rad=3,add=TRUE)

# estimate:
smfit1=secr.fit(capthist,model=list(D~s(x,y,k=10)),mask=mlmesh)
smfit1.Dhat = predictDsurface(smfit1)
names(covariates(smfit1.Dhat))
plotcovariate(smfit1.Dhat,covariate="D.0",contour=FALSE,col=gray.colors(40))
plot(dets,add=TRUE)
region.N(smfit1)


smfit2=secr.fit(capthist,model=list(D~s(x,y,k=25)),mask=mlmesh)
smfit2.Dhat = predictDsurface(smfit2)
names(covariates(smfit2.Dhat))
plotcovariate(smfit2.Dhat,covariate="D.0",contour=FALSE,col=gray.colors(40))
plot(dets,add=TRUE)
region.N(smfit2)

AIC(smfit1,smfit2)
#save(smfit1,smfit2,file="smfits1-2.RData")

# try with habitat suitability
h=covariates(mlmesh)$D
smfith=secr.fit(capthist,model=list(D~D),mask=mlmesh)
smfith.Dhat = predictDsurface(smfith)
plotcovariate(smfith.Dhat,covariate="D.0",contour=FALSE,col=gray.colors(40))
plot(dets,add=TRUE)
region.N(smfith)

plotcovariate(mlmesh,covariate="D",contour=FALSE,col=gray.colors(40))
plotcovariate(smfith.Dhat,covariate="D.0",contour=FALSE,col=gray.colors(40))

plot(covariates(smfith.Dhat)$D,covariates(smfith.Dhat)$D.0)







library(mvtnorm)
sigma=50


#x=1:100;y=1:120;sigma=10;ml=matrix(rep(1,12000),nrow=100);ml=ml/sum(ml)

Sigma=diag(2)*sigma^2
ml.sm=ml
nx=length(x)
ny=length(y)
nsig=3.5
# calculate total density:
ys=xs=-round(nsig*sigma):round(nsig*sigma)
nxs=length(xs)
nys=length(ys)
X=matrix(c(rep(xs,nys),rep(ys,times=rep(nxs,nys))),ncol=2)
d=dmvnorm(X,mean=c(0,0),sigma=Sigma)
dsum=sum(d)
# now smear matrix ml.sm
for(i in 1:nx) for(j in 1:ny){
  xmean=x[i]
  ymean=y[j]
  ilo=max(1,i-round(nsig*sigma))
  ihi=min(nx,i+round(nsig*sigma))
  jlo=max(1,j-round(nsig*sigma))
  jhi=min(ny,j+round(nsig*sigma))
  xs=x[ilo:ihi]
  ys=y[jlo:jhi]
  nxs=length(xs)
  nys=length(ys)
  X=matrix(c(rep(xs,nys),rep(ys,times=rep(nxs,nys))),ncol=2)
  d=dmvnorm(X,mean=c(xmean,ymean),sigma=Sigma)
  d=(d/dsum)*(ml[i,j]) # normalize and weight
  for(k in 1:nys) {
    ml.sm[xs,ys[k]]=ml.sm[xs,ys[k]]+d[(k-1)*nxs+(1:nxs)]
  }
}

image(x,y,ml.sm,col=gray(seq(0,1,length=100)),useRaster=TRUE,asp=1,
      xlab="",ylab="",xaxt="n",yaxt="n",bty="n")

#image.plot(x,y,ml.sm,asp=1,bty="n")

