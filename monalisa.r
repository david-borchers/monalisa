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
monalisaRGB=readJPEG("/Users/dlb/Research/SECR Book/monalisa.jpg")
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

image(x,y,ml,col=gray(seq(0,1,length=100)),useRaster=TRUE,asp=1,
      xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
addcircle(0,0,100,col="red")

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

