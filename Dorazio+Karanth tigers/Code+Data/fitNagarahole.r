require(secr)
require(scrmlebook)
require(pals)
require(plot3D)
require(sp)

tigerch = read.capthist("NagaraholeCHtimes.csv","Nagaraholetraps.csv",detector="count",noccasions=1,covnames=c("day","hour"))
cams = traps(tigerch)
quartz()
if (dim(tigerch)[2]>1) {
  par(mfrow=c(3,3))
  for(i in 1:dim(tigerch)[2]) plot(tigerch[[i]],border=0,tracks=TRUE)
} else {
  plot(tigerch,border=0,tracks=TRUE)
}
names(covariates(tigerch))
hist(covariates(tigerch)$day)
hist(covariates(tigerch)$hour,nclass=24)

# Get boundary polygon and create mask
sregion = readShapeSpatial("NHstatespace-utm.shp")
tigermask = make.mask(traps(tigerch),buffer=10000,spacing=500,type="trapbuffer",poly=sregion) # buffer and spacing consistent with Dorazio+Kranth (2017)
xlim=range(tigermask$x)
ylim=range(tigermask$y)

pdf(file="./keepfigure/NagaraholeMaskTraps.pdf",h=6,w=4)
par(mar=c(2,2,2,2))
plot(tigermask,dots=FALSE,border=0,xlim=xlim,ylim=ylim,ppoly=FALSE)
plotMaskEdge(tigermask,add=TRUE)
plot(cams, add=TRUE, detpar=list(col="black"))
dev.off()


fit0 = secr.fit(tigerch,mask=tigermask)
# Fit some density trend models
fitxy = secr.fit(tigerch,mask=tigermask,model=list(D~x+y),start=list(g0=detectpar(fit0)$g0,sigma=detectpar(fit0)$sigma))
fitx2 = secr.fit(tigerch,mask=tigermask,model=list(D~s(x,k=3)),start=list(g0=detectpar(fit0)$g0,sigma=detectpar(fit0)$sigma))
fity2 = secr.fit(tigerch,mask=tigermask,model=list(D~s(y,k=3)),start=list(g0=detectpar(fit0)$g0,sigma=detectpar(fit0)$sigma))
fitx2y2 = secr.fit(tigerch,mask=tigermask,model=list(D~s(x,k=3)+s(y,k=3)),start=list(g0=detectpar(fit0)$g0,sigma=detectpar(fit0)$sigma))
fitxy22 = secr.fit(tigerch,mask=tigermask,model=list(D~s(x,y,k=8)),start=list(g0=detectpar(fit0)$g0,sigma=detectpar(fit0)$sigma))
fitx = secr.fit(tigerch,mask=tigermask,model=list(D~x),start=list(g0=detectpar(fit0)$g0,sigma=detectpar(fit0)$sigma))
fity = secr.fit(tigerch,mask=tigermask,model=list(D~y),start=list(g0=detectpar(fit0)$g0,sigma=detectpar(fit0)$sigma))
fity3 = secr.fit(tigerch,mask=tigermask,model=list(D~s(y,k=4)),start=list(g0=detectpar(fit0)$g0,sigma=detectpar(fit0)$sigma))

# Put x and y into traps and then try model with  g0~y, sigma~y"
covariates(cams)$x = scale(cams$x)
covariates(cams)$y = scale(cams$y)
traps(tigerch) = cams
fity.g.s = secr.fit(tigerch,mask=tigermask,model=list(D~x+y, g0~y, sigma~y))

aics = AIC(fit0,fitx,fity,fitxy,fitx2,fity2,fitx2y2,fitxy22,fity3,fity.g.s)

save(fit0,fitx,fity,fitxy,fitx2,fity2,fitx2y2,fitxy22,fity3,fity.g.s,file="NagaraholeFits.RData")
     
model.average(fit0,fitx,fity,fitxy,fitx2,fity2,fitx2y2,fitxy22,fity3, realnames="D",criterion="AIC")
collate(fit0,fitx,fity,fitxy,fitx2,fity2,fitx2y2,fitxy22,fity3, realnames="D")

M=dim(tigermask)[1]
covdf = data.frame(D.0=rep(NA,M),SE.0=rep(NA,M),lcl.0=rep(NA,M),ucl.0=rep(NA,M),cv=rep(NA,M))
pred0D = predict(fit0)["D",]
covdf$cv = (rep(pred0D$SE.estimate,M)/rep(pred0D$estimate,M))*100
covdf[c("D.0","SE.0","lcl.0","ucl.0")] = 
  data.frame(D.0=rep(pred0D$estimate,M),SE.0=rep(pred0D$SE.estimate,M),lcl.0=rep(pred0D$lcl,M),ucl.0=rep(pred0D$ucl,M))*100^2
pred0 = tigermask
covariates(pred0) = covdf
predy2 = predictDsurface(fity2,mask=tigermask,se.D=TRUE,cl.D=TRUE)
covariates(predy2)$cv = (covariates(predy2)$SE.0/covariates(predy2)$D.0)*100
covariates(predy2)[c("D.0","SE.0","lcl.0","ucl.0")] = covariates(predy2)[c("D.0","SE.0","lcl.0","ucl.0")]*100^2
predy = predictDsurface(fity,mask=tigermask,se.D=TRUE,cl.D=TRUE)
covariates(predy)$cv = (covariates(predy)$SE.0/covariates(predy)$D.0)*100
covariates(predy)[c("D.0","SE.0","lcl.0","ucl.0")] = covariates(predy)[c("D.0","SE.0","lcl.0","ucl.0")]*100^2
predxy = predictDsurface(fitxy,mask=tigermask,se.D=TRUE,cl.D=TRUE)
covariates(predxy)$cv = (covariates(predxy)$SE.0/covariates(predxy)$D.0)*100
covariates(predxy)[c("D.0","SE.0","lcl.0","ucl.0")] = covariates(predxy)[c("D.0","SE.0","lcl.0","ucl.0")]*100^2

zlim=range(0,covariates(predxy)$D.0,covariates(predy)$D.0,covariates(predy2)$D.0,covariates(pred0)$D.0)
#quartz(h=8,w=6)
pdf(file="./keepfigure/NagaraholeDhats.pdf",h=8,w=6)
par(mfrow=c(2,2),mar=c(3,2,3,2))
scrmlebook:::plotcovariate(predy,covariate="D.0",contour=FALSE,col=parula(40),asp=1,zlim=zlim,main=paste("D~y (AICcwt=",round(aics[1,"AICcwt"]*100),"%)",sep=""),
                           bty="n",xaxt="n",yaxt="n")
scrmlebook:::plotcovariate(predy2,covariate="D.0",contour=FALSE,col=parula(40),asp=1,zlim=zlim,main=paste("D~s(y,3) (AICcwt=",round(aics[2,"AICcwt"]*100),"%)",sep=""),
                           bty="n",xaxt="n",yaxt="n")
scrmlebook:::plotcovariate(predxy,covariate="D.0",contour=FALSE,col=parula(40),asp=1,zlim=zlim,main=paste("D~x+y (AICcwt=",round(aics[3,"AICcwt"]*100),"%)",sep=""),
                           bty="n",xaxt="n",yaxt="n")
scrmlebook:::plotcovariate(pred0,covariate="D.0",contour=FALSE,col=parula(40),asp=1,zlim=zlim,main=paste("D~1 (AICcwt=",round(aics[4,"AICcwt"]*100),"%)",sep=""),
                           bty="n",xaxt="n",yaxt="n")
dev.off()

zlim=range(0,covariates(predxy)$SE.0,covariates(predy)$SE.0,covariates(predy2)$SE.0,covariates(pred0)$SE.0)
#quartz(h=8,w=6)
pdf(file="./keepfigure/NagaraholeSEs.pdf",h=8,w=6)
par(mfrow=c(2,2),mar=c(3,2,3,2))
scrmlebook:::plotcovariate(predy,covariate="SE.0",contour=FALSE,col=parula(40),asp=1,zlim=zlim,main=paste("D~y (AICcwt=",round(aics[1,"AICcwt"]*100),"%)",sep=""),
                           bty="n",xaxt="n",yaxt="n")
scrmlebook:::plotcovariate(predy2,covariate="SE.0",contour=FALSE,col=parula(40),asp=1,zlim=zlim,main=paste("D~s(y,3) (AICcwt=",round(aics[2,"AICcwt"]*100),"%)",sep=""),
                           bty="n",xaxt="n",yaxt="n")
scrmlebook:::plotcovariate(predxy,covariate="SE.0",contour=FALSE,col=parula(40),asp=1,zlim=zlim,main=paste("D~x+y (AICcwt=",round(aics[3,"AICcwt"]*100),"%)",sep=""),
                           bty="n",xaxt="n",yaxt="n")
scrmlebook:::plotcovariate(pred0,covariate="SE.0",contour=FALSE,col=parula(40),asp=1,zlim=zlim,main=paste("D~1 (AICcwt=",round(aics[4,"AICcwt"]*100),"%)",sep=""),
                           bty="n",xaxt="n",yaxt="n")
dev.off()

#quartz(h=8,w=6)
zlim=range(0,covariates(predxy)$cv,covariates(predy)$cv,covariates(predy2)$cv,covariates(pred0)$cv)
pdf(file="./keepfigure/NagaraholeCVs.pdf",h=8,w=6)
par(mfrow=c(2,2),mar=c(3,2,3,2))
scrmlebook:::plotcovariate(predy,covariate="cv",contour=FALSE,col=parula(40),asp=1,zlim=zlim,main=paste("D~y (AICcwt=",round(aics[1,"AICcwt"]*100),"%)",sep=""),
                           bty="n",xaxt="n",yaxt="n")
scrmlebook:::plotcovariate(predy2,covariate="cv",contour=FALSE,col=parula(40),asp=1,zlim=zlim,main=paste("D~s(y,3) (AICcwt=",round(aics[2,"AICcwt"]*100),"%)",sep=""),
                           bty="n",xaxt="n",yaxt="n")
scrmlebook:::plotcovariate(predxy,covariate="cv",contour=FALSE,col=parula(40),asp=1,zlim=zlim,main=paste("D~x+y (AICcwt=",round(aics[3,"AICcwt"]*100),"%)",sep=""),
                           bty="n",xaxt="n",yaxt="n")
scrmlebook:::plotcovariate(pred0,covariate="cv",contour=FALSE,col=parula(40),asp=1,zlim=zlim,main=paste("D~1 (AICcwt=",round(aics[4,"AICcwt"]*100),"%)",sep=""),
                           bty="n",xaxt="n",yaxt="n")
dev.off()

Nhat0 = region.N(fit0)
Nhaty2 = region.N(fity2)
Nhaty = region.N(fity)
Nhatxy = region.N(fitxy)

Nhat0
Nhaty2
Nhaty
Nhatxy

# predict in y dimension
ylim = range(tigermask$y)
xmean = mean(tigermask$x)
ymask = read.mask(data=data.frame(x=rep(xmean,200),y=seq(ylim[1],ylim[2],length=200)))
Dhat1 = covariates(predictDsurface(fity,mask=ymask,se.D=TRUE,cl.D=TRUE))*100^2
Dlim1 = range(Dhat1$lcl.0,Dhat1$ucl.0)
Dhat2 = covariates(predictDsurface(fity2,mask=ymask,se.D=TRUE,cl.D=TRUE))*100^2
Dlim2 = range(Dhat2$lcl.0,Dhat2$ucl.0)
Dhat3 = covariates(predictDsurface(fity3,mask=ymask,se.D=TRUE,cl.D=TRUE))*100^2
Dlim3 = range(Dhat3$lcl.0,Dhat3$ucl.0)
#Dlim = range(Dlim1,Dlim2,Dlim3)

#quartz(h=4,w=10)

pdf(file="./keepfigure/NagaraholeSmooths.pdf",h=4,w=10)
par(mfrow=c(1,3))
# df=2
plot(ymask$y,Dhat1$D.0,type="l",ylim=Dlim1,xlab="Northing",ylab="Density at mean Easting",lwd=1.5,main="D~y")
lines(ymask$y,Dhat1$lcl.0,lty=2)
lines(ymask$y,Dhat1$ucl.0,lty=2)
lines(ylim,rep(pred0D$estimate*100^2,2),col="gray")
lines(ymask$y,Dhat1$D.0,lwd=1.5)
# df=3
plot(ymask$y,Dhat2$D.0,type="l",ylim=Dlim2,xlab="Northing",ylab="Density at mean Easting",lwd=1.5,main="D~s(y,3)")
lines(ymask$y,Dhat2$lcl.0,lty=2)
lines(ymask$y,Dhat2$ucl.0,lty=2)
lines(ylim,rep(pred0D$estimate*100^2,2),col="gray")
lines(ymask$y,Dhat2$D.0,lwd=1.5)
# df=4
plot(ymask$y,Dhat3$D.0,type="l",ylim=Dlim3,xlab="Northing",ylab="Density at mean Easting",lwd=1.5,main="D~s(y,4)")
lines(ymask$y,Dhat3$lcl.0,lty=2)
lines(ymask$y,Dhat3$ucl.0,lty=2)
lines(ylim,rep(pred0D$estimate*100^2,2),col="gray")
lines(ymask$y,Dhat3$D.0,lwd=1.5)
dev.off()

# Spatial distribution of captures per trap, and best density surface
ncaps = apply(tigerch,3,sum) # capture frequency for each trap

pdf(file="./keepfigure/NagaraholeSurface.pdf",h=4,w=10)
par(mfrow=c(1,3))
par(mar=c(5, 4, 4, 2) + 0.1)
pcols = parula(12)
imf.y2=plotcovariate(predy,covariate="D.0",contour=FALSE,col=parula(44),asp=1,xlab="Easting",ylab="Northing",bty="n",zlim=zlim)
plot(sregion,col="gray",border="gray",add=TRUE)
points(cams$x,cams$y,cex=log(ncaps+1)/2+0.2,xlim=xlim,ylim=ylim,asp=1,col=pcols[ncaps+1],pch=19)
imf.y2=plotcovariate(predy,covariate="D.0",contour=FALSE,col=parula(44),asp=1,xlab="Easting",ylab="Northing",bty="n",zlim=zlim)
xlim=range(imf.y2$x)
ylim=range(imf.y2$y)
zlim = c(0,max(imf.y2$z[!is.na(imf.y2$z)]))
par(mar=c(3,2,3,2))
hist3D(imf.y2$x,imf.y2$y,imf.y2$z,col=parula(100),phi=45,theta=45,d=1.5,r=5, clim=zlim,
       xlim=xlim,ylim=ylim,zlim=zlim,xlab="Easting",ylab="Northing",zlab="Density",scale=FALSE,expand=3000)
#segments3D(cams$x,cams$y,rep(0,length(cams$x)),cams$x,cams$y,ncaps,col="gray",lty=2,lwd=1,phi=35,theta=45,d=1.5,r=5, 
#           xlim=xlim,ylim=ylim,xlab="Easting",ylab="Northing",zlab="Capture Frequency",scale=FALSE,expand=3000)
#points3D(cams$x,cams$y,z=ncaps,col=parula(11),pch=19,add=TRUE,cex=1,pch=19,phi=35,theta=45,d=1.5,r=5,scale=FALSE,expand=3000)
dev.off()


# Lat and lon of reserve (for GIS data, if I ever manage to get any!)
lat = c(12.1859,11.854158) # South
lon = c(76.0583,76.25061)  # East

