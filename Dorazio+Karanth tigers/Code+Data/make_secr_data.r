trapdat = read.csv(file="trap_2015_coord.csv",header=TRUE)
chdat = read.csv(file="Capt_date_time_noHrep.csv",header=TRUE)

# Convert date to days since start, and time to seconds since midnight:
dt = dmy(chdat$Date)
day = dt - min(dt)
datime = dmy_hms(paste(chdat$Date,chdat$Time))
totmins = as.numeric(datime)/60 # convert to minutes
mins = totmins%%(24*60) # modulo hours
hours = mins/60 # convert to decimal hours

pdf("NagraholeCapTimes.pdf",h=4,w=8)
hist(hours,nclass=24,main="",xlab="Hour")
dev.off()

# Create data frame suitable for import to secr:
n = length(mins)
sess = as.numeric(day)%/%6+1 # divide into 8 6-day sessions
sess = rep(1,n) # make a single session
tigerchdf = data.frame(Session=sess, ID=chdat$ANIMAL_ID, Occasion=rep(1,n), Detector=chdat$LOC_ID, day=day, hour=hours)
chnames = names(tigerchdf); chnames[1]=paste("#",chnames[1],sep="")
trnames = names(trapdat); trnames[1]=paste("#",trnames[1],sep="")
write.table(tigerchdf,file="NagaraholeCHtimes.csv",row.names=FALSE,col.names=chnames,sep=",",quote=FALSE)
write.table(trapdat,file="Nagaraholetraps.csv",row.names=FALSE,col.names=trnames,sep=",",quote=FALSE)
