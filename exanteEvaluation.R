# load relevant packages
library(sp)
library(gstat)
library(rgdal)
library(spcosa)
library(sampling)

####################################################
#Model input
####################################################

#Read shape file of study area
area <- readOGR(
  dsn = ".",
  layer = "esri_level1_german_states_3035"
)
spplot(area[, "PROV1NAME"])

#Select a square grid from the area. This grid is used as the sampling frame
grid<-spsample(area,type="regular",cellsize=2500)
mySamFrame<-data.frame(grid)

#Add stratum ids to grid
strat<-over(x=grid,y=as(area,"SpatialPolygons"))
mySamFrame$stratumId<-strat

#Merge stratum Bremen with surrounding stratum
mySamFrame$stratumId[mySamFrame$stratumId==1]<-3
mySamFrame$stratumId<-mySamFrame$stratumId-1

#Compute number of strata
nstrata<-length(unique(mySamFrame$stratumId))

#Compute stratum sizes and weights
Nh<-tapply(mySamFrame$x1,INDEX=mySamFrame$stratumId,FUN=length)

#NB order of Nh must be same as in mySamFrame, see help of function strata of R package sampling
Nh<-Nh[unique(mySamFrame$stratumId)]

#Compute stratum weights
wh<-Nh/sum(Nh)

###########################################################
#Parameter settings
###########################################################

#1. Design parameters

#1.1 Specify type of space-time design
#SS: Static-Synchronous
#IS: Independent-Synchronous
#SA: Serially Alternating
#SP: Supplemented Panel
STDesign<-"SA"

#1.2 Specify type of spatial design
#srswor: (stratified) simple random sampling without replacement
#srswr: (stratified) simple random sampling with replacement
SDesign<-"srswr"

#1.3 Set total number of monitoring locations per sampling time
n<-30

#Compute stratum sizes for proportional allocation
nh<-round(n*wh)

#1.4 Set period for serially alternating design.
#A period of 2 means that at the third time the locations of the first time are revisited et cetera
period<-2

#1.5 For supplemented panel: set proportion of locations that will be revisited
#NB Take care that the number of revisited locations per stratum, nhSS, must be >= 2 in order to estimate sampling covariances with argument experimental = FALSE
revisit<-0.5
nhIS<-floor(nh*revisit)
nhSS<-nh-nhIS

#1.6 Set sampling interval and end of monitoring period
interval<-1
end<-3

#Compute sampling times
stimes<-seq(from=0,to=end,by=interval)

#Compute number of sampling times
ntimes<-length(stimes)

#2. Parameters of the space--time variogram (see Table 1 Heuvelink and Griffith, 2010)
c0.s<-1.9E-5 #nugget variance spatial variogram
c.s<-13.3E-5 #sill variance of spatial variogram
a.s<-35000 #spatial range of spatial variogram in m.
c0.t<-0 #nugget variance of time variogram
c.t<-0.3E-5 #sill variance of time variogram
a.t<-6 #temporal range of time variogram in months 
c0.st<-0.1E-5 #nugget variance of space--time variogram
c.st<-0.7E-5 #sill variance of space--time variogram
a.st.s<-6000 #spatial range space--time variogram in m.
alpha<-1E4
a.st.t<-a.st.s/alpha #temporal range space--time variogram (computed as spatial range divided by alpha)

#2.1 Define sum-metric model; NB function vgmST cannot yet be used, as krigeST does not support yet simulation
vgm.1 = vgm(c0.s, "Exp", 1E12, 
  anis = c(0,90,0,1E-15,1E-15)) # spatial nugget
vgm.2 = vgm(c.s, "Bes", a.s*1E8,
  anis = c(0,90,0,1E-4,1E-4), add.to = vgm.1) # spatial psill
vgm.3 = vgm(c0.t, "Exp", 1E12,
  anis = c(0,0,0,1,1E-15), add.to = vgm.2) # temporal nugget
vgm.4 = vgm(c.t, "Exp", a.t*1E6, anis = c(0,0,0,1,1E-6),
  add.to = vgm.3) # temporal psill
vgm.5 = vgm(c0.st, "Exp", 1E-12, add.to = vgm.4) # space-time nugget
vgm.summetric = vgm(c.st, "Exp", a.st.s,
  anis = c(0,0,0,1,1/alpha), add.to = vgm.5) # space-time sill


#2.2 Hereafter, it is convenient to have the pure spatial, pure temporal variogram and the geometric anisotropic space-time variogram
vgm.space = vgm(psill=c.s, "Bes", range=a.s, nugget=c0.s)
vgm.time = vgm(psill=c.t, "Exp", range=a.t, nugget=c0.t)
vgm.spacetime = vgm(psill=c.st, "Exp", range=a.st.s,nugget=c0.st, anis = c(0,0,0,1,1/alpha))


#3. Parameters of Monte Carlo approximation
nsim<-1000 #number of simulations, 1000
nsam<-50 #number of samples per simulation
npairs<-1E6 #number of pairs of points to approximate mean semivariance

#set random seed (for reproduction of results)
set.seed(1415)

############################################################
#Now define functions
############################################################

#1. Function for drawing space-time samples
stSample<-function(nstrata,nh,mySamFrame,stimes,stdesign,sdesign,nhSS,period,nsamples) {
  ntimes<-length(stimes)
  nhIS<-nh-nhSS
  n<-sum(nh)
  if (stdesign=="SS") {
    xsall<-ysall<-tsall<-hsall<-NULL
    for (sam in 1:nsamples) {
      mySample<-strata(mySamFrame,stratanames="stratumId",size=nh,method=sdesign)
      xsall<-c(xsall,rep(mySamFrame$x1[mySample$ID_unit],ntimes))
      ysall<-c(ysall,rep(mySamFrame$x2[mySample$ID_unit],ntimes))
      tsall<-c(tsall,rep(stimes,each=n))
      hsall<-c(hsall,rep(mySample$stratumId,ntimes))
    }
    xytall<-data.frame(xsall,ysall,tsall,hsall)
  }
  if (stdesign=="IS") {
    xsall<-ysall<-tsall<-hsall<-NULL    
    for (sam in 1:nsamples) {
      mySample<-strata(mySamFrame,stratanames="stratumId",size=nh*ntimes,method=sdesign)
      xsall<-c(xsall,mySamFrame$x1[mySample$ID_unit])
      ysall<-c(ysall,mySamFrame$x2[mySample$ID_unit])
      ts<-NULL
      for (i in unique(mySamFrame$stratumId)) {
        ts<-rep(stimes,each=nh[names(nh)==i])
        tsall<-c(tsall,ts)
      }    
      hsall<-c(hsall,mySample$stratumId)
    }
    xytall<-data.frame(xsall,ysall,tsall,hsall)
  }
  if (stdesign=="SA") {
    xsall<-ysall<-tsall<-hsall<-NULL
    for (sam in 1:nsamples) {
    xs<-ys<-hs<-NULL    
      for (i in 1:period) {
        mySample<-strata(mySamFrame,stratanames="stratumId",size=nh,method=sdesign)
        xs<-c(xs,mySamFrame$x1[mySample$ID_unit])
        ys<-c(ys,mySamFrame$x2[mySample$ID_unit])
        hs<-c(hs,mySample$stratumId)
      }
      xsrep<-ysrep<-hsrep<-numeric(length=n*ntimes)
      xsrep[]<-rep(xs)
      ysrep[]<-rep(ys)
      hsrep[]<-rep(hs)
      xsall<-c(xsall,xsrep)
      ysall<-c(ysall,ysrep)
      tsall<-c(tsall,rep(stimes,each=n))
      hsall<-c(hsall,hsrep)
    }
    xytall<-data.frame(xsall,ysall,tsall,hsall)
  }
  if (stdesign=="SP") {
    xsall<-ysall<-tsall<-hsall<-panelall<-NULL    
    for (sam in 1:nsam) {
      mySample<-strata(mySamFrame,stratanames="stratumId",size=nhSS,method=sdesign)
      xsr<-mySamFrame$x1[mySample$ID_unit]
      ysr<-mySamFrame$x2[mySample$ID_unit]
      xsSS<-rep(xsr,ntimes)
      ysSS<-rep(ysr,ntimes)
      tsSS<-rep(stimes,each=sum(nhSS))
      hsSS<-rep(mySample$stratumId,ntimes)
      mySample<-strata(mySamFrame,stratanames="stratumId",size=nhIS*ntimes,method=sdesign)
      xsIS<-mySamFrame$x1[mySample$ID_unit]
      ysIS<-mySamFrame$x2[mySample$ID_unit]
      tIS<-tsIS<-NULL
      for (i in unique(mySamFrame$stratumId)) {
        tIS<-rep(stimes,each=nhIS[names(nhIS)==i])
        tsIS<-c(tsIS,tIS)
      }
      hsIS<-mySample$stratumId
      xsall<-c(xsall,xsSS,xsIS)
      ysall<-c(ysall,ysSS,ysIS)
      tsall<-c(tsall,tsSS,tsIS)
      hsall<-c(hsall,hsSS,hsIS)
      panelall<-c(panelall,c(rep(1,sum(nhSS)*ntimes),rep(2,sum(nhIS)*ntimes)))
    }
    xytall<-data.frame(xsall,ysall,tsall,hsall,panelall)
  }
  xytall
}


#2. Function for estimating spatial means at the sampling times.
EstimateMean<-function(dat,stdesign,wh){
  if (stdesign=="SP") {
    stratumMeans<-tapply(dat$z,INDEX=list(dat$sam,dat$t,dat$h,dat$panel),FUN=mean)
    wstratMeans<-array(dim=dim(stratumMeans))
    for (i in 1:nstrata){
      wstratMeans[,,i,]<-stratumMeans[,,i,]*wh[names(wh)==i]
    }
    spatialMeans<-apply(wstratMeans,MARGIN=c(1,2,4),FUN=sum)
  }else{
    stratumMeans<-tapply(dat$z,INDEX=list(dat$sam,dat$t,dat$h),FUN=mean)
    wstratMeans<-array(dim=dim(stratumMeans))
    for (i in 1:nstrata) {
      wstratMeans[,,i]<-stratumMeans[,,i]*wh[names(wh)==i]
    }
    spatialMeans<-apply(wstratMeans,MARGIN=c(1,2),FUN=sum)
  }
  spatialMeans
}


#3. Function for estimating sampling covariace matrix of the estimated spatial means
EstimateCp<-function(dat,stdesign,nh,nhSS,experimental,wh,period){
  ZeroCp<-function(Cp,stdesign,period){
    if (stdesign=="SS") {
      Cp0<-Cp
    }
    if (stdesign=="IS") {
      Cp0<-diag(nrow(Cp))
      diag(Cp0)<-diag(Cp)
    }
    if (stdesign=="SA") {
      delta<-abs(row(Cp)-col(Cp))
      Cp0<-Cp
      Cp0[delta%%period>0]<-0
    }
    if (stdesign=="SP") {
      ntimes<-nrow(Cp)/2
      Cp0<-Cp
      delta<-abs(row(Cp)-col(Cp))
      Cp0[delta>0 & row(Cp)>ntimes & col(Cp)>ntimes]<-0
    }
    Cp0
  }
  nsam<-length(unique(dat$sam))
  nstrata<-length(unique(dat$h))
  ntimes<-length(unique(dat$t))
  unique(dat$h)
  nhIS<-nh-nhSS
  Cph<-wCph<-CpISh<-wCpISh<-CpSSh<-wCpSSh<-corh<-array(dim=c(ntimes,ntimes,nstrata))
  if (experimental==TRUE) {
    spatialMeans<-EstimateMean(dat=dat,stdesign=stdesign,wh=wh)
    if (stdesign=="SP") {
      meansPanel1<-spatialMeans[,,1]
      meansPanel2<-spatialMeans[,,2]
      CpPanel1<-var(meansPanel1)
      CpPanel2<-var(meansPanel2)
      Cp<-matrix(data=0,ncol=2*ntimes,nrow=2*ntimes)
      Cp[1:ntimes,1:ntimes]<-CpPanel1
      Cp[(ntimes+1):(2*ntimes),(ntimes+1):(2*ntimes)]<-CpPanel2
    }else{
    Cp<-var(spatialMeans)
      if (stdesign!="SS") {
      Cp<-ZeroCp(Cp,stdesign,period)
      }
    }
  }else{ 
    if (stdesign=="SP") {
      for (j in 1:nstrata) {
        arraysimh<-array(dim=c(nh[names(nh)==j]*nsam,ntimes,nstrata))          
        array2simh<-array(dim=c(nhSS[names(nhSS)==j]*nsam,ntimes,nstrata))        
        for (k in 1:ntimes) {
          arraysimh[,k,j]<-dat$z[(dat$h==j & dat$t==stimes[k])]
          array2simh[,k,j]<-dat$z[(dat$h==j & dat$t==stimes[k] & dat$panel == 1)]
        }
        CpISh[,,j]<-var(arraysimh[,,j])/nhIS[names(nhIS)==j]
        wCpISh[,,j]<-CpISh[,,j]*wh[names(wh)==j]^2
        corh[,,j]<-cor(array2simh[,,j])
        Sh<-diag(sqrt(var(arraysimh[,,j])))
        S2h<-outer(Sh,Sh)
        CpSSh[,,j]<-corh[,,j]*S2h/nhSS[names(nhSS)==j]
        wCpSSh[,,j]<-CpSSh[,,j]*wh[names(wh)==j]^2
      }
      CpIS<-apply(wCpISh,MARGIN=c(1,2),FUN=sum)
      CpSS<-apply(wCpSSh,MARGIN=c(1,2),FUN=sum)
      Cp<-matrix(data=0,ncol=2*ntimes,nrow=2*ntimes)
      Cp[1:ntimes,1:ntimes]<-CpSS
      Cp[(ntimes+1):(2*ntimes),(ntimes+1):(2*ntimes)]<-CpIS      
    }else{
      for (j in 1:nstrata) {
        arraysimh<-array(dim=c(nh[names(nh)==j]*nsam,ntimes,nstrata))          
        for (k in 1:ntimes) {
          arraysimh[,k,j]<-dat$z[(dat$h==j & dat$t==stimes[k])]
        }
        Cph[,,j]<-var(arraysimh[,,j])/nh[names(nh)==j]
        wCph[,,j]<-Cph[,,j]*wh[names(wh)==j]^2
      }
      Cp<-apply(wCph,MARGIN=c(1,2),FUN=sum)      
    }
      if (stdesign!="SS") {
      Cp<-ZeroCp(Cp,stdesign,period)
      }    
  }
  Cp
}

#############################################################
#Computations
#############################################################

#1. Regularization of space--time variogram and computation of model covariance matrix of spatial means

#1.1 Approximate the mean semivariance on point support within 2D block (time-lag 0), gbar(A,A)
sample1<-as.data.frame(spsample(area,npairs,type="random"))
sample2<-as.data.frame(spsample(area,npairs,type="random"))
dx<-sqrt((sample1$x-sample2$x)^2)
dy<-sqrt((sample1$y-sample2$y)^2)
dxy<-sqrt(dx^2+dy^2)
gs <-variogramLine(vgm.space,dist_vector=dxy)
gst<-variogramLine(vgm.spacetime, dir=c(0,1,0), dist_vector=dxy)
gbarAA<-mean(gs$gamma+gst$gamma)
rm(sample1,sample2)

#1.2 Compute mean semivariance on point support between two 2D-blocks separated by time-lag u
gbarAAu<-numeric(length=ntimes)
for (i in 1:ntimes) {
  dt<-rep(stimes[i],times=npairs)
  gt<-variogramLine(vgm.time,dist_vector=dt)
  dxyt<-sqrt(dxy^2+(dt*alpha)^2)
  gst<-variogramLine(vgm.spacetime,dist_vector=dxyt,dir=c(1,0,0))
  gbarAAu[i]<-mean(gs$gamma+gt$gamma+gst$gamma)
  }

#1.3 Compute semivariances on block-support
g_A<-gbarAAu-gbarAA

#1.4 Compute covariances on 2D block support 
#Compute the sill of the space-time variogram on 2D block support as the mean semivariance between two spatial 2D blocks separated by an infinitely large time lag
sillst<-c0.t+c.t+c0.st+c.st+mean(gs$gamma)-gbarAA #see Eq. II.41, J&H, p.78
Cxiij<- sillst-g_A

#1.5 Fill the matrix Cxi
dum<-matrix(data=0,nrow=ntimes,ncol=ntimes)
tlag<-abs(row(dum)-col(dum))
Cxi<-matrix(data=0,nrow=ntimes,ncol=ntimes)
for (i in 1:ntimes) {
  Cxi[tlag==i-1]<-Cxiij[i]
}
save(Cxi,file="Cxi.RData")

#2. Selection of space--time samples and geostatistical simulation

#2.1 Select nsam (stratified) random samples 
stsamples<-stSample(nstrata=nstrata,nh=nh,mySamFrame=mySamFrame,stimes=stimes,stdesign=STDesign,sdesign=SDesign,nhSS=nhSS,period=period,nsamples=nsam)
coordinates(stsamples) <- ~xsall+ysall+tsall

#2.2 Check whether there are points with same coordinates and jitter if there are
js<-zerodist(stsamples, zero = 0.0, unique.ID = FALSE)
if (nrow(js)>0) {
  stsamples<-data.frame(stsamples)
  stsamples[js[,1],1]<-jitter(stsamples[js[,1],1],amount=1E-4)
  coordinates(stsamples) <- ~xsall+ysall+tsall  
}

#2.3 Gaussian simulation by means of simple kriging
simulation  <- krige(
  formula = dummy ~ 1,
#  data=,
  locations = NULL,
  newdata = stsamples,
  model = vgm.summetric,
  nmax = 100,
  nsim = nsim,
  beta = 0, 
  dummy = TRUE
)
simdf <- as.data.frame(simulation)


#3. Estimate determinant of covariance matrix of regression coefficients and variance of trend

#3.1 Compute design matrix
D<-matrix(nrow=ntimes,ncol=2)
D[,1]<-1
D[,2]<-stimes

#3.1 Compute design matrix for supplemented panel
Dsup<-matrix(nrow=2*ntimes,ncol=2)
Dsup[,1]<-1
Dsup[,2]<-c(stimes,stimes)

#3.2 Estimate model covariance matrix of regression coefficients
invCxi<-solve(Cxi)
Cxibeta<-solve(t(D)%*%invCxi%*%D)

#3.3 Start loop over simulations to estimate conditional sampling covariances matrix of estimated regression coefficients
Cpb<-array(data=0,dim=c(2,2,nsim))
samplenr<-as.numeric(rep(1:nsam,each=n*ntimes))
for (i in 1:nsim) {

#Make a dataframe containing an identifier for the space--time sample, the space--time coordinates, the simulated values, the stratum, and for the SP design, the panel
  if (STDesign=="SP") {
    dat<-data.frame(samplenr,simdf$xsall,simdf$ysall,simdf$tsall,simdf[[i+3]],data.frame(stsamples)$hsall,data.frame(stsamples)$panelall)
    names(dat)<-c("sam","x","y","t","z","h","panel")
  }else{
    dat<-data.frame(samplenr,simdf$xsall,simdf$ysall,simdf$tsall,simdf[[i+3]],data.frame(stsamples)$hsall)
    names(dat)<-c("sam","x","y","t","z","h")
  }

#Estimate sampling variance-covariance matrix of estimated means
  Cp<-EstimateCp(dat=dat,stdesign=STDesign,nh=nh,nhSS=nhSS,experimental=FALSE,wh=wh,period=period)

#Add xi-covariance matrix to sampling covariance matrix of estimated spatial means
  if (STDesign=="SP") {
  Xsup<-matrix(data=0,nrow=2*ntimes,ncol=ntimes)
  for (j in 1:ntimes) {
    Xsup[j,j]<-1
    Xsup[j+ntimes,j]<-1
  }
  XCxiX<-Xsup%*%Cxi%*%t(Xsup)
  Cxip<-XCxiX+Cp
  invCxip<-solve(Cxip)
  DCDinv<-solve(t(Dsup)%*%invCxip%*%Dsup)  
  
#Estimate spatial means and regression coefficients
  spatialMeans<-EstimateMean(dat=dat,stdesign=STDesign,wh=wh)
  b0GLS<-b1GLS<-numeric(length=nsam)
  for (samnr in 1:nsam) {
    spM<-as.numeric(spatialMeans[samnr,,])
    DCY<-t(Dsup)%*%invCxip%*%spM
    beta<-DCDinv%*%DCY
    b0GLS[samnr]<-beta[1]
    b1GLS[samnr]<-beta[2]
    }
  }else{
  Cxip<-Cxi+Cp
  invCxip<-solve(Cxip)
  DCDinv<-solve(t(D)%*%invCxip%*%D)  
  
#Estimate spatial means and regression coefficients
  spatialMeans<-EstimateMean(dat=dat,stdesign=STDesign,wh=wh)
  b0GLS<-b1GLS<-numeric(length=nsam)  
  for (samnr in 1:nsam) {
    DCY<-t(D)%*%invCxip%*%as.numeric(spatialMeans[samnr,])
    beta<-DCDinv%*%DCY
    b0GLS[samnr]<-beta[1]
    b1GLS[samnr]<-beta[2]
    }
  }
#Estimate conditional sampling variances of regression coefficients
  Cpb[1,1,i]<-var(b0GLS)
  Cpb[2,2,i]<-var(b1GLS)
  Cpb[1,2,i]<-Cpb[2,1,i]<-cov(b0GLS,b1GLS)
}

#########################################################
#Results
#########################################################

#Compute average of conditional sampling covariance matrices (averaged over simulations)
ExiCpb<-apply(Cpb,MARGIN=c(1,2),FUN=mean) 

#Add sampling- and model covariance matrices
Ctot<-ExiCpb+Cxibeta

#Compute determinants 
print(detCxip<-Ctot[1,1]*Ctot[2,2]-Ctot[1,2]^2)
print(detCxi<-Cxibeta[1,1]*Cxibeta[2,2]-Cxibeta[1,2]^2)
print(detCp<-ExiCpb[1,1]*ExiCpb[2,2]-ExiCpb[1,2]^2)

#Compute variances
print(Vxipb1<-Ctot[2,2])
print(Vxib1<-Cxibeta[2,2])
print(ExiVpb1<-ExiCpb[2,2])

#Make data.frame of results and write it to csv output file
print(result<-data.frame(STDesign,ntimes,n,detCxip,detCxi,detCp,Vxipb1,Vxib1,ExiVpb1))
#write.table(result,file="ExanteGermany_(4x30).csv",row.names=FALSE, sep=",")
write.table(result,file="ExanteGermany_(4x30).csv",row.names=FALSE, sep=",",append=TRUE,col.names=FALSE)

#Plot the results
results<-read.csv(file="ExanteGermany_(4x30).csv",header=TRUE)
names(results)<-c("Design","ntimes","n","detCxip","detCxi","detCp","Varxip","Varxi","Varp")
results$Design<-as.character(results$Design)
results$Design[results$Design=="SP"]<-"SuP"
results$Design<-as.factor(results$Design)
pdf(file = "ExanteGermany_(4x30).pdf", width = 7, height = 7)
print(
  ggplot() +
  geom_point(data=results,
     mapping = aes(
         x = detCxip*1E12,
         y = Varxip*1E7,
         shape=Design
     ),
     size = 3
 ) +
scale_x_continuous(name = "Determinant * 1E12") +
scale_y_continuous(name = "Variance trend *1E7"))
dev.off()