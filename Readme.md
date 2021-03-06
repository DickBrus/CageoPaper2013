exanteEvaluation.R: an R script for ex ante evaluation of space-time sampling designs for estimating the temporal trend of spatial means
==========================================================================================

Dick J. Brus and Jaap J. de Gruijter
------------------------------------

### Description

This ```R``` script can be used for ex ante evaluation of space-time sampling designs for trend monitoring. An a priori geostatistical model describing the variation of the variable of interest in space and time is required. The parameters of interest are the two regression coefficients (intercept and slope) of a time-series model of the spatial means with added linear trend. The output is the predicted variance of the estimated linear trend (slope coefficient), and the predicted determinant of the variance-covariance matrix of the two coefficients. The method is described in the paper D.J. Brus and J.J. de Gruijter (2013). Effects of spatial pattern persistence on the performance of sampling designs for regional trend monitoring analyzed by simulation
of space–time fields. Computers and Geosciences 61, 75-83. (http://dx.doi.org/10.1016/j.cageo.2013.09.001)

### Dependencies
This ```R``` script uses several ```R``` packages:


```r
library(sp)
library(gstat)
library(rgdal)
```

```
## rgdal: version: 0.8-11, (SVN revision 479M) Geospatial Data Abstraction
## Library extensions to R successfully loaded Loaded GDAL runtime: GDAL
## 1.9.2, released 2012/10/08 Path to GDAL shared files: d:/R/lib/rgdal/gdal
## GDAL does not use iconv for recoding strings. Loaded PROJ.4 runtime: Rel.
## 4.7.1, 23 September 2009, [PJ_VERSION: 470] Path to PROJ.4 shared files:
## d:/R/lib/rgdal/proj
```

```r
library(sampling)
```

```
## Loading required package: MASS Loading required package: lpSolve
```


### Model input
The ```R``` script needs as input a SpatialGrid. When a SpatialPolygonsDataFrame is read, first a square grid is sampled using function ```spsample```. Finally the SpatialGrid is coerced to a data.frame, that serves as the sampling frame. The spatial sampling design implemented is (stratified) simple random sampling with or without replacement. A column must be added to the data.frame containing the stratum identifiers of all grid nodes. The number of strata, and the number of grid nodes (sizes) per stratum are computed. These sizes are needed in design-based estimation of the spatial means (totals) per sampling time. 

```r
# Read shape file of study area
area <- readOGR(dsn = ".", layer = "esri_level1_german_states_3035")
```

```
## OGR data source with driver: ESRI Shapefile 
## Source: ".", layer: "esri_level1_german_states_3035"
## with 4 features and 13 fields
## Feature type: wkbPolygon with 2 dimensions
```

```r

# Select a square grid from the area. This grid is used as the sampling
# frame
grid <- spsample(area, type = "regular", cellsize = 2500)
mySamFrame <- data.frame(grid)

# Add stratum ids to grid
strat <- over(x = grid, y = as(area, "SpatialPolygons"))
mySamFrame$stratumId <- strat

# Merge stratum Bremen with surrounding stratum
mySamFrame$stratumId[mySamFrame$stratumId == 1] <- 3
mySamFrame$stratumId <- mySamFrame$stratumId - 1

# Compute number of strata
nstrata <- length(unique(mySamFrame$stratumId))

# Compute stratum sizes
Nh <- tapply(mySamFrame$x1, INDEX = mySamFrame$stratumId, FUN = length)

# NB order of Nh must be same as in mySamFrame, see help of function strata
# of R package sampling
Nh <- Nh[unique(mySamFrame$stratumId)]

# Compute stratum weights
wh <- Nh/sum(Nh)
```


### Parameter settings
The next step is to specify:
* sampling design parameters
* parameters of the space-time model
* parameters of the Monte Carlo approximation

#### Sampling design parameters
The following sampling design parameters must be specified
* type of space--time design $STDesign$
* type of spatial sampling design $SDesign$
* number of sampling locations per time $n$
* sampling interval $interval$
* end of the monitoring period $end$
* for serially alternating design the period $period$
* for supplemented panel design the proportion of locations that will be revisited $revisit$

The spatial sample sizes per stratum $nh$ are computed by multiplying the total number of sampling locations per time $n$ by relative sizes of the strata $wh$ (proportional allocation).

Take care that for the supplemented panel design the number of revisited locations per stratum, $nhSS$, equals 2 or larger. This is needed for estimation of the sampling covariances with argument ```experimental == FALSE```.


```r
# Specify type of space-time design SS: Static-Synchronous IS:
# Independent-Synchronous SA: Serially Alternating SP: Supplemented Panel
STDesign <- "SP"

# Specify type of spatial design srswor: (stratified) simple random sampling
# without replacement srswr: (stratified) simple random sampling with
# replacement
SDesign <- "srswr"

# Set total number of sampling locations per sampling time
n <- 30

# Compute stratum sample sizes for proportional allocation
nh <- round(n * wh)
print(sum(nh))
```

```
## [1] 30
```

```r

# Set period for serially alternating design. A period of 2 means that at
# the third time the locations of the first time are revisited et cetera
period <- 2

# For supplemented panel: set proportion of locations that will be revisited
# NB Take care that the number of revisited locations per stratum, nhSS is
# at least 2
revisit <- 0.5
nhIS <- floor(nh * revisit)
nhSS <- nh - nhIS
print(nhSS)
```

```
## 3 1 2 
## 3 5 7
```

```r

# Set sampling interval and end of monitoring period
interval <- 1
end <- 3

# Compute sampling times
stimes <- seq(from = 0, to = end, by = interval)

# Compute number of sampling times
ntimes <- length(stimes)
```


#### Parameters of space-time model
The space-time model is a sum-metric model, i.e. the sum of a pure spatial model, pure temporal model, and a geometric anisotropic space-time model. The sum-metric model is specified using function ```vgm```, with arguments ```add.to``` and ```anis```. This is not easy, we hope that in the near future this can be done with the function ```vgmST```. We could not make use of thus fucntion yet because function ```krigeST``` does not support yet geostatistical simulation (only prediction).


```r
# Set parameters of the space--time variogram (see Table 1 Heuvelink and
# Griffith, 2010)
c0.s <- 1.9e-05  #nugget variance spatial variogram
c.s <- 0.000133  #sill variance of spatial variogram
a.s <- 35000  #spatial range of spatial variogram in m.
c0.t <- 0  #nugget variance of time variogram
c.t <- 3e-06  #sill variance of time variogram
a.t <- 6  #temporal range of time variogram in months 
c0.st <- 1e-06  #nugget variance of space--time variogram
c.st <- 7e-06  #sill variance of space--time variogram
a.st.s <- 6000  #spatial range space--time variogram in m.
alpha <- 10000
a.st.t <- a.st.s/alpha  #temporal range space--time variogram (computed as spatial range divided by alpha)

# Define sum-metric model; NB function vgmST cannot yet be used, as krigeST
# does not support yet simulation
vgm.1 = vgm(c0.s, "Exp", 1e+12, anis = c(0, 90, 0, 1e-15, 1e-15))  # spatial nugget
vgm.2 = vgm(c.s, "Bes", a.s * 1e+08, anis = c(0, 90, 0, 1e-04, 1e-04), add.to = vgm.1)  # spatial psill
vgm.3 = vgm(c0.t, "Exp", 1e+12, anis = c(0, 0, 0, 1, 1e-15), add.to = vgm.2)  # temporal nugget
vgm.4 = vgm(c.t, "Exp", a.t * 1e+06, anis = c(0, 0, 0, 1, 1e-06), add.to = vgm.3)  # temporal psill
vgm.5 = vgm(c0.st, "Exp", 1e-12, add.to = vgm.4)  # space-time nugget
vgm.summetric = vgm(c.st, "Exp", a.st.s, anis = c(0, 0, 0, 1, 1/alpha), add.to = vgm.5)  # space-time sill

# Hereafter it is convenient to have the pure spatial, pure temporal
# variogram and the geometric anisotropic space-time variogram
vgm.space = vgm(psill = c.s, "Bes", range = a.s, nugget = c0.s)
vgm.time = vgm(psill = c.t, "Exp", range = a.t, nugget = c0.t)
vgm.spacetime = vgm(psill = c.st, "Exp", range = a.st.s, nugget = c0.st, anis = c(0, 
    0, 0, 1, 1/alpha))
```

#### Monte Carlo parameters
Three parameters must be set
* number of geostatsical simulations $nsim$
* number of space-time samples per simulation $nsam$
* number of pairs of points for approximation of mean semivariance between two 2D blocks (with geometry equal to that of study area) separated in time


```r
nsim <- 1000  #number of simulations
nsam <- 50  #number of samples per simulation
npairs <- 1e+06  #number of pairs of points to approximate mean semivariance

# set random seed (for reproduction of results)
set.seed(1415)
```


### Functions
Three functions are used:

* ```stSample```
* ```EstimateMean```
* ```EstimateCp}```

#### ```stSample```
The function ```stSample``` draws repeatedly stratified random samples. This function makes use of the function ```strata``` of the ```R``` package ```sampling```. This function can also be used for selecting simple random samples by setting all values in column stratumId of data.frame mySamFrame to the same value. The output of the function  ```stSample``` is a data.frame containing the coordinates in space-time, the stratum, and for space-time design SP the panel.


```r
stSample <- function(nstrata, nh, mySamFrame, stimes, stdesign, sdesign, nhSS, 
    period, nsamples) {
    ntimes <- length(stimes)
    nhIS <- nh - nhSS
    n <- sum(nh)
    if (stdesign == "SS") {
        xsall <- ysall <- tsall <- hsall <- NULL
        for (sam in 1:nsamples) {
            mySample <- strata(mySamFrame, stratanames = "stratumId", size = nh, 
                method = sdesign)
            xsall <- c(xsall, rep(mySamFrame$x1[mySample$ID_unit], ntimes))
            ysall <- c(ysall, rep(mySamFrame$x2[mySample$ID_unit], ntimes))
            tsall <- c(tsall, rep(stimes, each = n))
            hsall <- c(hsall, rep(mySample$stratumId, ntimes))
        }
        xytall <- data.frame(xsall, ysall, tsall, hsall)
    }
    if (stdesign == "IS") {
        xsall <- ysall <- tsall <- hsall <- NULL
        for (sam in 1:nsamples) {
            mySample <- strata(mySamFrame, stratanames = "stratumId", size = nh * 
                ntimes, method = sdesign)
            xsall <- c(xsall, mySamFrame$x1[mySample$ID_unit])
            ysall <- c(ysall, mySamFrame$x2[mySample$ID_unit])
            ts <- NULL
            for (i in unique(mySamFrame$stratumId)) {
                ts <- rep(stimes, each = nh[names(nh) == i])
                tsall <- c(tsall, ts)
            }
            hsall <- c(hsall, mySample$stratumId)
        }
        xytall <- data.frame(xsall, ysall, tsall, hsall)
    }
    if (stdesign == "SA") {
        xsall <- ysall <- tsall <- hsall <- NULL
        for (sam in 1:nsamples) {
            xs <- ys <- hs <- NULL
            for (i in 1:period) {
                mySample <- strata(mySamFrame, stratanames = "stratumId", size = nh, 
                  method = sdesign)
                xs <- c(xs, mySamFrame$x1[mySample$ID_unit])
                ys <- c(ys, mySamFrame$x2[mySample$ID_unit])
                hs <- c(hs, mySample$stratumId)
            }
            xsrep <- ysrep <- hsrep <- numeric(length = n * ntimes)
            xsrep[] <- rep(xs)
            ysrep[] <- rep(ys)
            hsrep[] <- rep(hs)
            xsall <- c(xsall, xsrep)
            ysall <- c(ysall, ysrep)
            tsall <- c(tsall, rep(stimes, each = n))
            hsall <- c(hsall, hsrep)
        }
        xytall <- data.frame(xsall, ysall, tsall, hsall)
    }
    if (stdesign == "SP") {
        xsall <- ysall <- tsall <- hsall <- panelall <- NULL
        for (sam in 1:nsam) {
            mySample <- strata(mySamFrame, stratanames = "stratumId", size = nhSS, 
                method = sdesign)
            xsr <- mySamFrame$x1[mySample$ID_unit]
            ysr <- mySamFrame$x2[mySample$ID_unit]
            xsSS <- rep(xsr, ntimes)
            ysSS <- rep(ysr, ntimes)
            tsSS <- rep(stimes, each = sum(nhSS))
            hsSS <- rep(mySample$stratumId, ntimes)
            mySample <- strata(mySamFrame, stratanames = "stratumId", size = nhIS * 
                ntimes, method = sdesign)
            xsIS <- mySamFrame$x1[mySample$ID_unit]
            ysIS <- mySamFrame$x2[mySample$ID_unit]
            tIS <- tsIS <- NULL
            for (i in unique(mySamFrame$stratumId)) {
                tIS <- rep(stimes, each = nhIS[names(nhIS) == i])
                tsIS <- c(tsIS, tIS)
            }
            hsIS <- mySample$stratumId
            xsall <- c(xsall, xsSS, xsIS)
            ysall <- c(ysall, ysSS, ysIS)
            tsall <- c(tsall, tsSS, tsIS)
            hsall <- c(hsall, hsSS, hsIS)
            panelall <- c(panelall, c(rep(1, sum(nhSS) * ntimes), rep(2, sum(nhIS) * 
                ntimes)))
        }
        xytall <- data.frame(xsall, ysall, tsall, hsall, panelall)
    }
    xytall
}
```


#### Function ```EstimateMean```
The function ```EstimateMean``` estimates the spatial means at the sampling times by design-based inference.


```r
EstimateMean <- function(dat, stdesign, wh) {
    if (stdesign == "SP") {
        stratumMeans <- tapply(dat$z, INDEX = list(dat$sam, dat$t, dat$h, dat$panel), 
            FUN = mean)
        wstratMeans <- array(dim = dim(stratumMeans))
        for (i in 1:nstrata) {
            wstratMeans[, , i, ] <- stratumMeans[, , i, ] * wh[names(wh) == 
                i]
        }
        spatialMeans <- apply(wstratMeans, MARGIN = c(1, 2, 4), FUN = sum)
    } else {
        stratumMeans <- tapply(dat$z, INDEX = list(dat$sam, dat$t, dat$h), FUN = mean)
        wstratMeans <- array(dim = dim(stratumMeans))
        for (i in 1:nstrata) {
            wstratMeans[, , i] <- stratumMeans[, , i] * wh[names(wh) == i]
        }
        spatialMeans <- apply(wstratMeans, MARGIN = c(1, 2), FUN = sum)
    }
    spatialMeans
}
```

#### Function ```EstimateCp```
This function estimates the sampling covariace matrix of the estimated spatial means at the sampling times. There are two options for estimating these covariances. When the argument ```experimental``` is set to TRUE, the $ntimes$ spatial means are estimated for each of the $nsam$ space-time samples. The variances and covariances of these estimated spatial means are used as an estimate of the sampling (co)variances. When this argument is set to FALSE, the sampling covariance matrix is estimated by the Horvitz-Thompson estimator of the sampling variance for (stratified) simple random sampling. The implemented estimator does not account for a finite population correction, so for (stratified) simple random sampling from finite populations with large sampling fractions (within strata), this variance estimator must be adapted. The *spatial* variances (population variances) required for this estimator are estimated from *all* $nsam$ space-time samples, which is of course a much better estimate than from a single sample. For stratified random sampling with 1 point per stratum, the sampling covariance matrix must be estimated with ```experimental==TRUE```. Our experience is that quite a few samples are required to obtain reliable estimates of the sampling covariance matrix. We recommend to set $nsam$ to 1000 or larger.


```r
EstimateCp <- function(dat, stdesign, nh, nhSS, experimental, wh, period) {
    ZeroCp <- function(Cp, stdesign, period) {
        if (stdesign == "SS") {
            Cp0 <- Cp
        }
        if (stdesign == "IS") {
            Cp0 <- diag(nrow(Cp))
            diag(Cp0) <- diag(Cp)
        }
        if (stdesign == "SA") {
            delta <- abs(row(Cp) - col(Cp))
            Cp0 <- Cp
            Cp0[delta%%period > 0] <- 0
        }
        if (stdesign == "SP") {
            ntimes <- nrow(Cp)/2
            Cp0 <- Cp
            delta <- abs(row(Cp) - col(Cp))
            Cp0[delta > 0 & row(Cp) > ntimes & col(Cp) > ntimes] <- 0
        }
        Cp0
    }
    nsam <- length(unique(dat$sam))
    nstrata <- length(unique(dat$h))
    ntimes <- length(unique(dat$t))
    unique(dat$h)
    nhIS <- nh - nhSS
    Cph <- wCph <- CpISh <- wCpISh <- CpSSh <- wCpSSh <- corh <- array(dim = c(ntimes, 
        ntimes, nstrata))
    if (experimental == TRUE) {
        spatialMeans <- EstimateMean(dat = dat, stdesign = stdesign, wh = wh)
        if (stdesign == "SP") {
            meansPanel1 <- spatialMeans[, , 1]
            meansPanel2 <- spatialMeans[, , 2]
            CpPanel1 <- var(meansPanel1)
            CpPanel2 <- var(meansPanel2)
            Cp <- matrix(data = 0, ncol = 2 * ntimes, nrow = 2 * ntimes)
            Cp[1:ntimes, 1:ntimes] <- CpPanel1
            Cp[(ntimes + 1):(2 * ntimes), (ntimes + 1):(2 * ntimes)] <- CpPanel2
        } else {
            Cp <- var(spatialMeans)
            if (stdesign != "SS") {
                Cp <- ZeroCp(Cp, stdesign, period)
            }
        }
    } else {
        if (stdesign == "SP") {
            for (j in 1:nstrata) {
                arraysimh <- array(dim = c(nh[names(nh) == j] * nsam, ntimes, 
                  nstrata))
                array2simh <- array(dim = c(nhSS[names(nhSS) == j] * nsam, ntimes, 
                  nstrata))
                for (k in 1:ntimes) {
                  arraysimh[, k, j] <- dat$z[(dat$h == j & dat$t == stimes[k])]
                  array2simh[, k, j] <- dat$z[(dat$h == j & dat$t == stimes[k] & 
                    dat$panel == 1)]
                }
                CpISh[, , j] <- var(arraysimh[, , j])/nhIS[names(nhIS) == j]
                wCpISh[, , j] <- CpISh[, , j] * wh[names(wh) == j]^2
                corh[, , j] <- cor(array2simh[, , j])
                Sh <- diag(sqrt(var(arraysimh[, , j])))
                S2h <- outer(Sh, Sh)
                CpSSh[, , j] <- corh[, , j] * S2h/nhSS[names(nhSS) == j]
                wCpSSh[, , j] <- CpSSh[, , j] * wh[names(wh) == j]^2
            }
            CpIS <- apply(wCpISh, MARGIN = c(1, 2), FUN = sum)
            CpSS <- apply(wCpSSh, MARGIN = c(1, 2), FUN = sum)
            Cp <- matrix(data = 0, ncol = 2 * ntimes, nrow = 2 * ntimes)
            Cp[1:ntimes, 1:ntimes] <- CpSS
            Cp[(ntimes + 1):(2 * ntimes), (ntimes + 1):(2 * ntimes)] <- CpIS
        } else {
            for (j in 1:nstrata) {
                arraysimh <- array(dim = c(nh[names(nh) == j] * nsam, ntimes, 
                  nstrata))
                for (k in 1:ntimes) {
                  arraysimh[, k, j] <- dat$z[(dat$h == j & dat$t == stimes[k])]
                }
                Cph[, , j] <- var(arraysimh[, , j])/nh[names(nh) == j]
                wCph[, , j] <- Cph[, , j] * wh[names(wh) == j]^2
            }
            Cp <- apply(wCph, MARGIN = c(1, 2), FUN = sum)
        }
        if (stdesign != "SS") {
            Cp <- ZeroCp(Cp, stdesign, period)
        }
    }
    Cp
}
```

### Computations

#### Regularization of space--time variogram
The matrix with model covariances of the spatial means is estimated by regularization of the space--time variogram on point-support.

First, the semivariogram on 2D-block support is computed by $\gamma_A(u) = \bar{\gamma}(A,A_u) - \bar{\gamma}(A,A)$, with $u$ the time-lag,  For this the mean semivariance *on point support* within 2D block (time-lag 0), $\bar{\gamma}(A,A)$ is computed first by selecting two simple random samples of size $npairs$. These two samples are used to construct $npairs$ of pairs of points. These are then used to compute $npairs$ semivariances from the spatial variogram and from the space--time variogram. After adding these two semivariance components, the mean semivariance is computed.

```r
sample1 <- as.data.frame(spsample(area, npairs, type = "random"))
sample2 <- as.data.frame(spsample(area, npairs, type = "random"))
dx <- sqrt((sample1$x - sample2$x)^2)
dy <- sqrt((sample1$y - sample2$y)^2)
dxy <- sqrt(dx^2 + dy^2)
gs <- variogramLine(vgm.space, dist_vector = dxy)
gst <- variogramLine(vgm.spacetime, dir = c(0, 1, 0), dist_vector = dxy)
gbarAA <- mean(gs$gamma + gst$gamma)
rm(sample1, sample2)
print(gbarAA)
```

```
## [1] 0.0001471
```

Next the mean semivariance *on point support* between two 2D-blocks separated by time-lag $u$, $\bar{\gamma}(A,A_u)$ is computed.

```r
gbarAAu <- numeric(length = ntimes)
for (i in 1:ntimes) {
    dt <- rep(stimes[i], times = npairs)
    gt <- variogramLine(vgm.time, dist_vector = dt)
    gst <- c0.st + c.st * (1 - exp(-1 * (sqrt((dxy/a.st.s)^2 + (dt/a.st.t)^2))))
    gbarAAu[i] <- mean(gs$gamma + gt$gamma + gst)
}
print(gbarAAu)
```

```
## [1] 0.0001471 0.0001476 0.0001480 0.0001483
```

Now the semivariance *on block-support* is computed, see equation above, and from this the covariance on block-support. The sill of the space-time variogram on 2D block support is computed as the mean semivariance between two spatial 2D blocks separated by an infinitely large time lag $u$. Finally the model covariance matrix is filled.

```r
g_A <- gbarAAu - gbarAA

# compute sill
sillst <- c0.t + c.t + c0.st + c.st + mean(gs$gamma) - gbarAA  #see Eq. II.41, Journel and Huijbregts, p.78
Cxiij <- sillst - g_A

# Fill the matrix Cxi
dum <- matrix(data = 0, nrow = ntimes, ncol = ntimes)
tlag <- abs(row(dum) - col(dum))
Cxi <- matrix(data = 0, nrow = ntimes, ncol = ntimes)
for (i in 1:ntimes) {
    Cxi[tlag == i - 1] <- Cxiij[i]
}
print(Cxi)
```

```
##           [,1]      [,2]      [,3]      [,4]
## [1,] 3.014e-06 2.547e-06 2.152e-06 1.820e-06
## [2,] 2.547e-06 3.014e-06 2.547e-06 2.152e-06
## [3,] 2.152e-06 2.547e-06 3.014e-06 2.547e-06
## [4,] 1.820e-06 2.152e-06 2.547e-06 3.014e-06
```


#### Selection of space--time samples and geostatistical simulation
Now $nsam$ space-time samples are selected with the specified design. At the sampling events $z$-values are then simulated using function ````krige```` of ```R``` package ```gstat```. First $nsam$ (stratified) simple random samples are selected

```r
stsamples <- stSample(nstrata = nstrata, nh = nh, mySamFrame = mySamFrame, stimes = stimes, 
    stdesign = STDesign, sdesign = SDesign, nhSS = nhSS, period = period, nsamples = nsam)
coordinates(stsamples) <- ~xsall + ysall + tsall
```

By chance there can be sampling events with exactly the same space--time coordinates. This leads to problems with the geostatistical simulation. To avoid this problem duplicates are jittered, i.e. a randomly chosen very small value is added to the first spatial coordinate.

```r
# Check whether there are points with same coordinates
js <- zerodist(stsamples, zero = 0, unique.ID = FALSE)

# Jitter first coordinate of duplicate
if (nrow(js) > 0) {
    stsamples <- data.frame(stsamples)
    stsamples[js[, 1], 1] <- jitter(stsamples[js[, 1], 1], amount = 1e-04)
    coordinates(stsamples) <- ~xsall + ysall + tsall
}
```

Now simulate values at the selectEd points in space-time

```r
# Gaussian simulation by means of simple kriging
simulation <- krige(formula = dummy ~ 1, locations = NULL, newdata = stsamples, 
    model = vgm.summetric, nmax = 100, nsim = nsim, beta = 0, dummy = TRUE)
```

```
## [using unconditional Gaussian simulation]
```

```r
simdf <- as.data.frame(simulation)
```

#### Estimate determinant of covariance matrix and variance of trend
Finally for each simulation, the covariance matrix of the estimated model parameters (intercept and slope /trend) is estimated by GLS. In GLS the covariance matrix is the sum of the model covariance matrix and the sampling covariance matrix.


```r
# Compute design matrix
D <- matrix(nrow = ntimes, ncol = 2)
D[, 1] <- 1
D[, 2] <- stimes

# Compute design matrix for supplemented panel
Dsup <- matrix(nrow = 2 * ntimes, ncol = 2)
Dsup[, 1] <- 1
Dsup[, 2] <- c(stimes, stimes)

# Estimate model covariance matrix of regression coefficients
invCxi <- solve(Cxi)
Cxibeta <- solve(t(D) %*% invCxi %*% D)

# Estimate model covariance matrix of regression coefficients
invCxi <- solve(Cxi)
Cxibeta <- solve(t(D) %*% invCxi %*% D)

# Start loop over simulations to estimate matrix with conditional sampling
# variances of estimated regression coefficients
Cpb <- array(data = 0, dim = c(2, 2, nsim))
samplenr <- as.numeric(rep(1:nsam, each = n * ntimes))
for (i in 1:nsim) {
    
    # Make a dataframe containing an identifier for the space--time sample, the
    # space--time coordinates, the simulated values, the stratum, and for the SP
    # design, the panel
    if (STDesign == "SP") {
        dat <- data.frame(samplenr, simdf$xsall, simdf$ysall, simdf$tsall, simdf[[i + 
            3]], data.frame(stsamples)$hsall, data.frame(stsamples)$panelall)
        names(dat) <- c("sam", "x", "y", "t", "z", "h", "panel")
    } else {
        dat <- data.frame(samplenr, simdf$xsall, simdf$ysall, simdf$tsall, simdf[[i + 
            3]], data.frame(stsamples)$hsall)
        names(dat) <- c("sam", "x", "y", "t", "z", "h")
    }
    
    # Estimate sampling variance-covariance matrix of estimated means
    Cp <- EstimateCp(dat = dat, stdesign = STDesign, nh = nh, nhSS = nhSS, experimental = FALSE, 
        wh = wh, period = period)
    
    # Add xi-covariance matrix to sampling covariance matrix of estimated
    # spatial means
    if (STDesign == "SP") {
        Xsup <- matrix(data = 0, nrow = 2 * ntimes, ncol = ntimes)
        for (j in 1:ntimes) {
            Xsup[j, j] <- 1
            Xsup[j + ntimes, j] <- 1
        }
        XCxiX <- Xsup %*% Cxi %*% t(Xsup)
        Cxip <- XCxiX + Cp
        invCxip <- solve(Cxip)
        DCDinv <- solve(t(Dsup) %*% invCxip %*% Dsup)
        
        # Estimate spatial means and regression coefficients for SP
        spatialMeans <- EstimateMean(dat = dat, stdesign = STDesign, wh = wh)
        b0GLS <- b1GLS <- numeric(length = nsam)
        for (samnr in 1:nsam) {
            spM <- as.numeric(spatialMeans[samnr, , ])
            DCY <- t(Dsup) %*% invCxip %*% spM
            beta <- DCDinv %*% DCY
            b0GLS[samnr] <- beta[1]
            b1GLS[samnr] <- beta[2]
        }
    } else {
        Cxip <- Cxi + Cp
        invCxip <- solve(Cxip)
        DCDinv <- solve(t(D) %*% invCxip %*% D)
        
        # Estimate spatial means and regression coefficients for IS, SS, SA
        spatialMeans <- EstimateMean(dat = dat, stdesign = STDesign, wh = wh)
        b0GLS <- b1GLS <- numeric(length = nsam)
        for (samnr in 1:nsam) {
            DCY <- t(D) %*% invCxip %*% as.numeric(spatialMeans[samnr, ])
            beta <- DCDinv %*% DCY
            b0GLS[samnr] <- beta[1]
            b1GLS[samnr] <- beta[2]
        }
    }
    # Estimate conditional sampling variances of regression coefficients
    Cpb[1, 1, i] <- var(b0GLS)
    Cpb[2, 2, i] <- var(b1GLS)
    Cpb[1, 2, i] <- Cpb[2, 1, i] <- cov(b0GLS, b1GLS)
}
```

### Results
Now compute the determinant of the covariance matrix of the model parameters, and the variance of the linear trend parameter.

```r
# Compute average of conditional sampling covariance matrices (averaged over
# simulations)
ExiCpb <- apply(Cpb, MARGIN = c(1, 2), FUN = mean)

# Add sampling- and model covariance matrices
Ctot <- ExiCpb + Cxibeta

# Compute determinants
print(detCxip <- Ctot[1, 1] * Ctot[2, 2] - Ctot[1, 2]^2)
```

```
## [1] 1.001e-12
```

```r
print(detCxi <- Cxibeta[1, 1] * Cxibeta[2, 2] - Cxibeta[1, 2]^2)
```

```
## [1] 6.39e-13
```

```r
print(detCp <- ExiCpb[1, 1] * ExiCpb[2, 2] - ExiCpb[1, 2]^2)
```

```
## [1] 3.535e-14
```

```r

# Compute variances
print(Vxipb1 <- Ctot[2, 2])
```

```
## [1] 3.599e-07
```

```r
print(Vxib1 <- Cxibeta[2, 2])
```

```
## [1] 2.654e-07
```

```r
print(ExiVpb1 <- ExiCpb[2, 2])
```

```
## [1] 9.454e-08
```

```r

print(result <- data.frame(STDesign, ntimes, n, detCxip, detCxi, detCp, Vxipb1, 
    Vxib1, ExiVpb1))
```

```
##   STDesign ntimes  n   detCxip   detCxi     detCp    Vxipb1     Vxib1
## 1       SP      4 30 1.001e-12 6.39e-13 3.535e-14 3.599e-07 2.654e-07
##     ExiVpb1
## 1 9.454e-08
```

