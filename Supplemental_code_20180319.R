##############################################################################
#
# A probabilistic model for soil populations of Fusarium culmorum in agricultural soil across the Inland Pacific Northwest based on local climate and future predictions under climate change.
# 
#   Section 1:    Loading in dataframes
#   Section 2:    Reference model code
#   Section 3:    Climate model code
#   Section 4:    Generating Figure 2
#   Section 5:    Generating Figure 3
#   Section 6:    Generating Figure 4
#   Section 7:    Generating Figure 5
#   Section 8:    Downloading downscaled GCM datasets
#   Section 9:    Calculating potential evapotranspiration
#   Section 10:   Generating Figure 6
#
# Author: Andrew Lloyd Robinson
#
##############################################################################

#dependent packages
library(rethinking)
library(zoo)
library(lubridate)
library(RNetCDF)

#Section 1: Loading in dataframes----

#climate frame of historical data for each field
af = read.csv(file = "hist_precip_evap_diff.csv", stringsAsFactors = FALSE)

#sampling date records
sampling.dates = as.Date(c("2016-06-11","2016-06-18","2016-06-25",
                           "2016-09-03","2016-09-10","2016-09-17",
                           "2016-12-03","2016-12-10","2016-12-15",
                           "2017-03-04","2017-03-11","2017-03-18",
                           "2017-06-03","2017-06-10","2017-06-17",
                           "2017-09-02","2017-09-09","2017-09-16",
                           "2017-12-02","2017-12-09","2017-12-16",
                           "2018-03-03","2018-03-10","2018-03-17"))

#converting into zoo object
df.zoo = read.zoo(af, format = "%Y-%m-%d")

#calculateing  rolling sums for seasonal term (prior 90 days)
sum.90 <- rollapply(df.zoo, 90, sum, align = c("right"))
df.final <- sum.90
df.final <- df.final[complete.cases(df.final),]
final.index = index(df.final)

#starting what will become the date sequence vector
ds <- 0

#looking up index values and adding to planting date vector
for (d in sampling.dates){
  pointer = which.min(abs(as.Date(d) - final.index))
  ds <- append(ds,pointer,after=length(ds))
}

#removing the 0 used at the beginning
ds <- ds[-1]

#converting into dataframe
bf <- data.frame(df.final)

#total precip 90 d prior
p.90 <- c(bf$p.1[ds[1]],bf$p.2[ds[1]],bf$p.3[ds[1]],bf$p.4[ds[2]],bf$p.5[ds[2]],bf$p.6[ds[2]],bf$p.7[ds[3]],bf$p.8[ds[3]],bf$p.9[ds[3]],
          bf$p.1[ds[4]],bf$p.2[ds[4]],bf$p.3[ds[4]],bf$p.4[ds[5]],bf$p.5[ds[5]],bf$p.6[ds[5]],bf$p.7[ds[6]],bf$p.8[ds[6]],bf$p.9[ds[6]],
          bf$p.1[ds[7]],bf$p.2[ds[7]],bf$p.3[ds[7]],bf$p.4[ds[8]],bf$p.5[ds[8]],bf$p.6[ds[8]],bf$p.7[ds[9]],bf$p.8[ds[9]],bf$p.9[ds[9]],
          bf$p.1[ds[10]],bf$p.2[ds[10]],bf$p.3[ds[10]],bf$p.4[ds[11]],bf$p.5[ds[11]],bf$p.6[ds[11]],bf$p.7[ds[12]],bf$p.8[ds[12]],bf$p.9[ds[12]],
          bf$p.1[ds[13]],bf$p.2[ds[13]],bf$p.3[ds[13]],bf$p.4[ds[14]],bf$p.5[ds[14]],bf$p.6[ds[14]],bf$p.7[ds[15]],bf$p.8[ds[15]],bf$p.9[ds[15]],
          bf$p.1[ds[16]],bf$p.2[ds[16]],bf$p.3[ds[16]],bf$p.4[ds[17]],bf$p.5[ds[17]],bf$p.6[ds[17]],bf$p.7[ds[18]],bf$p.8[ds[18]],bf$p.9[ds[18]],
          bf$p.1[ds[19]],bf$p.2[ds[19]],bf$p.3[ds[19]],bf$p.4[ds[20]],bf$p.5[ds[20]],bf$p.6[ds[20]],bf$p.7[ds[21]],bf$p.8[ds[21]],bf$p.9[ds[21]],
          bf$p.1[ds[22]],bf$p.2[ds[22]],bf$p.3[ds[22]],bf$p.4[ds[23]],bf$p.5[ds[23]],bf$p.6[ds[23]],bf$p.7[ds[24]],bf$p.8[ds[24]],bf$p.9[ds[24]])
#total evap 90 d prior
e.90 <- c(bf$evap.1[ds[1]],bf$evap.2[ds[1]],bf$evap.3[ds[1]],bf$evap.4[ds[2]],bf$evap.5[ds[2]],bf$evap.6[ds[2]],bf$evap.7[ds[3]],bf$evap.8[ds[3]],bf$evap.9[ds[3]],
          bf$evap.1[ds[4]],bf$evap.2[ds[4]],bf$evap.3[ds[4]],bf$evap.4[ds[5]],bf$evap.5[ds[5]],bf$evap.6[ds[5]],bf$evap.7[ds[6]],bf$evap.8[ds[6]],bf$evap.9[ds[6]],
          bf$evap.1[ds[7]],bf$evap.2[ds[7]],bf$evap.3[ds[7]],bf$evap.4[ds[8]],bf$evap.5[ds[8]],bf$evap.6[ds[8]],bf$evap.7[ds[9]],bf$evap.8[ds[9]],bf$evap.9[ds[9]],
          bf$evap.1[ds[10]],bf$evap.2[ds[10]],bf$evap.3[ds[10]],bf$evap.4[ds[11]],bf$evap.5[ds[11]],bf$evap.6[ds[11]],bf$evap.7[ds[12]],bf$evap.8[ds[12]],bf$evap.9[ds[12]],
          bf$evap.1[ds[13]],bf$evap.2[ds[13]],bf$evap.3[ds[13]],bf$evap.4[ds[14]],bf$evap.5[ds[14]],bf$evap.6[ds[14]],bf$evap.7[ds[15]],bf$evap.8[ds[15]],bf$evap.9[ds[15]],
          bf$evap.1[ds[16]],bf$evap.2[ds[16]],bf$evap.3[ds[16]],bf$evap.4[ds[17]],bf$evap.5[ds[17]],bf$evap.6[ds[17]],bf$evap.7[ds[18]],bf$evap.8[ds[18]],bf$evap.9[ds[18]],
          bf$evap.1[ds[19]],bf$evap.2[ds[19]],bf$evap.3[ds[19]],bf$evap.4[ds[20]],bf$evap.5[ds[20]],bf$evap.6[ds[20]],bf$evap.7[ds[21]],bf$evap.8[ds[21]],bf$evap.9[ds[21]],
          bf$evap.1[ds[22]],bf$evap.2[ds[22]],bf$evap.3[ds[22]],bf$evap.4[ds[23]],bf$evap.5[ds[23]],bf$evap.6[ds[23]],bf$evap.7[ds[24]],bf$evap.8[ds[24]],bf$evap.9[ds[24]])

#adding historic onservations to dataframe
df$p.h.90 <- rep(p.90,each=9)
df$e.h.90 <- rep(e.90,each=9)

#calculating difference
df$p.d.90 <- df$p.h.90-df$e.h.90

#historic grand mean and sd
historic_diff_mu <- mean(as.matrix(bf[,19:27]))
historic_diff_sd <- sd(as.matrix(bf[,19:27]))

#standardizing atmospheric water balance using grand mean and standard deviation
df$prior_diff <- (df$p.d.90-historic_diff_mu)/historic_diff_sd

#adding a reference column
df$number <- as.integer(rep( seq( 1 , 72 , 1 ), each=9 ))

#replacing 1 quadrat in field 9 with the LOD (50 F. culmorum ppg) to include field 9 in the analysis
#field 9 December 2016 had no detectable F. culmorum in all 9 quadrats sampled
df$f.c.ppg[235] = 50
#field 9 December 2017 had no detectable F. culmorum in all 9 quadrats sampled
df$f.c.ppg[559] = 50

#Section 2: Reference model code----

#list for Stan
test.list.1 <- list(ppg = df$f.c.ppg , field = df$field , sample = df$number)

#model code
m.field.sample <- map2stan(
  alist(
    ppg ~ dexp(lambda),
    log(lambda) <- a + a_field[field] + a_sample[sample],
    a_field[field] ~ dnorm(0,sigma_field),
    a_sample[sample] ~ dnorm(0,sigma_sample),
    a ~ dnorm(0,1),
    sigma_field ~ dexp(1),
    sigma_sample ~ dexp(1)
  ) ,
  data = test.list.1 ,
  warmup = 1000 ,
  iter = 3500 ,
  chains = 4 ,
  cores = 4 ,
  control = list(adapt_delta = 0.99)
)

#parameter summary
precis(m.field.sample,depth=2)
#visualization
plot(precis(m.field.sample,depth=2))

#extracting samples
e.1 <- extract.samples(m.field.sample,n=10000)

#Section 3: Climate model code----

#list for model
test.list.2 <- list(
  ppg = df$f.c.ppg,
  field = as.numeric(df$field),
  prior_diff = df$prior_diff)

#model equation
m.1 <- map2stan(
  alist(
    ppg ~ dexp(lambda),
    log(lambda) <- A + B*(prior_diff) ,
    A <- a + a_field[field] ,
    B <- b + b_field[field] ,
    c(a_field,b_field)[field] ~ dmvnorm2(0, sigma_field , Rho ),
    a ~ dnorm(-2.87,0.99) ,
    b ~ dnorm(0,1) ,
    sigma_field ~ dcauchy(0, 1) ,
    Rho ~ dlkjcorr(7)
  ) ,
  data =  test.list.2 ,
  warmup = 1000 ,
  iter = 3500 ,
  chains = 4 ,
  cores = 4 ,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

#parameter summary
precis(m.1,prob=0.95,depth=2)

#extracting samples
e.2 <- extract.samples(m.1,n=10000)

#Section 4: Generating Figure 2----

#collecting reference model output
link.rm <- link(m.field.sample,n=10000)
link.cm <- link(m.1,n=10000)

#only keeping one set of rates per sampling iteration
sim.rates.rm <- link.rm[,seq(1,648,9)]
sim.rates.cm <- link.cm$lambda[,seq(1,648,9)]

#name sequence for fields
ns <- c("Field 1","Field 2","Field 3","Field 4","Field 5","Field 6","Field 7","Field 8","Field 9",
        rep(" ",63))

#season sequence
seasons <- as.vector(c("Jun 2016",rep(" ",8),
                       "Sep 2016",rep(" ",8),
                       "Dec 2016",rep(" ",8),
                       "Mar 2017",rep(" ",8),
                       "Jun 2017",rep(" ",8),
                       "Sep 2017",rep(" ",8),
                       "Dec 2017",rep(" ",8),
                       "Mar 2018",rep(" ",8)))

#logical sequence for x axis labels
xs <- as.vector(c(rep(FALSE,8),TRUE,
                  rep(FALSE,8),TRUE,
                  rep(FALSE,8),TRUE,
                  rep(FALSE,8),TRUE,
                  rep(FALSE,8),TRUE,
                  rep(FALSE,8),TRUE,
                  rep(FALSE,8),TRUE,
                  rep(FALSE,8),TRUE,
                  rep(FALSE,8),TRUE))

#logical sequence for y axis labels
ys <- as.vector(c(rep(TRUE,9),rep(FALSE,63)))

#index sequence for raw data
fs <- seq(1,648,9)

#plotting
par(oma=c(5,5,1,1))
par(mar=c(0.5,0.5,0.5,0.5))
par(mfcol=c(9,8))
for(n in 1:72){
  plot(0,type="n",xlim=c(0,1),ylim=c(0,7500),xlab=" ",ylab=" ",axes=FALSE)
  box()
  axis(1,at=seq(0,1,0.25),las=2,labels=xs[n])
  axis(2,at=seq(0,7500,2500),las=1,labels=ys[n])
  for(i in 1:100){
    lines(rev(seq(0.01,0.99,0.01)),-log(1-seq(0.01,0.99,0.01))/sim.rates.rm[i,n],col=col.alpha("red",0.15))
  }
  for(i in 1:100){
    lines(rev(seq(0.01,0.99,0.01)),-log(1-seq(0.01,0.99,0.01))/sim.rates.cm[i,n],col=col.alpha("green3",0.15))
  }
  x= df$f.c.ppg[fs[n]:(fs[n]+8)]
  points((1:length(x)-0.5)/length(x),sort(x,decreasing=TRUE),type="p",pch=16,col="black")
  mtext(ns[n],line=-1.5,side=3)
  mtext(seasons[n],line=0,side=3)
}
mtext("Fraction of field exceeded", side = 1, outer = TRUE, line = 3)
mtext("PPG", side = 2, outer = TRUE, line = 3)

#Section 5: Generating Figure 3----

#index sequence for raw data
fs <- seq(1,648,9)

for(n in 1:72){
  x= df$f.c.ppg[fs[n]:(fs[n]+8)]
  sample <- sort(x,decreasing=FALSE)
  df$sorted.ppg[fs[n]:(fs[n]+8)] <- sample
}


df$quantile <- rep((seq(1:9)-0.5)/9,72)

ref_sample_1 <- link(m.reference,10000,data=test.list.1)
ref_sample_2 <- link(m.1,10000,data=test.list.2)

post_ppg_1 <- matrix(NA,10000,648)
post_ppg_2 <- matrix(NA,10000,648)

for(n in 1:648){
  post_ppg_1[,n] <- -log(1-df$quantile[n])/ref_sample_1[,n]
}
for(n in 1:648){
  post_ppg_2[,n] <- -log(1-df$quantile[n])/ref_sample_2$lambda[,n]
}

post_ppg_quant_1 <- apply(post_ppg_1,2,QUANT)
post_ppg_quant_2 <- apply(post_ppg_2,2,QUANT)

par(mfrow=c(1,2))
plot(0,type="n",xlim=c(0,5500),ylim=c(0,5500),xlab="Measured PPG",ylab="Predicted PPG at same frequency of exceedance",main="A",axes=FALSE)
box()
axis(1,at=seq(0,6000,1000),las=1)
axis(2,at=seq(0,6000,1000),las=2)
arrows(df$sorted.ppg, post_ppg_quant_1[2,], df$sorted.ppg, post_ppg_quant_1[4,], length=0.05, angle=90, code=3,col="red")
abline(a=0,b=1)


plot(0,type="n",xlim=c(0,5500),ylim=c(0,5500),xlab="Measured PPG",ylab="Predicted PPG at same frequency of exceedance",main="B",axes=FALSE)
box()
axis(1,at=seq(0,6000,1000),las=1)
axis(2,at=seq(0,6000,1000),las=2)
arrows(df$sorted.ppg, post_ppg_quant_2[2,], df$sorted.ppg, post_ppg_quant_2[4,], length=0.05, angle=90, code=3,col="green3")
abline(a=0,b=1)

measured_median <- sapply(split(df$f.c.ppg,df$number),median)


#Section 6: Generating Figure 4----

#distribution of expected values for the population average
meta_ev <- log(2)/exp(e.2$a)

#calculating distributions of posterior medians for each field (N=9)
for(n in 1:9){
  assign(paste0("field_",n,"_ev"),log(2)/exp(e.2$a + e.2$a_field[,n]))
}

#adding a year column to the historical climate data dataframe
af$Y <- lubridate::year(af$Date)

#annual precipitation by year
p.1 <- sapply(split(af$p.1,af$Y),sum)
p.2 <- sapply(split(af$p.2,af$Y),sum)
p.3 <- sapply(split(af$p.3,af$Y),sum)
p.4 <- sapply(split(af$p.4,af$Y),sum)
p.5 <- sapply(split(af$p.5,af$Y),sum)
p.6 <- sapply(split(af$p.6,af$Y),sum)
p.7 <- sapply(split(af$p.7,af$Y),sum)
p.8 <- sapply(split(af$p.8,af$Y),sum)
p.9 <- sapply(split(af$p.9,af$Y),sum)

#annual potential exapotranspiration by year
evap.1 <- sapply(split(af$evap.1,af$Y),sum)
evap.2 <- sapply(split(af$evap.2,af$Y),sum)
evap.3 <- sapply(split(af$evap.3,af$Y),sum)
evap.4 <- sapply(split(af$evap.4,af$Y),sum)
evap.5 <- sapply(split(af$evap.5,af$Y),sum)
evap.6 <- sapply(split(af$evap.6,af$Y),sum)
evap.7 <- sapply(split(af$evap.7,af$Y),sum)
evap.8 <- sapply(split(af$evap.8,af$Y),sum)
evap.9 <- sapply(split(af$evap.9,af$Y),sum)

#plotting the results
par(mfrow=c(1,3))
par(cex=1.25)

#Figure 5A
plot(0,type="n",xlim=c(1,10),ylim=c(0,1000),xlab="Field",ylab="PPG",main="A",axes=FALSE)
box()
axis(1,at=seq(1,10,1),las=1,labels=c("A",1,2,3,4,5,6,7,8,9))
axis(2,at=seq(0,1000,100),las=1)
points(x=rep(1,100),y=sample(meta_ev,100),col=(col.alpha("black",0.1)),pch=15,cex=2)
points(x=rep(2,100),y=sample(field_1_ev,100),col=(col.alpha("green3",0.1)),pch=15,cex=2)
points(x=rep(3,100),y=sample(field_2_ev,100),col=(col.alpha("green3",0.1)),pch=15,cex=2)
points(x=rep(4,100),y=sample(field_3_ev,100),col=(col.alpha("green3",0.1)),pch=15,cex=2)
points(x=rep(5,100),y=sample(field_4_ev,100),col=(col.alpha("green3",0.1)),pch=15,cex=2)
points(x=rep(6,100),y=sample(field_5_ev,100),col=(col.alpha("green3",0.1)),pch=15,cex=2)
points(x=rep(7,100),y=sample(field_6_ev,100),col=(col.alpha("green3",0.1)),pch=15,cex=2)
points(x=rep(8,100),y=sample(field_7_ev,100),col=(col.alpha("green3",0.1)),pch=15,cex=2)
points(x=rep(9,100),y=sample(field_8_ev,100),col=(col.alpha("green3",0.1)),pch=15,cex=2)
points(x=rep(10,100),y=sample(field_9_ev,100),col=(col.alpha("green3",0.1)),pch=15,cex=2)

#Figure 5B
plot(0,type="n",xlim=c(1,9),ylim=c(350,900),
     xlab="Field",ylab="Annual precipitation (mm)",main="B",axes=FALSE)
box()
axis(1,at=seq(1,9,1),las=1,labels=c(1,2,3,4,5,6,7,8,9))
axis(2,at=seq(350,900,50),las=1)
for (n in 1:39){
  lines(x=seq(1,9,1),y=c(p.1[n],p.2[n],p.3[n],p.4[n],p.5[n],p.6[n],p.7[n],p.8[n],p.9[n]),
        col=(col.alpha("blue",0.25)),pch=15,cex=2)
}

#Figure 5C
plot(0,type="n",xlim=c(1,9),ylim=c(-1150,-900),
     xlab="Field",ylab="Annual potential evapotranspiration (-mm)",main="C",axes=FALSE)
box()
axis(1,at=seq(1,9,1),las=1,labels=c(1,2,3,4,5,6,7,8,9))
axis(2,at=seq(-1150,900,50),las=1)
for (n in 1:39){
  lines(x=seq(1,9,1),y=-c(evap.1[n],evap.2[n],evap.3[n],evap.4[n],evap.5[n],evap.6[n],evap.7[n],evap.8[n],evap.9[n]),
        col=(col.alpha("red",0.25)),pch=15,cex=2)
}

#Section 7: Generating Figure 5----


#function to plot predicted trends by field (n)
TREND_SIMULATOR <- function(n){
  #finding relevant pieces of model
  intercept <- sample(e.2$a+e.2$a_field[,n],100)
  slope <- sample(e.2$b+e.2$b_field[,n],100)
  #sample sequence
  range <- seq(-2,2,0.5)
  reference <- historic_diff_mu+range*historic_diff_sd
  ev_matrix <- matrix(0,9,100)
  for (i in 1:9){
    for (j in 1:100){
      ev_matrix[i,j] <- log(2)/exp(intercept[j] + slope[j]*range[i])
    }
  }
  plot(0,type="n",xlim=c(-2,2),ylim=c(0,1500),xlab="Standardized prior 90 day P - PET",ylab="PPG",
       main=paste("Field",n),axes=FALSE)
  box()
  #axis(1,at=seq(-2,2,0.5),las=1,labels=round(reference,digits=0))
  axis(1,at=seq(-2,2,0.5),las=1)
  axis(2,at=seq(0,1500,250),las=2)
  for(j in 1:100){
    lines(x=range,y=ev_matrix[,j],col=col.alpha("green3",0.2))
  }
  slope_mu <- round(mean(e.2$b+e.2$b_field[,n]),digits=2)
  slope_sd <- round(sd(e.2$b+e.2$b_field[,n]),digits=2)
  mtext(paste("B =",slope_mu,"±",slope_sd),side=3,line=-2)
}

#generating plot
par(mfrow=c(3,3))
for (n in 1:9){
  TREND_SIMULATOR(n)
}



#Section 8: Downloading downscaled GCM datasets----

#dataframe containing GPS coordinates of all sampled fields
locations <- read.csv("gps_ref.csv")

#column for naming files
locations$town <- c("FIELD_1","FIELD_2","FIELD_3","FIELD_4","FIELD_5","FIELD_6","FIELD_7","FIELD_8","FIELD_9")

#number of sites in location file
N <- nrow(locations)

#adding elevation data
elevation_grid <- open.nc("metdata_elevationdata.nc")
elevation_ref <- var.get.nc(elevation_grid,variable=2)

#adding elevation column
locations$elevation <- rep(NA,N)

#adding elevation values
for (n in 1:N){
  x <- locations$lon[n]
  y <- locations$lat[n]
  lat <- var.get.nc(elevation_grid,"lat")
  lon <- var.get.nc(elevation_grid,"lon")
  flat = match(abs(lat - y) < 1/48, 1)
  latindex = which(flat %in% 1)
  flon = match(abs(lon - x) < 1/48, 1)
  lonindex = which(flon %in% 1)
  locations$elevation[n] <- 0.1*elevation_ref[lonindex,latindex]
}

#base URL
url_1 <- "http://thredds.northwestknowledge.net:8080/thredds/dodsC/NWCSC_INTEGRATED_SCENARIOS_ALL_CLIMATE/macav2livneh/"

#model list
model <- c("bcc-csm1-1", "bcc-csm1-1-m","BNU-ESM","CanESM2",
           "CSIRO-Mk3-6-0","GFDL-ESM2G","GFDL-ESM2M","HadGEM2-CC365",
           "HadGEM2-ES365","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",
           "IPSL-CM5B-LR","MIROC5","MIROC-ESM","MIROC-ESM-CHEM",
           "MRI-CGCM3","NorESM1-M","CNRM-CM5","CCSM4")

#model condition list becuase CCSM4 is different
condition <- c(rep("_r1i1p1_",19),"_r6i1p1_")

#timestep lists
timestep_historical <- c("historical_1950_1969_CONUS_daily.nc",
                         "historical_1970_1989_CONUS_daily.nc",
                         "historical_1990_2005_CONUS_daily.nc")
timestep_45 <- c("rcp45_2006_2025_CONUS_daily.nc",
                 "rcp45_2026_2045_CONUS_daily.nc",
                 "rcp45_2046_2065_CONUS_daily.nc",
                 "rcp45_2066_2085_CONUS_daily.nc",
                 "rcp45_2086_2099_CONUS_daily.nc")
timestep_85 <- c("rcp85_2006_2025_CONUS_daily.nc",
                 "rcp85_2026_2045_CONUS_daily.nc",
                 "rcp85_2046_2065_CONUS_daily.nc",
                 "rcp85_2066_2085_CONUS_daily.nc",
                 "rcp85_2086_2099_CONUS_daily.nc")

#variable list
variable_1 <- c("_tasmax_","_tasmin_","_pr_","_huss_","_was_","_rsds_")

#reference longitude and latitude sequence
lon <- seq(235.40625,292.96875,0.0625)
lat <- seq(25.15625,52.84375,0.0625)

#building 3-d arrays [variable, model, time] for 9 fields for 1950-2005
for (n in 1:N){
  #assinging lon and lat from csv
  x <- locations$lon[n]
  y <- locations$lat[n]
  #changing into coordinates 
  coord <- c(360+x,y)
  #locating appropriate latitude index
  flat = match(abs(lat - coord[2]) < 1/32, 1)
  latindex = which(flat %in% 1)
  #locating appropriate longitude index
  flon = match(abs(lon - coord[1]) < 1/32, 1)
  lonindex = which(flon %in% 1)
  #creating a blank 3-d array
  db_hist <- array(dim=c(20,20,20454))
  #looping through variables and models
  for(i in 1:6){
    for(j in 1:20){
      nc_1 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_historical[1]))
      var_1 <- as.numeric(var.get.nc(nc_1, variable=4,
                                     start=c(lonindex, latindex, 1),
                                     count=c(1,1,7305)))
      nc_2 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_historical[2]))
      var_2 <- as.numeric(var.get.nc(nc_2, variable=4,
                                     start = c(lonindex, latindex, 1),
                                     count=c(1,1,7305)))
      nc_3 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_historical[3]))
      var_3 <- as.numeric(var.get.nc(nc_3, variable=4,
                                     start=c(lonindex, latindex, 1),
                                     count=c(1,1,5844)))
      #combining results
      var <- c(var_1,var_2,var_3)
      #filling the 3-d array
      db_hist[i,j,] <- var
      #crudely printing the progress
      print(paste(locations$town[n],variable_1[i],model[j],"historic simulation completed"))
    }
  }
  #saving the 3-d array with the appropriate town name
  assign(paste0(locations$town[n],"_hist"),db_hist)
}

#building 3-d arrays [variable, model, time] for 9 fields for RCP 4.5
for (n in 1:1.5){
  #assinging lon and lat from csv
  x <- locations$lon[n]
  y <- locations$lat[n]
  #changing into coordinates 
  coord <- c(360+x,y)
  #locating appropriate latitude index
  flat = match(abs(lat - coord[2]) < 1/32, 1)
  latindex = which(flat %in% 1)
  #locating appropriate longitude index
  flon = match(abs(lon - coord[1]) < 1/32, 1)
  lonindex = which(flon %in% 1)
  #creating a blank 3-d array
  db_45 <- array(dim=c(20,20,34333))
  #looping through variables and models
  for(i in 1:6){
    for(j in 1:20){
      nc_1 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_45[1]))
      var_1 <- as.numeric(var.get.nc(nc_1, variable=4,
                                     start=c(lonindex, latindex, 1),
                                     count=c(1,1,7305)))
      nc_2 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_45[2]))
      var_2 <- as.numeric(var.get.nc(nc_2, variable=4,
                                     start=c(lonindex, latindex, 1),
                                     count=c(1,1,7305)))
      nc_3 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_45[3]))
      var_3 <- as.numeric(var.get.nc(nc_3, variable=4,
                                     start=c(lonindex, latindex, 1),
                                     count=c(1,1,7305)))
      nc_4 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_45[4]))
      var_4 <- as.numeric(var.get.nc(nc_4, variable=4,
                                     start=c(lonindex, latindex, 1),
                                     count=c(1,1,7305)))
      nc_5 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_45[5]))
      var_5 <- as.numeric(var.get.nc(nc_5, variable=4,
                                     start=c(lonindex, latindex, 1),
                                     count=c(1,1,5113)))
      #combining results
      var <- c(var_1,var_2,var_3,var_4,var_5)
      #filling the 3-d array
      db_45[i,j,] <- var
      #crudely printing the progress
      print(paste(locations$town[n],variable_1[i],model[j],"RCP 4.5 completed"))
    }
  }
  #saving the 3-d array with the appropriate town name
  assign(paste0(locations$town[n],"_45"),db_45)
}

#building 3-d arrays [variable, model, time] for 9 fields for RCP 8.5
for (n in 1:N){
  #assinging lon and lat from csv
  x <- locations$lon[n]
  y <- locations$lat[n]
  #changing into coordinates 
  coord <- c(360+x,y)
  #locating appropriate latitude index
  flat = match(abs(lat - coord[2]) < 1/32, 1)
  latindex = which(flat %in% 1)
  #locating appropriate longitude index
  flon = match(abs(lon - coord[1]) < 1/32, 1)
  lonindex = which(flon %in% 1)
  #creating a blank 3-d array
  db_85 <- array(dim=c(20,20,34333))
  #looping through variables and models
  for(i in 1:6){
    for(j in 1:20){
      nc_1 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_85[1]))
      var_1 <- as.numeric(var.get.nc(nc_1, variable=4,
                                     start=c(lonindex, latindex, 1),
                                     count=c(1,1,7305)))
      nc_2 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_85[2]))
      var_2 <- as.numeric(var.get.nc(nc_2, variable=4,
                                     start=c(lonindex, latindex, 1),
                                     count=c(1,1,7305)))
      nc_3 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_85[3]))
      var_3 <- as.numeric(var.get.nc(nc_3, variable=4,
                                     start=c(lonindex, latindex, 1),
                                     count=c(1,1,7305)))
      nc_4 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_85[4]))
      var_4 <- as.numeric(var.get.nc(nc_4, variable=4,
                                     start=c(lonindex, latindex, 1),
                                     count=c(1,1,7305)))
      nc_5 <- open.nc(paste0(url_1,model[j],"/macav2livneh",variable_1[i],
                             model[j],condition[j],timestep_85[5]))
      var_5 <- as.numeric(var.get.nc(nc_5, variable=4,
                                     start=c(lonindex, latindex, 1),
                                     count=c(1,1,5113)))
      #combining results
      var <- c(var_1,var_2,var_3,var_4,var_5)
      #filling the 3-d array
      db_85[i,j,] <- var
      #crudely printing the progress
      print(paste(locations$town[n],variable_1[i],model[j],"RCP 8.5 completed"))
    }
  }
  #saving the 3-d array with the appropriate town name
  assign(paste0(locations$town[n],"_85"),db_85)
}

#saving the 3-d arrays (historical simulation)
for(n in 1:N){
  saveRDS(get(paste0(locations$town[n],"_hist")),
          paste0(locations$town[n],"_hist.rdata"))
  print(paste(locations$town[n],"saved"))
}

#saving the 3-d arrays (RCP 4.5)
for(n in 1:N){
  saveRDS(get(paste0(locations$town[n],"_45")),
          paste0(locations$town[n],"_45.rdata"))
  print(paste(locations$town[n],"saved"))
}

#saving the 3-d arrays (RCP 8.5)
for(n in 1:N){
  saveRDS(get(paste0(locations$town[n],"_85")),
          paste0(locations$town[n],"_85.rdata"))
  print(paste(locations$town[n],"saved"))
}

#loading the 3-d arrays (historical simulation)
for(n in 1:N){
  assign(paste0(locations$town[n],"_hist"),
         readRDS(paste0(locations$town[n],"_hist.rdata")))
  print(paste(locations$town[n],"historical simulation loaded"))
}

#loading the 3-d arrays (RCP 4.5)
for(n in 1:N){
  assign(paste0(locations$town[n],"_45"),
         readRDS(paste0(locations$town[n],"_45.rdata")))
  print(paste(locations$town[n],"RCP 4.5 loaded"))
}

#loading the 3-d arrays (RCP 8.5)
for(n in 1:N){
  assign(paste0(locations$town[n],"_85"),
         readRDS(paste0(locations$town[n],"_85.rdata")))
  print(paste(locations$town[n],"RCP 8.5 loaded"))
}

#Section 9: Calculating potential evapotranspiration----

#Reference grass potential evapotranspiration
ETo <- function(elevation,max_temp,min_temp,specific_humidity,wind_speed,rad,J,lat){
  #atmospheric pressure based on site elevation (kPa)
  atmo_pressure <- 101.3*((293-0.0065*elevation)/293)^5.26
  #psychrometric constant
  psy_constant<- 0.000665*atmo_pressure
  #mean air temperature (C)
  avg_temp <- ((max_temp-273.15)+(min_temp-273.15))/2
  #converting surface downwelling shortwave radiation from W/m^2 to MJ/m^2 per day
  solar_rad <- 0.0864*rad
  #slope of the saturation vapor pressure-temperature curve
  delta <- (2503*exp((17.27*avg_temp)/(avg_temp+237.3)))/(avg_temp+237.3)^2
  #calculating dew temp
  dew_temp <- (((1/273.15)-(1.844*10^-4)*log((specific_humidity*atmo_pressure/0.622)/0.6113))^-1)-273.15
  #wind speed at 2 m above ground (measured from 10 m)
  wind <- wind_speed*(4.87/(log(67.8*10-5.42)))
  #inverse relative distance factor
  dist_factor <- 1+0.033*cos((2*pi/365)*J)
  #solar declination
  solar_dec <- 0.409*sin(((2*pi/365)*J)-1.39)
  #sunset hour angle
  sunset_hour <- acos(-tan(lat)*tan(solar_dec))
  #extraterrestrial radiation for 24 hour period
  extra_rad <- (24/pi)*4.92*dist_factor*(sunset_hour*sin(lat)*sin(solar_dec)+cos(lat)*cos(solar_dec)*sin(sunset_hour))
  #clear sky radiation
  clrsky_solar_rad <- (0.75+(2*10^-5)*elevation)*(extra_rad)
  #cloudiness function
  #bound (solar_rad/clrsky_solar_rad) ratio between 0.3 and 1
  cloudy <- 1.35*(ifelse((solar_rad/clrsky_solar_rad)<0.3,0.3,ifelse((solar_rad/clrsky_solar_rad)>1,1,(solar_rad/clrsky_solar_rad))))-0.35
  #saturation vapor pressure
  Es <- ((0.6108*exp((17.27*(max_temp-273.15))/((max_temp-273.15)+237.3)))+(0.6108*exp((17.27*(min_temp-273.15))/((min_temp-273.15)+237.3))))/2
  #actual vapor pressure
  Ea <- (0.6108*exp((17.27*(dew_temp))/((dew_temp)+237.3)))
  #net long wave radiation
  nlw_rad <- (4.901*10^-9)*cloudy*(0.34-0.14*sqrt(Ea))*((max_temp^4+min_temp^4)/2)
  #net shortwave radiation: albedo fixed at 0.23
  nsw_rad <- (1-0.23)*solar_rad
  #net radiation
  net_rad <- nsw_rad-nlw_rad
  #final equation, if negative output is 0
  ifelse((0.408*delta*net_rad + psy_constant*(900/(avg_temp+273))*wind*(Es-Ea))/(delta + psy_constant*(1+0.34*wind))<0,0,
         (0.408*delta*net_rad + psy_constant*(900/(avg_temp+273))*wind*(Es-Ea))/(delta + psy_constant*(1+0.34*wind)))
}

#Section 10: Generating Figure 6----

#creating the array
historical_diff_array <- array(NA,dim=c(9,20,12,56))

#looping through fields (n) and months (m)
for(n in 1:9){
  #reference date sequence
  future_dates <- as.Date(seq(38716,73048,1),origin="1900-01-01")
  historical_dates <- as.Date(seq(18262,38715,1),origin="1900-01-01")
  #date sequences as day of the year
  historical_doy <- as.integer(format(historical_dates,"%j"))
  future_doy <- as.integer(format(future_dates,"%j"))
  #month label for plots
  month_lab <- c("January","February","March","April","May","June","July","August","September","October","November","December")
  #copying the array of interest
  A <- get(paste0(locations$town[n],"_hist"))
  #adding ETo dimension
  for(j in 1:20){
    #calculating Eto by model
    A[7,j,] <- ETo(locations$elevation[n],
                   A[1,j,], A[2,j,], A[4,j,], A[5,j,], A[6,j,],
                   J=historical_doy,lat=(pi/180)*locations$lat[n])
    #placing the results back into the array for each town
    assign(paste0(locations$town[n],"_hist"),A)
    print(paste(locations$town[n],"historical simulation now has ETo for",model[j]))
  }
  
  #calculating sums for precip [var=8]
  for (j in 1:20){
    model_zoo <- as.zoo(A[3,j,])
    model_sum <- rollapply(model_zoo,90,sum,align=c("right"))
    as.vector(model_sum)
    A[8,j,] <- c(rep(0,89),as.vector(model_sum))
    print(paste("Model",j,"completed for precip"))
  }
  #calculating sums for precip [var=9]
  for (j in 1:20){
    model_zoo <- as.zoo(A[7,j,])
    model_sum <- rollapply(model_zoo,90,sum,align=c("right"))
    as.vector(model_sum)
    A[9,j,] <- c(rep(0,89),as.vector(model_sum))
    print(paste("Model",j,"completed for evap"))
  }
  #looping through months
  for(m in 1:12){
    #finding relevant pieces of model
    #intercept <- sample(e.2$a + e.2$a_field[,n],100)
    #slope <- sample(e.2$b + e.2$b_field[,n],100)
    #finding relevant pieces of model using same samples throughout
    intercept <- e.2$a[1:1000] + e.2$a_field[1:1000,n]
    slope <- e.2$b[1:1000] + e.2$b_field[1:1000,n]
    #looping through models
    for (j in 1:20){
      #starting data frame
      hist_frame <- as.data.frame(historical_dates)
      #lubridating by month
      hist_frame$M <- lubridate::month(hist_frame$historical_dates)
      #lubridating by year
      hist_frame$Y <- lubridate::year(hist_frame$historical_dates)
      hist_frame$P <- A[8,j,]
      hist_frame$ETo <- A[9,j,]
      #only keeping relevant month
      #only keeping relevant month
      hist_frame <- hist_frame[hist_frame$M==m,]
      #summarizing results by year
      mean_precip <- sapply(split(hist_frame$P,hist_frame$Y),mean)
      mean_evap <- sapply(split(hist_frame$ETo,hist_frame$Y),mean)
      mean_diff <- ((mean_precip-mean_evap)-historic_diff_mu)/historic_diff_sd
      
      historical_diff_array[n,j,m,] <- mean_diff
    }
    #printing progress
    print(paste(n,m,"done"))
  }
}

diff_45_array <- array(NA,dim=c(9,20,12,93))
diff_85_array <- array(NA,dim=c(9,20,12,93))

#looping through fields (n) and months (m)
for(n in 1:9){
  #reference date sequence
  future_dates <- as.Date(seq(38716,73048,1),origin="1900-01-01")
  historical_dates <- as.Date(seq(18262,38715,1),origin="1900-01-01")
  #date sequences as day of the year
  historical_doy <- as.integer(format(historical_dates,"%j"))
  future_doy <- as.integer(format(future_dates,"%j"))
  #month label for plots
  month_lab <- c("January","February","March","April","May","June","July","August","September","October","November","December")
  #copying the array of interest
  B <- get(paste0(locations$town[n],"_45"))
  C <- get(paste0(locations$town[n],"_85"))
  #looping through the models
  for(j in 1:20){
    #calculating GDD by model
    B[7,j,] <- ETo(locations$elevation[n],
                   B[1,j,], B[2,j,], B[4,j,], B[5,j,], B[6,j,],
                   J=future_doy,lat=(pi/180)*locations$lat[n])
    #placing the results back into the array for each town
    assign(paste0(locations$town[n],"_45"),A)
    print(paste(locations$town[n],"RCP 4.5 now has ETo for",model[j]))
  }
  #looping through the models
  for(j in 1:20){
    #calculating GDD by model
    C[7,j,] <- ETo(locations$elevation[n],
                   C[1,j,], C[2,j,], C[4,j,], C[5,j,], C[6,j,],
                   J=future_doy,lat=(pi/180)*locations$lat[n])
    #placing the results back into the array for each town
    assign(paste0(locations$town[n],"_45"),C)
    print(paste(locations$town[n],"RCP 8.5 now has ETo for",model[j]))
  }
  #calculating sums for precip [var=8]
  for (j in 1:20){
    model_zoo <- as.zoo(B[3,j,])
    model_sum <- rollapply(model_zoo,90,sum,align=c("right"))
    as.vector(model_sum)
    B[8,j,] <- c(rep(0,89),as.vector(model_sum))
    model_zoo <- as.zoo(C[3,j,])
    model_sum <- rollapply(model_zoo,90,sum,align=c("right"))
    as.vector(model_sum)
    C[8,j,] <- c(rep(0,89),as.vector(model_sum))
    print(paste("Model",j,"completed for precip"))
  }
  #calculating sums for precip [var=9]
  for (j in 1:20){
    model_zoo <- as.zoo(B[7,j,])
    model_sum <- rollapply(model_zoo,90,sum,align=c("right"))
    as.vector(model_sum)
    B[9,j,] <- c(rep(0,89),as.vector(model_sum))
    model_zoo <- as.zoo(C[7,j,])
    model_sum <- rollapply(model_zoo,90,sum,align=c("right"))
    as.vector(model_sum)
    C[9,j,] <- c(rep(0,89),as.vector(model_sum))
    print(paste("Model",j,"completed for evap"))
  }
  #looping through months
  for(m in 1:12){
    #finding relevant pieces of model
    #intercept <- sample(e.2$a + e.2$a_field[,n],100)
    #slope <- sample(e.2$b + e.2$b_field[,n],100)
    #finding relevant pieces of model using same samples throughout
    intercept <- e.2$a[1:1000] + e.2$a_field[1:1000,n]
    slope <- e.2$b[1:1000] + e.2$b_field[1:1000,n]
    #looping through models
    for (j in 1:20){
      #starting data frame
      fut_frame_45 <- as.data.frame(future_dates)
      #lubridating by month
      fut_frame_45$M <- lubridate::month(fut_frame_45$future_dates)
      #lubridating by year
      fut_frame_45$Y <- lubridate::year(fut_frame_45$future_dates)
      fut_frame_45$P <- B[8,j,]
      fut_frame_45$ETo <- B[9,j,]
      #only keeping relevant month
      #only keeping relevant month
      fut_frame_45 <- fut_frame_45[fut_frame_45$M==m,]
      #summarizing results by year
      mean_precip_45 <- sapply(split(fut_frame_45$P,fut_frame_45$Y),mean)
      mean_evap_45 <- sapply(split(fut_frame_45$ETo,fut_frame_45$Y),mean)
      mean_diff_45 <- ((mean_precip_45-mean_evap_45)-historic_diff_mu)/historic_diff_sd
      
      diff_45_array[n,j,m,] <- as.vector(mean_diff_45[1:93])
      
      #starting data frame
      fut_frame_85 <- as.data.frame(future_dates)
      #lubridating by month
      fut_frame_85$M <- lubridate::month(fut_frame_85$future_dates)
      #lubridating by year
      fut_frame_85$Y <- lubridate::year(fut_frame_85$future_dates)
      fut_frame_85$P <- C[8,j,]
      fut_frame_85$ETo <- C[9,j,]
      #only keeping relevant month
      #only keeping relevant month
      fut_frame_85 <- fut_frame_85[fut_frame_85$M==m,]
      #summarizing results by year
      mean_precip_85 <- sapply(split(fut_frame_85$P,fut_frame_85$Y),mean)
      mean_evap_85 <- sapply(split(fut_frame_85$ETo,fut_frame_85$Y),mean)
      mean_diff_85 <- ((mean_precip_85-mean_evap_85)-historic_diff_mu)/historic_diff_sd
      
      diff_85_array[n,j,m,] <- as.vector(mean_diff_85[1:93])
      
    }
    #printing progress
    print(paste(n,m,"done"))
  }
}

RISK_PLOTTER <- function(n){
  
  threshold <- mean(log(2)/exp(e.2$a + e.2$a_field[,n]))
  
  intercept <- mean(e.2$a + e.2$a_field[,n])
  slope <- mean(e.2$b + e.2$b_field[,n])
  
  array_1 <- exp(intercept + slope*(historical_diff_array))
  array_2 <- exp(intercept + slope*(diff_45_array))
  array_3 <- exp(intercept + slope*(diff_85_array))
  
  final_1 <- exp(-threshold*array_1)
  final_2 <- exp(-threshold*array_2)
  final_3 <- exp(-threshold*array_3)
  
  
  quant_hist <- apply(final_1,c(1,3),QUANT)
  quant_45 <- apply(final_2,c(1,3),QUANT)
  quant_85 <- apply(final_3,c(1,3),QUANT)
  
  plot(0,type="n",xlim=c(1,12),ylim=c(0.2,0.8),xlab="Month",ylab="P > M (%)",
       main=paste("Field",n),axes=FALSE)
  box()
  axis(1,at=seq(1,12,1),las=1,lab=c("J","F","M","A","M","J","J","A","S","O","N","D"))
  axis(2,at=seq(0.2,0.8,0.1),labels=seq(20,80,10))
  abline(h=0.5,lty=2)
  shade(quant_85[c(1,5),n,],lim=seq(1,12,1),col=col.alpha("red",0.75))
  #shade(quant_85[c(2,4),n,],lim=seq(1,12,1),col=col.alpha("red",0.75))
  shade(quant_45[c(1,5),n,],lim=seq(1,12,1),col=col.alpha("blue",0.75))
  #shade(quant_45[c(2,4),n,],lim=seq(1,12,1),col=col.alpha("blue",0.75))
  shade(quant_hist[c(1,5),n,],lim=seq(1,12,1),col=col.alpha("black",0.75))
  #shade(quant_hist[c(2,4),n,],lim=seq(1,12,1),col=col.alpha("black",0.75))
  #lines(y=quant_85[3,n,],x=seq(1,12,1),col="red",lwd=2)
  #lines(y=quant_45[3,n,],x=seq(1,12,1),col="blue",lwd=2)
  #lines(y=quant_hist[3,n,],x=seq(1,12),col="black",lwd=2)
  
  mtext(paste("M =",round(threshold,digits=0),"PPG"),side=3,line=-1.5)
}

par(mfrow=c(3,3))
for(n in 1:9){
  RISK_PLOTTER(n)
}
