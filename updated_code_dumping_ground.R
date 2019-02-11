#dependent packages
library(rethinking)
library(zoo)
library(RNetCDF)

#increasing max print for summaries of model parameters
options(max.print = 1000000)

#Section 1: Loading in dataframe & some preprocessing----

#loading in data frame
df = read.csv(file = "soilpop.csv", stringsAsFactors = FALSE)

#adding a reference column
df$number <- as.integer(rep( seq( 1 , 72 , 1 ), each=9 ))

#replacing 1 quadrat in field 9 with the LOD (50 F. culmorum ppg) to include field 9 in the analysis
#field 9 December 2016 had no detectable F. culmorum in all 9 quadrats sampled
df$f.c.ppg[235] = 50
#field 9 December 2017 had no detectable F. culmorum in all 9 quadrats sampled
df$f.c.ppg[559] = 50

#adding total soil Fusaria
df$tf.1 <- df$total.1*df$dilution.1
df$tf.2 <- df$total.2*df$dilution.2
df$tf.ppg <- ((df$tf.1+df$tf.2)/2)*df$solids

#Section 2: Hierarchical model treating each collection of 9 ppg estimates as equivalently random (N=72)----

#list for Stan
test.list <- list(ppg = df$f.c.ppg , field = df$field , sample = df$number)

#code for model, used as a reference
m.reference <- map2stan(
  alist(
    ppg ~ dexp(lambda),
    log(lambda) <- a_sample[sample],
    a_sample[sample] ~ dnorm(a,sigma_sample),
    a ~ dnorm(0,1),
    sigma_sample ~ dexp(1)
  ) ,
  data=test.list,
  warmup=1000,
  iter=3500,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.99)
)

#viewing chains
plot(m.reference)

#parameter summary
precis(m.reference , depth = 2)
#visualization of parameter summary
plot(precis(m.reference , depth = 2))

#extracting samples
e.1 <- extract.samples(m.reference,n=10000)

#Section 3: Quantile plot comparison: empirical vs simulated----

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
    lines(rev(seq(0.01,0.99,0.01)),-log(1-seq(0.01,0.99,0.01))/exp(e.1$a_sample[i,n]),col=col.alpha("red",0.05))
  }
  x= df$f.c.ppg[fs[n]:(fs[n]+8)]
  points((1:length(x)-0.5)/length(x),sort(x,decreasing=TRUE),type="p",pch=16,col="black")
  mtext(ns[n],line=-1.5,side=3)
  mtext(seasons[n],line=0,side=3)
}
mtext("Fraction of field exceeded", side = 1, outer = TRUE, line = 3)
mtext("PPG", side = 2, outer = TRUE, line = 3)

#field simulator

field_simulator <- function(N,M,S){
  L <- rnorm(N^2,mean=M,sd=S)
  P <- rexp(N^2,rate=exp(L))
  field <- matrix(data=P,nrow=N,ncol=N)
  field <- round(field)
  par(mar=c(0,0,0,0))
  plot(0,type="n",xlim=c(0,N),ylim=c(0,N),axes=FALSE)
  for( i in 1:N){
    for( j in 1:N){
      rect(i-0.5,j-0.5,i+0.5,j+0.5,col=col.alpha("red",alpha=(field[i,j]/max(field))/2))
      text(x=i,y=j,labels=field[i,j],cex=0.5)
    }
  }
}

field_simulator(30,-4.8,0.22)

#Section 4: Hierarchical model with offsets for field (N=9) and each collection of 9 ppg estimates (N=72)----

#model code, uses same list as above
m.field.sample <- map2stan(
  alist(
    ppg ~ dexp(lambda),
    log(lambda) <- a + a_field[field] + a_sample[sample],
    a_field[field] ~ dnorm(0,sigma_field),
    a_sample[sample] ~ dnorm(0,sigma_sample),
    a ~ dnorm(-4.8,0.22),
    sigma_field ~ dexp(1),
    sigma_sample ~ dexp(1)
  ) ,
  data = test.list ,
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
e.2 <- extract.samples(m.field.sample,n=10000)

#plotting field differences in expected value

meta_ev <- 1/exp(e.2$a[1:10000])

for(n in 1:9){
  assign(paste0("field_",n,"_ev"),1/exp(e.2$a[1:10000] + e.2$a_field[,n]))
}



plot(0,type="n",xlim=c(1,10),ylim=c(0,2000),xlab="Field",ylab="PPG",axes=FALSE)
box()
axis(1,at=seq(1,10,1),las=1,labels=c("META",1,2,3,4,5,6,7,8,9))
axis(2,at=seq(0,2000,100),las=1)
points(x=rep(1,100),y=sample(meta_ev,100),col=(col.alpha("black",0.1)),pch=15,cex=2)
points(x=rep(2,100),y=sample(field_1_ev,100),col=(col.alpha("red",0.1)),pch=15,cex=2)
points(x=rep(3,100),y=sample(field_2_ev,100),col=(col.alpha("red",0.1)),pch=15,cex=2)
points(x=rep(4,100),y=sample(field_3_ev,100),col=(col.alpha("red",0.1)),pch=15,cex=2)
points(x=rep(5,100),y=sample(field_4_ev,100),col=(col.alpha("red",0.1)),pch=15,cex=2)
points(x=rep(6,100),y=sample(field_5_ev,100),col=(col.alpha("red",0.1)),pch=15,cex=2)
points(x=rep(7,100),y=sample(field_6_ev,100),col=(col.alpha("red",0.1)),pch=15,cex=2)
points(x=rep(8,100),y=sample(field_7_ev,100),col=(col.alpha("red",0.1)),pch=15,cex=2)
points(x=rep(9,100),y=sample(field_8_ev,100),col=(col.alpha("red",0.1)),pch=15,cex=2)
points(x=rep(10,100),y=sample(field_9_ev,100),col=(col.alpha("red",0.1)),pch=15,cex=2)

#Section 5: Adding climate variables of interest----

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

df$prior_diff <- (df$p.d.90-historic_diff_mu)/historic_diff_sd

#list for model
test.list <- list(
  ppg = df$f.c.ppg,
  field = as.numeric(df$field),
  prior_diff = df$prior_diff)

m.1 <- map2stan(
  alist(
    ppg ~ dexp(lambda),
    log(lambda) <- A + B*(prior_diff) ,
    A <- a + a_field[field] ,
    B <- b + b_field[field] ,
    c(a_field,b_field)[field] ~ dmvnormNC( sigma_field , Rho ),
    a ~ dnorm(-4.80,0.22) ,
    b ~ dnorm(0,1) ,
    sigma_field ~ dcauchy(0, 1) ,
    Rho ~ dlkjcorr(7)
  ) ,
  data =  test.list ,
  warmup = 1000 ,
  iter = 3500 ,
  chains = 4 ,
  cores = 4 ,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)


precis(m.1,prob=0.95,depth=2)

#extracting samples
e.3 <- extract.samples(m.1,n=10000)

#correlating prior 90 day diff with soil moisture

#function to standardize values
standardizer <- function(x) (x-mean(x))/(sd(x))

#soil moisture reference frame
gf = read.csv(file = "soilmoist.csv", stringsAsFactors = FALSE)

gf$sm.1.s <- standardizer(gf$sm.1)
gf$sm.2.s <- standardizer(gf$sm.2)
gf$sm.3.s <- standardizer(gf$sm.3)
gf$sm.4.s <- standardizer(gf$sm.4)
gf$sm.5.s <- standardizer(gf$sm.5)
gf$sm.6.s <- standardizer(gf$sm.6)
gf$sm.7.s <- standardizer(gf$sm.7)
gf$sm.8.s <- standardizer(gf$sm.8)
gf$sm.9.s <- standardizer(gf$sm.9)

#calculateing  rolling sums
sum.10 <- rollapply(df.zoo, 10, sum, align = c("right"))
sum.20 <- rollapply(df.zoo, 20, sum, align = c("right"))
sum.30 <- rollapply(df.zoo, 30, sum, align = c("right"))
sum.40 <- rollapply(df.zoo, 40, sum, align = c("right"))
sum.50 <- rollapply(df.zoo, 50, sum, align = c("right"))
sum.60 <- rollapply(df.zoo, 60, sum, align = c("right"))
sum.70 <- rollapply(df.zoo, 70, sum, align = c("right"))
sum.80 <- rollapply(df.zoo, 80, sum, align = c("right"))
sum.90 <- rollapply(df.zoo, 90, sum, align = c("right"))

test.frame <- data.frame(df.final)
#only keeping past 2010-01-01
test.frame <- test.frame[11235:14246,]

test.frame$diff.1.sum.90.s <- standardizer(test.frame$diff.1.sum.90)
test.frame$diff.2.sum.90.s <- standardizer(test.frame$diff.2.sum.90)
test.frame$diff.3.sum.90.s <- standardizer(test.frame$diff.3.sum.90)
test.frame$diff.4.sum.90.s <- standardizer(test.frame$diff.4.sum.90)
test.frame$diff.5.sum.90.s <- standardizer(test.frame$diff.5.sum.90)
test.frame$diff.6.sum.90.s <- standardizer(test.frame$diff.6.sum.90)
test.frame$diff.7.sum.90.s <- standardizer(test.frame$diff.7.sum.90)
test.frame$diff.8.sum.90.s <- standardizer(test.frame$diff.8.sum.90)
test.frame$diff.9.sum.90.s <- standardizer(test.frame$diff.9.sum.90)

plot(type="l",gf$sm.1.s)
lines(test.frame$diff.1.sum.90.s,col="red")

ys <- c(2010,2011,2012,2013,2014,2015,2016,2017,2018)

par(oma=c(5,5,1,1))
par(mar=c(0.5,0.5,0.5,0.5))

par(mfrow=c(9,1))

plot(0,type="n",xlim=c(0,3012),ylim=c(-3,3),xlab="Year",ylab="Standard deviation",main=NULL,axes=FALSE)
box()
axis(1,at=seq(0,3012,365),las=1,labels=rep(" ",9))
axis(2,at=seq(-3,3,1),las=2)
lines(gf$sm.1.s,col="black",lwd=1)
lines(test.frame$diff.1.sum.10.s,col="red")
text(x=0,y=2.5,labels = "Field 1",cex=1.5)
abline(h=0)

plot(0,type="n",xlim=c(0,3012),ylim=c(-3,3),xlab="Year",ylab="Standard deviation",main=NULL,axes=FALSE)
box()
axis(1,at=seq(0,3012,365),las=1,labels=rep(" ",9))
axis(2,at=seq(-3,3,1),las=2)
lines(gf$sm.2.s,col="black",lwd=1)
lines(test.frame$diff.2.sum.10.s,col="red")
text(x=0,y=2.5,labels = "Field 2",cex=1.5)
abline(h=0)

plot(0,type="n",xlim=c(0,3012),ylim=c(-3,3),xlab="Year",ylab="Standard deviatinn",main=NULL,axes=FALSE)
box()
axis(1,at=seq(0,3012,365),las=1,labels=rep(" ",9))
axis(2,at=seq(-3,3,1),las=2)
lines(gf$sm.3.s,col="black",lwd=1)
lines(test.frame$diff.3.sum.10.s,col="red")
text(x=0,y=2.5,labels = "Field 3",cex=1.5)
abline(h=0)

plot(0,type="n",xlim=c(0,3012),ylim=c(-3,3),xlab="Year",ylab="Standard deviation",main=NULL,axes=FALSE)
box()
axis(1,at=seq(0,3012,365),las=1,labels=rep(" ",9))
axis(2,at=seq(-3,3,1),las=2)
lines(gf$sm.4.s,col="black",lwd=1)
lines(test.frame$diff.4.sum.10.s,col="red")
text(x=0,y=2.5,labels = "Field 4",cex=1.5)
abline(h=0)

plot(0,type="n",xlim=c(0,3012),ylim=c(-3,3),xlab="Year",ylab="Standard deviation",main=NULL,axes=FALSE)
box()
axis(1,at=seq(0,3012,365),las=1,labels=rep(" ",9))
axis(2,at=seq(-3,3,1),las=2)
lines(gf$sm.5.s,col="black",lwd=1)
lines(test.frame$diff.5.sum.10.s,col="red")
text(x=0,y=2.5,labels = "Field 5",cex=1.5)
abline(h=0)

plot(0,type="n",xlim=c(0,3012),ylim=c(-3,3),xlab="Year",ylab="Standard deviation",main=NULL,axes=FALSE)
box()
axis(1,at=seq(0,3012,365),las=1,labels=rep(" ",9))
axis(2,at=seq(-3,3,1),las=2)
lines(gf$sm.6.s,col="black",lwd=1)
lines(test.frame$diff.6.sum.10.s,col="red")
text(x=0,y=2.5,labels = "Field 6",cex=1.5)
abline(h=0)

plot(0,type="n",xlim=c(0,3012),ylim=c(-3,3),xlab="Year",ylab="Standard deviation",main=NULL,axes=FALSE)
box()
axis(1,at=seq(0,3012,365),las=1,labels=rep(" ",9))
axis(2,at=seq(-3,3,1),las=2)
lines(gf$sm.7.s,col="black",lwd=1)
lines(test.frame$diff.7.sum.10.s,col="red")
text(x=0,y=2.5,labels = "Field 7",cex=1.5)
abline(h=0)

plot(0,type="n",xlim=c(0,3012),ylim=c(-3,3),xlab="Year",ylab="Standard deviation",main=NULL,axes=FALSE)
box()
axis(1,at=seq(0,3012,365),las=1,labels=rep(" ",9))
axis(2,at=seq(-3,3,1),las=2)
lines(gf$sm.8.s,col="black",lwd=1)
lines(test.frame$diff.8.sum.10.s,col="red")
text(x=0,y=2.5,labels = "Field 8",cex=1.5)
abline(h=0)

plot(0,type="n",xlim=c(0,3012),ylim=c(-3,3),xlab="Year",ylab="Standard deviation",main=NULL,axes=FALSE)
box()
axis(1,at=seq(0,3012,365),las=1,labels=ys)
axis(2,at=seq(-3,3,1),las=2)
lines(gf$sm.9.s,col="black",lwd=1)
lines(test.frame$diff.9.sum.10.s,col="red")
text(x=0,y=2.5,labels = "Field 9",cex=1.5)
abline(h=0)

mtext("Year", side = 1, outer = TRUE, line = 3)
mtext("Standard deviation", side = 2, outer = TRUE, line = 3)

#downloading climate data for each field

locations <- read.csv("gps_ref.csv")

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

#adding a ETo dimension [variable = 7] to all the historical arrays
for (n in 1:N){
  #copying the array of interest
  A <- get(paste0(locations$town[n],"_hist"))
  #looping through the models
  for(j in 1:20){
    #calculating GDD by model
    A[7,j,] <- ETo(locations$elevation[n],
                    A[1,j,], A[2,j,], A[4,j,], A[5,j,], A[6,j,],
                    J=historical_doy,lat=(pi/180)*locations$lat[n])
    #placing the results back into the array for each town
    assign(paste0(locations$town[n],"_hist"),A)
    print(paste(locations$town[n],"historical simulation now has ETo for",model[j]))
  }
}

for (n in 1:N){
  #copying the array of interest
  A <- get(paste0(locations$town[n],"_45"))
  #looping through the models
  for(j in 1:20){
    #calculating GDD by model
    A[7,j,] <- ETo(locations$elevation[n],
                    A[1,j,], A[2,j,], A[4,j,], A[5,j,], A[6,j,],
                    J=future_doy,lat=(pi/180)*locations$lat[n])
    #placing the results back into the array for each town
    assign(paste0(locations$town[n],"_45"),A)
    print(paste(locations$town[n],"RCP 4.5 now has ETo for",model[j]))
  }
}

for (n in 1:N){
  #copying the array of interest
  A <- get(paste0(locations$town[n],"_85"))
  #looping through the models
  for(j in 1:20){
    #calculating GDD by model
    A[7,j,] <- ETo(locations$elevation[n],
                    A[1,j,], A[2,j,], A[4,j,], A[5,j,], A[6,j,],
                    J=future_doy,lat=(pi/180)*locations$lat[n])
    #placing the results back into the array for each town
    assign(paste0(locations$town[n],"_85"),A)
    print(paste(locations$town[n],"RCP 8.5 now has ETo for",model[j]))
  }
}


#just plotting with model mu and sd----

TREND_SIMULATOR <- function(n){
  #finding relevant pieces of model
  intercept <- sample(e.3$a+e.3$a_field[,n],100)
  slope <- sample(e.3$b+e.3$b_field[,n],100)
  #sample sequence
  range <- seq(-2,2,0.5)
  reference <- historic_diff_mu+range*historic_diff_sd
  ev_matrix <- matrix(0,9,100)
  for (i in 1:9){
    for (j in 1:100){
      ev_matrix[i,j] <- 1/exp(intercept[j] + slope[j]*range[i])
    }
  }
  plot(0,type="n",xlim=c(-2,2),ylim=c(0,2000),xlab="Prior 90 day P - PET (mm)",ylab="PPG",
       main=NULL,axes=FALSE)
  box()
  axis(1,at=seq(-2,2,0.5),las=1,labels=round(reference,digits=0))
  axis(3,at=seq(-2,2,0.5),las=1)
  axis(2,at=seq(0,2000,500),las=2)
  for(j in 1:100){
    lines(x=range,y=ev_matrix[,j],col=col.alpha("red",0.1))
  }
  mtext(paste(" Field",n),side=3,line=-2,adj=0,cex=1.5)
}

TREND_SIMULATOR(4)

par(mfrow=c(3,3))
for (n in 1:9){
  TREND_SIMULATOR(n)
}

#forecasting----

n <- 1
m <- 6

#differences in scale related to 100 random samples from posterior
CULMORUM_FORECASTER <- function(n,m,s){
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
  B <- get(paste0(locations$town[n],"_45"))
  C <- get(paste0(locations$town[n],"_85"))
  #adding ETo dimension
  for(j in 1:20){
    #calculating GDD by model
    A[7,j,] <- ETo(locations$elevation[n],
                   A[1,j,], A[2,j,], A[4,j,], A[5,j,], A[6,j,],
                   J=historical_doy,lat=(pi/180)*locations$lat[n])
    #placing the results back into the array for each town
    assign(paste0(locations$town[n],"_hist"),A)
    print(paste(locations$town[n],"historical simulation now has ETo for",model[j]))
  }
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
    model_zoo <- as.zoo(A[3,j,])
    model_sum <- rollapply(model_zoo,90,sum,align=c("right"))
    as.vector(model_sum)
    A[8,j,] <- c(rep(0,89),as.vector(model_sum))
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
    model_zoo <- as.zoo(A[7,j,])
    model_sum <- rollapply(model_zoo,90,sum,align=c("right"))
    as.vector(model_sum)
    A[9,j,] <- c(rep(0,89),as.vector(model_sum))
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
  #finding relevant pieces of model
  intercept <- sample(e.3$a+e.3$a_field[,n],100)
  slope <- sample(e.3$b+e.3$b_field[,n],100)
  #field scale
  field_mu <- mean(1/exp(intercept))
  field_sd <- sd(1/exp(intercept))
  #historical plot
  plot(0,type="n",xlim=c(0,56),ylim=c(field_mu-s*field_sd,field_mu+s*field_sd),xlab="Year",
       ylab="PPG",
       main=paste(month_lab[m],"- Historical"),axes=FALSE)
  box()
  axis(1,at=seq(0,56,10),las=2,labels=seq(1950,2006,10))
  axis(2,at=seq(field_mu-s*field_sd,field_mu+s*field_sd,length.out=9),las=2,labels=round(seq(field_mu-s*field_sd,field_mu+s*field_sd,length.out=9),digits=0))
  #adding field reference
  rect(-5,field_mu-1*field_sd,65,field_mu+1*field_sd,density=NA,col=col.alpha("black",0.15))
  rect(-5,field_mu-2*field_sd,65,field_mu+2*field_sd,density=NA,col=col.alpha("black",0.15))
  #adding model predictions
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
    model_predictions <- matrix(0,100,56)
    for (k in 1:56){
      model_predictions[,k] <- 1/exp(intercept +slope*mean_diff[k])
    }
    model_track <- apply(model_predictions,2,mean)
    lines(y=model_track,x=seq(1:56),col=col.alpha("red",1),lwd=1.5)
    abline(h=field_mu)
  }
  #RCP 4.5 plot
  plot(0,type="n",xlim=c(0,94),ylim=c(field_mu-s*field_sd,field_mu+s*field_sd),xlab="Year",
       ylab="PPG",
       main=paste(month_lab[m],"- RCP 4.5"),axes=FALSE)
  box()
  axis(1,at=seq(4,94,10),las=2,labels=seq(2010,2100,10))
  axis(2,at=seq(field_mu-s*field_sd,field_mu+s*field_sd,length.out=9),las=2,labels=round(seq(field_mu-s*field_sd,field_mu+s*field_sd,length.out=9),digits=0))
  #adding field reference
  rect(-5,field_mu-1*field_sd,105,field_mu+1*field_sd,density=NA,col=col.alpha("black",0.15))
  rect(-5,field_mu-2*field_sd,105,field_mu+2*field_sd,density=NA,col=col.alpha("black",0.15))
  #adding model predictions
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
    mean_precip <- sapply(split(fut_frame_45$P,fut_frame_45$Y),mean)
    mean_evap <- sapply(split(fut_frame_45$ETo,fut_frame_45$Y),mean)
    mean_diff <- ((mean_precip-mean_evap)-historic_diff_mu)/historic_diff_sd
    model_predictions <- matrix(0,100,94)
    for (k in 1:94){
      model_predictions[,k] <- 1/exp(intercept +slope*mean_diff[k])
    }
    model_track <- apply(model_predictions,2,mean)
    lines(y=model_track,x=seq(1:94),col=col.alpha("red",1),lwd=1.5)
    abline(h=field_mu)
  }
  #RCP 8.5 plot
  plot(0,type="n",xlim=c(0,94),ylim=c(field_mu-s*field_sd,field_mu+s*field_sd),xlab="Year",
       ylab="PPG",
       main=paste(month_lab[m],"- RCP 8.5"),axes=FALSE)
  box()
  axis(1,at=seq(4,94,10),las=2,labels=seq(2010,2100,10))
  axis(2,at=seq(field_mu-s*field_sd,field_mu+s*field_sd,length.out=9),las=2,labels=round(seq(field_mu-s*field_sd,field_mu+s*field_sd,length.out=9),digits=0))
  #adding field reference
  rect(-5,field_mu-1*field_sd,105,field_mu+1*field_sd,density=NA,col=col.alpha("black",0.15))
  rect(-5,field_mu-2*field_sd,105,field_mu+2*field_sd,density=NA,col=col.alpha("black",0.15))
  #adding model predictions
  for (j in 1:20){
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
    mean_precip <- sapply(split(fut_frame_85$P,fut_frame_85$Y),mean)
    mean_evap <- sapply(split(fut_frame_85$ETo,fut_frame_85$Y),mean)
    mean_diff <- ((mean_precip-mean_evap)-historic_diff_mu)/historic_diff_sd
    model_predictions <- matrix(0,100,94)
    for (k in 1:94){
      model_predictions[,k] <- 1/exp(intercept +slope*mean_diff[k])
    }
    model_track <- apply(model_predictions,2,mean)
    lines(y=model_track,x=seq(1:94),col=col.alpha("red",0.5),lwd=1.5)
    abline(h=field_mu)
  }
}

sd_seq <- c(4,4,8,12,4,8,14,4,4)

for (n in 1:9){
  par(mfrow=c(4,3))
  CULMORUM_FORECASTER(n,3,sd_seq[n])
  CULMORUM_FORECASTER(n,6,sd_seq[n])
  CULMORUM_FORECASTER(n,9,sd_seq[n])
  CULMORUM_FORECASTER(n,12,sd_seq[n])
}
