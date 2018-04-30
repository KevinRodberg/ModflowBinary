library(reshape2)
library(raster)
library(dplyr)
library(rasterVis)
library(data.table)
library(sp)
library(rgdal)
library(maps)
library(maptools)
library(animation)
library(rgeos)

geomSeries <- function(base, max) {
  (base^(0:floor(log(max, base)))-1)/100
}

clipToExtent <- function( sp, extent ) {
  require(rgeos)
  keep <- gContains( extent, sp,byid=TRUE ) | gOverlaps( extent, sp,byid=TRUE )
  stopifnot( ncol(keep)==1 )
  sp[drop(keep),]
}
gClip <- function(shp, bb){
  if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else b_poly <- as(extent(bb), "SpatialPolygons")
  print (shp)
  print (b_poly)
  gIntersection(shp, b_poly, byid = T)
}
DistrictBnds <- function() {
#  color=rgb(0,0,0,alpha=.01)
#  plot(clpBnds2,col=color )
  plot(clpBnds2, bg="transparent", add=TRUE)
}
# Clear raster Stack
rasStack <-stack()

# Clear raster Stack
# Calculate raster extents
UTM17m = CRS("+init=epsg:26917") # UTM Zone 17 Meters
llwgs84  = CRS("+init=epsg:4326") # lat-long projection
HARNSP17ft  = CRS("+init=epsg:2881") # NAD_1983_HARN_StatePlane_Florida_East_FIPS_0901_Feet
#  Meters 381 = 1250 Feet

res = 1250
xmin =24352.000
ymin =983097.000

nrows = 603  
ncols = 740
xmax = xmin+ (ncols*res)
ymax = ymin + (nrows*res)
ECFTXras <- raster(resolution=1250, xmn=xmin, xmx=xmax,ymn=ymin, ymx=ymax,crs=HARNSP17ft)

WMDbnd.Path<-"//whqhpc01p/hpcc_shared/krodberg/NexRadTS"
WMDbnd.Shape<-"CntyBnds.shp"
setwd(WMDbnd.Path)
WMDbnd <- readShapePoly(WMDbnd.Shape,proj4string=HARNSP17ft)

rasPoly <- as(extent(ECFTXras),'SpatialPolygons')
proj4string(rasPoly) <- HARNSP17ft
clpdBnds<-clipToExtent(WMDbnd, rasPoly)
clpBnds2<-gClip(WMDbnd, rasPoly)

# Eliminate Scientific Notation
options(scipen=10000)
myTheme=rasterTheme(region=brewer.pal('Oranges', n=9))
yourTheme=rasterTheme(region=brewer.pal('Blues', n=9))

DaysPerMonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
#rchMonthly <-DiffVectorBySP
rchMonthly <-ETVectorBySP

SPknt = dim(DiffVectorBySP)[4]
# create points from each column and convert 
# them to rasters and add each to the raster stack
#====================================
printingOn <- TRUE
#printingOn <- FALSE
#====================================
RCHrasStack <- stack()
for (StressPeriod in seq(1:SPknt))
#for (StressPeriod in seq(1:12))
  {
  # swap rows for columns with transpose function t()
  year <- 1996 + floor(StressPeriod/12)
  if (StressPeriod%%12>0){
    imonth <- StressPeriod%%12
#    rchMon_1 <- t(rchMonthly[,,1,StressPeriod])*as.numeric(12*DaysPerMonth[imonth])
    rchMon_1 <- t(rchMonthly[,,1,StressPeriod])
    print(paste(StressPeriod, " ", mean(rchMon_1)))
    
  } else {
    imonth <- 12
    year <- year-1
#    rchMon_1 <- t(rchMonthly[,,1,StressPeriod])*as.numeric(12*DaysPerMonth[imonth])
    rchMon_1 <- t(rchMonthly[,,1,StressPeriod])
    print(paste(StressPeriod, " ", mean(rchMon_1)))
  }
  RCHras <-raster(rchMon_1)
  e<-extent(c(xmin,xmax,ymin,ymax))
  extent(RCHras) <-e
  if (printingOn){
    # layout provides options to plot c(cols,rows) figures per page
    rchFile=paste("//ad.sfwmd.gov/dfsroot/data/wsd/MOD/ECFTX/NetRecharge/plots/ET",year,"_",imonth,".png",sep="")
    png( file=rchFile, width=3000,height=2400,units="px",res=300)
    print(levelplot(RCHras,contour=FALSE,par.settings = yourTheme, scales=list(x=list(rot=90))) + layer(sp.polygons(WMDbnd)))
    dev.off()  
  }
  RCHrasStack <- stack(RCHrasStack, RCHras)
}

exit('Stop')
rasNames <- sprintf("RCH (inches Per Month) %s",format(seq(as.Date("1996/1/1"), by = "month", length.out = SPknt),format = "%b %Y"))
names(RCHrasStack)<-rasNames
#==============================
Animform<-'AVI'
# AnimForm<-'GIF'
#==============================
#Aimations require ImageMagick conversion tool
ani.options(convert = 'h:\\magick.exe')
ani.options(ffmpeg  = 'h:\\ffmpeg.exe')
#ani.options(ffmpeg = 'C:/Program Files/ImageMagick-7.0.7-Q16/ffmpeg.exe')
#loop = 0 will enable auto replay, .3 is the delay interval
ani.options(nmax = 30, loop=1, ani.width =1000, ani.height = 700, interval = .3)
LegendBreaks=seq(0,10,.5)
Orangecols <- brewer.pal(9, 'Oranges')
bluecols <- brewer.pal(9, 'Blues')
setwd('h:\\')
if (Animform=='GIF'){
  newcol <- colorRampPalette(bluecols)
  ramp <-newcol(100)
  RCHgif <- "//ad.sfwmd.gov/dfsroot/data/wsd/MOD/ECFTX/NetRecharge/netRCH.gif"
  saveGIF(animate(RCHrasStack, addfun=DistrictBnds, col=ramp ), movie.name = RCHgif,convert="convert")
} else if (Animform =='AVI') {
  newcol <- colorRampPalette(Orangecols)
  ramp <-newcol(100)
  RCHavi <- "//ad.sfwmd.gov/dfsroot/data/wsd/MOD/ECFTX/NetRecharge/netRCH.avi"
  saveVideo(animate(RCHrasStack, addfun=DistrictBnds, col=ramp ),video.name=RCHavi, convert="ffmpeg" )
}
#print(levelplot(ETrasStack,contour=FALSE,par.settings = myTheme, layout=c(4,1),scales=list(x=list(rot=90))) + layer(sp.polygons(WMDbnd)))

