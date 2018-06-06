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
library(rModflow)
options(warn=-1)
options(scipen = 999)

source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ReusableFunctions/rasterFuncs.R")
source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ReusableFunctions/gisFuncs.R")
source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ReusableFunctions/tclFuncs.R")



#=================================================================

# Beginning of Script

MFmodel.Params <- defineMFmodel()
if (!exists("MFmodel")) {
  MFmodel <- 'x'
}

if (!(MFmodel %in% (c('ECFTX', 'NPALM', 'LWCSIM')))) {
  MFmodel <- chooseModel()
  if ((!MFmodel %in% (c('ECFTX', 'NPALM', 'LWCSIM')))) {
    exit('Incorrect model choice')
  }
}

res     = MFmodel.Params[MFmodel,]$res
xmin    = MFmodel.Params[MFmodel,]$xmin
ymin    = MFmodel.Params[MFmodel,]$ymin
nrows   = MFmodel.Params[MFmodel,]$nrows
ncols   = MFmodel.Params[MFmodel,]$ncols
nsp     = MFmodel.Params[MFmodel,]$nsp
startYr = MFmodel.Params[MFmodel,]$startYr
freq    = MFmodel.Params[MFmodel,]$freq
if (!exists("SP_rng")) {  SP_rng <- seq(1, nsp)}
SPknt <- length(SP_rng)

cat(paste("Plot functions will provide output for the following stress periods: \n\t",list(SP_rng)),"\n")

# Calculate raster extents

xmax = xmin + (ncols * res)
ymax = ymin + (nrows * res)
modelras <-  raster(resolution = res,xmn = xmin,xmx = xmax,ymn = ymin,ymx = ymax,
                    crs = HARNSP17ft  )

WMDbnd.Path <- "//whqhpc01p/hpcc_shared/krodberg/NexRadTS"
WMDbnd.Shape <- "CntyBnds.shp"
setwd(WMDbnd.Path)
WMDbnd <- readShapePoly(WMDbnd.Shape, proj4string = HARNSP17ft)
rasPoly <- as(extent(modelras), 'SpatialPolygons')
proj4string(rasPoly) <- HARNSP17ft
clpBnds2 <- gClip(WMDbnd, rasPoly)

# Eliminate Scientific Notation
options(scipen = 10000)
myTheme = rasterTheme(region = brewer.pal('Oranges', n = 9))
yourTheme = rasterTheme(region = brewer.pal('Blues', n = 9))
newTheme = rasterTheme(region = brewer.pal('RdBu', n = 9))

DaysPerMonth <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


pltOpts<-definePlotOpts()

matrixSource=pltOpts$matrixSource
printingOn=pltOpts$printingOn
printAnnualOn=pltOpts$printAnnualOn
Animform=pltOpts$Animform

datSrc <- chooseDataSource(matrixSource)

# create points from each column and convert
# them to rasters and add each to the raster stack

cat (paste("Raster files will be saved by Stress Period is", printingOn, "\n"))

#==============================
# Clear raster Stack
#==============================
rasStack <- stack()
MFrasType <- datSrc$rasType
printingOn = TRUE
rasNames <- sprintf("%s  %s", MFrasType,format(seq(as.Date(paste0(startYr, "/1/1")),
            by = "month",length.out = max(SP_rng)),format = "%b %Y"))
indices <-  as.numeric(sprintf("%s",format(seq(as.Date(paste0(startYr, "/1/1")),
            by = "month",length.out = max(SP_rng)),format = "%Y")))

rasNamesSubset <- vector()
idxList <-vector()
#============================================================
# Display raster creation progress as % complete
#============================================================
if (printingOn) {  cat(paste0("Raster plots are being created: \n",
      "//ad.sfwmd.gov/dfsroot/data/wsd/MOD/",MFmodel,"/NetRecharge/plots \n"))
} else{
  cat(paste0("Rasters are being processed\n"))
}
cat(paste("0%.."))
rampScale <-  equalIntSeries(min(datSrc$Array3D) * 12 / (res * res), max(datSrc$Array3D) * 12 /(res * res))
for (iter in seq(1:length(SP_rng)))  {
  StressPeriod = SP_rng[iter]
  
  year <- startYr + floor(StressPeriod / 12)
  if (StressPeriod %% 12 > 0) {
    imonth <- StressPeriod %% 12
  } else {
    imonth <- 12
    year <- year - 1
  }

#============================================================
# swap rows for columns with transpose function: t()
# raster from array converted to inches per day from cuft per Day
#============================================================
  Array2dBySP <- t(datSrc$Array3D[, , iter] * 12 / (res * res))
  ras <- raster(Array2dBySP)
  crs(ras) <- HARNSP17ft
  e <- extent(c(xmin, xmax, ymin, ymax))
  extent(ras) <- e
  
  #============================================================
  #  plot rasters if printing Flag is set to TRUE
  #============================================================
  if (printingOn) {
    # use layout option to plot c(cols,rows) figures per page
    rchFile = paste("//ad.sfwmd.gov/dfsroot/data/wsd/MOD/",MFmodel,
                    "/NetRecharge/plots/",MFrasType,year,"_",imonth,".png",sep = "")
    png(file = rchFile,width = 3000,height = 2400,units = "px",res = 300)
    print(levelplot(ras,contour = FALSE,par.settings = yourTheme,at=rampScale,scales = list(x = list(rot = 90))) 
          + layer(sp.polygons(WMDbnd))
    )
    dev.off()
  }
  rasStack <- stack(rasStack, ras)
  cat(paste('\r', format(as.numeric(iter) / length(SP_rng) * 100,digits = 2,nsmall = 2), "%"))
  rasNamesSubset[iter] <- rasNames[StressPeriod]
  idxList[iter]<- year-startYr + 1
}
rm(datSrc)
gc(verbose=TRUE)

#============================================================
#  Assign names to each raster based on Stress Period
#============================================================
if (printAnnualOn) {
  cat(
    paste0("\nAnnual Raster plots are being created: \n",
      "//ad.sfwmd.gov/dfsroot/data/wsd/MOD/",MFmodel,"/NetRecharge/plots \n")
  )
} else{
  cat(paste0("\nAnnual Rasters are being processed\n"))
}
cat(paste("0%.."))
names(rasStack) <- rasNamesSubset
s.sum <- stackApply(rasStack, indices=idxList, fun = sum, na.rm=TRUE)
rasMin <- min(slot(slot(s.sum, "data"), 'min'))
rasMax <- max(slot(slot(s.sum, "data"), 'max'))
rampScale <-  equalIntSeries(rasMin,rasMax)

if (printAnnualOn) {
  # use layout option to plot c(cols,rows) figures per page
  idx = 0
  for (year in seq(min(idxList), max(idxList))) {
    idx = idx + 1
    {
      rchFile = paste("//ad.sfwmd.gov/dfsroot/data/wsd/MOD/",
        MFmodel,"/NetRecharge/plots/",MFrasType,year+startYr-1,"_total.png",sep = "")
      png(file = rchFile,width = 3000,height = 2400,units = "px",res = 300)
      
      print(levelplot(s.sum[[idx]],contour = FALSE,par.settings = newTheme,at=rampScale,
          main = paste(MFrasType, '_', slot(s.sum[[idx]]@data, 'names')),
          scales = list(x = list(rot = 90))) + layer(sp.polygons(clpBnds2)))
      dev.off()
    }
    cat(paste('\r', format((year-min(idxList)+1) / (max(idxList)-min(idxList)+1) * 100,digits = 2,nsmall = 2), "%"))
  }
  cat('\n')
}

#============================================================
# Animations require ImageMagick conversion tool installed on
# system R is running on.

ani.options(convert = 'h:\\magick.exe')
#ani.options(ffmpeg  = 'h:\\ffmpeg.exe')
ani.options(ffmpeg  = '//ad.sfwmd.gov/dfsroot/data/wsd/sup/devel/source/ReusableFunctions/ffmpeg.exe')
#ani.options(ffmpeg = 'C:/Program Files/ImageMagick-7.0.7-Q16/ffmpeg.exe')

#============================================================

#loop = 0 will enable auto replay, .3 is the delay interval
ani.options(
  nmax = SPknt,
  loop = FALSE,
  ani.width = 1000,
  ani.height = 700,
  interval = .3
)
LegendBreaks = seq(0, 10, .5)
Orangecols <- brewer.pal(9, 'Oranges')
bluecols <- brewer.pal(9, 'Blues')
setwd('h:\\')
timeSec <- 5 * SPknt
#============================================================
# Create animation as if flag set to GIF or AVI
#============================================================
if (Animform == 'off'){cat(paste("Animations will not be saved \n"))}
if (Animform == 'GIF') {
  newcol <- colorRampPalette(bluecols)
  ramp <- newcol(100)
  RCHgif <-paste0("//ad.sfwmd.gov/dfsroot/data/wsd/MOD/",
      MFmodel,"/NetRecharge/plots/netRCH.gif")
  cat("Creating GIF amimation... this may take ",timeSec / 60," minutes, ...\n\n")
  saveGIF(animate(rasStack, addfun = DistrictBnds, col = ramp),
          movie.name = RCHgif,convert = "convert"  )
} else if (Animform == 'AVI') {
  newcol <- colorRampPalette(bluecols)
  ramp <- newcol(100)
  RCHavi <-
    paste0("//ad.sfwmd.gov/dfsroot/data/wsd/MOD/",MFmodel,"/NetRecharge/netRCH.avi")
  cat("\n\n\nCreating AVI amimation... this may take ",timeSec / 60," minutes,...\n\n")
  saveVideo(animate(rasStack, addfun = DistrictBnds, col = ramp),
            video.name = RCHavi,convert = "ffmpeg"
  )
}
options(warn=0)

cat ("open //ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ModflowBinary/netRCH_anim.R \n to continue raster plotting")
