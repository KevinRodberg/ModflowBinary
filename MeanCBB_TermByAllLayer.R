
#=================================================================
# R:\ModflowBinary\MeanCBC_TermByLayer.R
#=================================================================
# Created by Kevin A. Rodberg - October 2019
#
# Purpose: Create matrix of mean 1 CBC term and 1 layer for all simulation 
#          stress periods from Modflow Binary data
#=================================================================

#=================================================================
# package management
#=================================================================
pkgChecker <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}
list.of.packages <-c( "rModflow","future.apply","tictoc",
                      "future","listenv","rasterVis","devtools","tcltk","tcltk2")
suppressWarnings(pkgChecker(list.of.packages))

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (!'githubinstall' %in% installed.packages()[,"Package"]){
  install.packages('githubinstall')
}
if ("rModflow" %in% new.packages) devtools::install_github("KevinRodberg/rModflow")

options(scipen = 999)
#=================================================================
# Set up multiprocessing environment variables
#=================================================================
cat('Setting up multiprocessor environment \n')
plan(multiprocess, .skip = TRUE)
data <- listenv()

#=================================================================
# Choose Pre-defined Modflow Model Parameters
#=================================================================
MFmodel.Params <- rModflow::defineMFmodel()
model <- rModflow::chooseModel()
M <- as.data.frame(MFmodel.Params[model,])

#=================================================================
# Choose Modflow Binary CBC file
#=================================================================
winA <- tktoplevel(bg="yellow")
msg = paste('Identify Binary Cell by Cell Budget file for :', model)
lbl.message <- tk2label(winA, text = msg, font = fontHeading, background="yellow")
tkgrid(lbl.message, padx = 30)
tkraise(winA)

mpath <- toString(MFmodel.Params[model,]$mpath)
cbbFile<-choose.files(default=mpath)

tkdestroy(winA)

if (length(cbbFile) == 0) {
  rModflow::exit("User cancelled CellxCell choice")
}
to.read = file(cbbFile, "rb")

#===============================================
# Scan CBC file for avaialble Budget Terms
# and model characteristics
#===============================================
CBCTermSet <- listBinHeaders(to.read)
ncols <- as.numeric(CBCTermSet[4])
nrows <- as.numeric(CBCTermSet[3])
nlays <- as.numeric(CBCTermSet[2])
CBCterms <- unlist(CBCTermSet[5])

#=================================================================
# Calculate number of stress periods in CBC file
#=================================================================
fileSz <- file.info(cbbFile)$size
TtlStrPd = fileSz / ((length(CBCterms) * ((ncols * nrows * nlays * 4) + 36)))

#=================================================================
# Define range of Stress Periods to read
#=================================================================
SP_rng <- readRange()
if (max(SP_rng) > TtlStrPd || min(SP_rng) < 1) {
  rModflow::exit('Out of Range')
}

#=================================================================
# Choose from Display List of available Budget Terms in CBC file
#=================================================================
close(to.read)

options <- rModflow::chooseBudgetTerms1Col(CBCterms)
n1=options$n1
term <- CBCterms[[n1]]

#=================================================================
# Function definition used to create mean rasters from binary 
# cell by cell budget terms for all stress periods by layer
#=================================================================
makeMeanLayer<- function(cbbFile,selectLayer,SP_rng,term){
  print (paste("Retrieving", trimws(CBCterms[[n1]])))
  term2Read <-trimws(CBCterms[[n1]])
  maxSP <- as.integer(TtlStrPd)
  tictoc::tic(paste('Read CBC for for layer ',selectLayer))
  #=================================================================
  #  call rModflow function to read from CBC binary file.
  #=================================================================
  filPtr <- file(cbbFile, "rb")
  CBCdata1 <- rModflow::readCBCbinByTerm(filPtr, term, SP_rng, selectLayer)
  tictoc::toc()
  #=================================================================
  # Reformat CBCdata1 as 3D array using col, row, StressPeriod dimensions
  # create dataframe of Budget values for this Layer 
  #=================================================================
  CBCsMatrix<- array(CBCdata1,c(M$ncols,M$nrows,maxSP))
  
  tictoc::tic(paste('Creating mean raster for layer ',selectLayer))
  xf <-future.apply::future_apply(CBCsMatrix,MARGIN=c(1,2),FUN=mean,na.rm=T)
  MeanSim <- t(xf[,])
  
  #=================================================================
  # NAD83 HARN StatePlane Florida East FIPS 0901 Feet
  #=================================================================
  HARNSP17ft  = CRS("+init=epsg:2881")
  HARNUTM17Nm  = CRS("+init=epsg:3747")
  latlongs = CRS("+proj=longlat +datum=WGS84")
  #=================================================================
  # calculate number of rows and columns
  #=================================================================
  res=MFmodel.Params[model,]$res
  xmin=MFmodel.Params[model,]$xmin
  ymin=MFmodel.Params[model,]$ymin
  rasRows=MFmodel.Params[model,]$nrows
  rasCols=MFmodel.Params[model,]$ncols
  xmax=xmin+(res*rasCols)
  ymax=ymin+(res*rasRows)
  #=================================================================
  # define raster and map extents using MFmodel data extents
  #=================================================================
  cellsize=c(res,res)
  ras <- raster::raster(res=cellsize, xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,crs=HARNSP17ft)
  rasExt <- raster::extent(ras)
  #=================================================================
  # create raster from meanSim matrix
  #=================================================================
  MeanRasL1<-raster::raster(MeanSim,xmin,xmax,ymin,ymax, crs=HARNSP17ft)
  
  #=================================================================
  # Export raster as tiff file 
  #=================================================================  
  basePath <- dirname(cbbFile)
  filename = paste0(basePath,'/Mean_',term2Read,'_Lay',selectLayer,'.tif')
  raster::writeRaster(MeanRasL1, filename, format="GTiff", overwrite=TRUE)
}

ix=0
tictoc::tic("Raster Creation")

for (selectLayer in seq(1,M$nlays)){
  #=================================================================
  # Retrieve CBC by Layer using multiple processors
  #=================================================================
  cat(paste("Start making rasters for",term,"Layer",selectLayer, '\n'))
  ix = ix + 1
  data[[ix]] %<-%  makeMeanLayer(cbbFile,selectLayer,SP_rng,term)
}
cat(paste("Rasters will be created in ",dirname(cbbFile)),"in a couple minutes \n")
lastOne <- future::futureOf(data[[ix]])
while (!future::resolved(lastOne)) {
     cat ('-')
          Sys.sleep(0.5)
      }
tictoc::toc()
