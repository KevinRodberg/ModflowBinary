pkgChecker <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}

list.of.packages <-c( "data.table","tcltk2","rModflow","future.apply","tictoc",
                      "future","listenv","rasterVis","RColorBrewer")
suppressWarnings(pkgChecker(list.of.packages))

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if (!'githubinstall' %in% installed.packages()[,"Package"]){
  install.packages('githubinstall')
}
library(devtools)
if(length(new.packages)) install.packages(new.packages)
if ("rModflow" %in% new.packages) devtools::install_github("KevinRodberg/rModflow")

#=================================================================
# R:\ModflowBinary\MeanCBC_TermByLayer.R
#=================================================================
# Beginning of P50 Modflow CBC
#
# Created by Kevin A. Rodberg - October 2018
#
# Purpose: Create matrix of mean 1 CBC term and 1 layer for all simulation 
#          stress periods from Modflow Binary data

#=================================================================
# Choose Modflow Binary CBC file
#=================================================================
source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ReusableFunctions/tclFuncs.R")

#' @title Read budget term by layer
#' @description \code{readCBCbinByTerm} search for a defined budget term
#'      and returns a vector of values by stress periods identified in range of values
#' @param filPtr, term, SP_rng, lay
#' @return bigVector
#' @export

readCBCbinByTerm <- function(filPtr, term, SP_rng, lay) {
  bigVector <- vector('numeric')
  HeaderRead <- readCBCHeader(filPtr)
  kntFloats <- HeaderRead$K * HeaderRead$NR * HeaderRead$NC
  Lay1floats <- HeaderRead$NR * HeaderRead$NC
  cbcBlock <- readBin(filPtr, double(), n = kntFloats, size = 4)
  strt<-1+((lay-1)*Lay1floats)
  end <- lay*Lay1floats
  i <- 1
  cat(paste("0%.."))
  if (HeaderRead$TEXT == term  &&
      is.element(HeaderRead$KPER, SP_rng &&
                 thisHeader$KPER <= max(SP_rng))) {
    bigVector <- c(bigVector , cbcBlock)
    i <- i + 1
  } else{
    repeat {
      thisHeader <- readCBCHeader(filPtr)
      # Don't read past EOF
      if (length(thisHeader) > 0) {
        if (term ==thisHeader$TEXT   &&
            is.element(thisHeader$KPER, SP_rng) &&
            thisHeader$KPER <= max(SP_rng)) {
          i <- i + 1
          cbcBlock <-
            readBin(filPtr, double(), n = kntFloats, size = 4)
          
          bigVector <- c(bigVector, cbcBlock[strt:end])
          # bigVector <- c(bigVector, cbcBlock[1:Lay1floats])
        } else {
          seek(filPtr, (kntFloats * 4), origin = 'current')
        }
      }
      # don't read everything unless necessary
      if (length(thisHeader) == 0) {
        cat('\n')
        break
      }
      
      if (thisHeader$KPER > max(SP_rng)) {
        cat('\n')
        break
      }
      # Display % complete
      cat(paste('\r',format(as.numeric(thisHeader$KPER) / max(SP_rng) * 100,digits = 2,nsmall = 2),"%"))
    }
  }
  return(bigVector)
}
#' @title Choose Modflow Cell-by-cell Budget Terms
#' @description GUI choices for Modflow Budget Terms
#' @param CBCterms vector created by \code{listBinHeaders}
#' @export
chooseBudgetTerms1Col <- function(CBCterms) {
  win2 <- tktoplevel()
  frame1 <-  tk2frame(win2,borderwidth = 3,relief = "sunken",padding = 10)
  frame3 <-  tk2frame(win2,borderwidth = 3,relief = "sunken",padding = 10)
  lbl.CBCSelect <- tk2label(win2, text = "Select Budget Terms to Extract from CBC file", font = fontHeading)
  tkpack(lbl.CBCSelect,  side = "top",  expand = FALSE,  ipadx = 5,  ipady = 5,  fill = "x")
  tkpack(frame3,side = "bottom",expand = TRUE,fill = "both")
  tkpack(frame1,side = "left",expand = TRUE,fill = "both")
  rBtnVal1 = tclVar(trimws(CBCterms[[1]]))
  
  # create 1 columns of CBCterms radioButtons
  btns.f1 = vector()
  for (num in seq(1, length(CBCterms))) {
    btn <- tk2radiobutton(frame1)
    tkconfigure(btn, variable = rBtnVal1, value = num)
    tkgrid(tk2label(frame1, text = trimws(CBCterms[[num]])),btn,padx = 10,pady =  5)
    btns.f1 = append(btns.f1, btn)
  }
  
  tkgrid(tk2button(frame3,text ="Cancel",width = -6,command = fnCncl),
         tk2button(frame3,text ="OK",    width = -6,command = fnOK  ), padx = 10,pady = c(5, 15))
  tkbind(win2, "<Return>", fnOK)
  
  tkraise(win2)
  tkwait.variable(done)
  tkdestroy(win2)
  if (tclvalue(done) != 1) {
    exit("User canceled Model Selection")
  }
  n1  <- as.integer((tclvalue(rBtnVal1)))
  return(list(n1=n1))
}

#=================================================================
# Choose Modflow Cell by Cell Budget Term file
#=================================================================
MFmodel.Params <- defineMFmodel()
model <- chooseModel()
M <- as.data.frame(MFmodel.Params[model,])

winA <- tktoplevel()
msg = paste('Identify Binary Cell by Cell Budget file for :', model)
lbl.message <- tk2label(winA, text = msg, font = fontHeading)
tkgrid(lbl.message, padx = 30)
tkraise(winA)

mpath <- toString(MFmodel.Params[model,]$mpath)
cbbFile<-choose.files(default=mpath)

tkdestroy(winA)

if (length(cbbFile) == 0) {
  exit("User cancelled CellxCell choice")
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

#===============================================
# Calculate number of stress periods in CBC file
#===============================================
fileSz <- file.info(cbbFile)$size
TtlStrPd = fileSz / ((length(CBCterms) * ((ncols * nrows * nlays * 4) + 36)))

#===============================================
# Define range of Stress Periods to read
#===============================================
SP_rng <- readRange()
if (max(SP_rng) > TtlStrPd || min(SP_rng) < 1) {
  exit('Out of Range')
}

#===============================================
# Display List of available Budget Terms in CBC file
#===============================================
close(to.read)

options <- chooseBudgetTerms1Col(CBCterms)
n1=options$n1
term <- CBCterms[[n1]]

makeMeanLayer<- function(cbbFile,selectLayer,SP_rng,term){
  to.read <- file(cbbFile, "rb")
  filPtr<- to.read 
  print (paste("Retrieving", trimws(CBCterms[[n1]])))
  term2Read <-trimws(CBCterms[[n1]])
  
  maxSP <- as.integer(TtlStrPd)
  tic(paste('Reading CBC for for layer ',selectLayer))
  CBCdata1 <- readCBCbinByTerm(to.read, term, SP_rng, selectLayer)
  
  close(to.read)
  toc()
  yourTheme = rasterTheme(region = brewer.pal('BrBG', n = 9))
  # Reformat CBCdata1 as 3D array using col, row, StressPeriod dimensions
  # create dataframe of Head values for this Layer 
  CBCsMatrix<- array(CBCdata1,c(M$ncols,M$nrows,maxSP))
  
  tic(paste('Creating mean raster for layer ',selectLayer))
  #xf <-future_apply (CBCsMatrix,MARGIN=c(1,2),FUN=stats::quantile,probs=c(.1),na.rm=T)
  xf <-future_apply (CBCsMatrix,MARGIN=c(1,2),FUN=mean,na.rm=T)
  # my.at = seq(-80,100,15)
  # levelplot(raster(t(xf[,])),par.settings = yourTheme,at=my.at)
  
  MeanSim <- t(xf[,])
  
  #--------------------------------------------------------------------------------------------------
  # NAD83 HARN StatePlane Florida East FIPS 0901 Feet
  #--------------------------------------------------------------------------------------------------
  HARNSP17ft  = CRS("+init=epsg:2881")
  HARNUTM17Nm  = CRS("+init=epsg:3747")
  latlongs = CRS("+proj=longlat +datum=WGS84")
  #--------------------------------------------------------------------------------------------------
  # calculate number of rows and columns
  #--------------------------------------------------------------------------------------------------
  res=MFmodel.Params[model,]$res
  xmin=MFmodel.Params[model,]$xmin
  ymin=MFmodel.Params[model,]$ymin
  rasRows=MFmodel.Params[model,]$nrows
  rasCols=MFmodel.Params[model,]$ncols
  xmax=xmin+(res*rasCols)
  ymax=ymin+(res*rasRows)
  
  cellsize=c(res,res)
  ras <- raster::raster(res=cellsize, xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,crs=HARNSP17ft)
  #--------------------------------------------------------------------------------------------------
  # define raster and map extents using MFmodel data extents
  #--------------------------------------------------------------------------------------------------
  rasExt <- raster::extent(ras)
  
  MeanRasL1<-raster::raster(MeanSim,xmin,xmax,ymin,ymax, crs=HARNSP17ft)
  my.at = c(minValue(MeanRasL1),-1,1,maxValue(MeanRasL1))
  basePath <- dirname(cbbFile)
  levelplot(MeanRasL1,par.settings = yourTheme,at=my.at)
  Sys.sleep(0)
  filename = paste0(basePath,'/Mean_',term2Read,'_Lay',selectLayer,'.tif')
  raster::writeRaster(MeanRasL1, filename, format="GTiff", overwrite=TRUE)
}
plan(multiprocess)

data <- listenv()
ix=0
tic("Raster Creation")
for (selectLayer in seq(1,M$nlays)){
  #===============================================
  # Retrieve CBC by Layer
  #===============================================
  #----------------------------------------------------------------------------
  cat(paste("Start making rasters for",term,"Layer",selectLayer, '\n'))
      ix = ix + 1
      data[[ix]] %<-%  makeMeanLayer(cbbFile,selectLayer,SP_rng,term)
}
cat(paste("Rasters will be created in ",dirname(cbbFile)),"in a couple minutes \n")
lastOne <- futureOf(data[[ix]])
while (!resolved(lastOne)) {
     cat ('-')
          Sys.sleep(0.5)
      }
toc()
