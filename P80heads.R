list.of.packages <-c( "data.table","tcltk2","rModflow","future.apply","tictoc")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (!'githubinstall' %in% installed.packages()[,"Package"]){
  install.packages('githubinstall')
}
library(devtools)
if(length(new.packages)) install.packages(new.packages)
if ("rModflow" %in% new.packages) devtools::install_github("KevinRodberg/rModflow")

library (data.table)
library(tcltk2)
library(rModflow)
library(future.apply)
library(rasterVis)
library(RColorBrewer)
library(tictoc)
#=================================================================
# R:\ModflowBinary\P80heads.R
#=================================================================
# Beginning of P80 Modflow Heads
#
# Created by Kevin A. Rodberg - October 2018
#
# Purpose: Create matrix of Layer 1 heads for all simulation 
#          stress periods from Modflow Binary data
#          and calculate P80

#=================================================================
# Choose Modflow Binary Heads file
#=================================================================
source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ReusableFunctions/tclFuncs.R")

readHeadsbinByLay <- function(filPtr, selectLayer,maxSP) {
  bigVector <- vector('numeric')
  HeaderRead <- readHeadsHeader(filPtr)
  kntFloats <- HeaderRead$K * HeaderRead$NR * HeaderRead$NC
  Lay1floats <- HeaderRead$NR * HeaderRead$NC
  HeadBlock <- readBin(filPtr, double(), n = Lay1floats, size = 4)
  bigVector <- c(bigVector, HeadBlock[1:Lay1floats])
  i <- 1
  cat(paste("0%.."))
  SP_rng <- maxSP-1
  repeat {
    HeaderRead <- readHeadsHeader(filPtr)
    # Don't read past EOF
    if (length(HeaderRead) > 0) {
      if (HeaderRead$K == selectLayer) {
        i <- i + 1
        HeadBlock <-
          readBin(filPtr, double(), n = Lay1floats, size = 4)
        bigVector <- c(bigVector, HeadBlock[1:Lay1floats])
      } else {
        seek(filPtr, (Lay1floats * 4), origin = 'current')
      }
    }
    # don't read everything unless necessary
    if (length(HeaderRead) == 0) {
      cat('\n')
      break
    }

    if (HeaderRead$KPER > max(SP_rng)) {
      cat('\n')
      break
    }
    # Display % complete
    cat(paste('\r',format(as.numeric(HeaderRead$KPER) / max(SP_rng) * 100,digits = 2,nsmall = 2),"%"))
  }

  return(bigVector)
}

MFmodel.Params <- defineMFmodel()
model <- chooseModel()
M <- as.data.frame(MFmodel.Params[model,])


winA <- tktoplevel()
msg = paste('Identify Binary Heads file for :', model)
lbl.message <- tk2label(winA, text = msg, font = fontHeading)
tkgrid(lbl.message, padx = 30)
tkraise(winA)

mpath <- toString(MFmodel.Params[model,]$mpath)
headsFile<-choose.files(default=mpath)

tkdestroy(winA)

if (length(headsFile) == 0) {
  exit("User cancelled HeadsFile choice")
}
to.read = file(headsFile, "rb")

#===============================================
# Estimate number of stress periods in Heads file
#===============================================
fileSz <- file.info(headsFile)$size
TtlStrPd = fileSz / ( M$nlays * ((M$ncols * M$nrows * 4) + 44))

#===============================================
# Define range of Stress Periods to read
#===============================================
SP_rng <- readRange()
if (max(SP_rng) > TtlStrPd || min(SP_rng) < 1) {
  exit('Out of Range')
}

#===============================================
# Retrieve Heads by Layer
#===============================================
to.read <- file(headsFile, "rb")
selectLayer = 1
maxSP <- as.integer(TtlStrPd)
Layer1 <- readHeadsbinByLay(to.read, selectLayer, maxSP)
close(to.read)

yourTheme = rasterTheme(region = brewer.pal('BrBG', n = 9))
# Reformat Layer1 as 3D array using col, row, StressPeriod dimensions
# create dataframe of Head values for this Layer (layer 1)
HeadsMatrix<- array(Layer1,c(M$ncols,M$nrows,maxSP))
plan(multiprocess)
#system.time(x<-apply(HeadsMatrix,c(1,2),quantile,probs=c(.5),na.rm=T))
tic()
xf <-future_apply (HeadsMatrix,MARGIN=c(1,2),FUN=stats::quantile,probs=c(.2),na.rm=T)
toc()
my.at = seq(-1,150,15)
levelplot(raster(t(xf[,])),par.settings = yourTheme,at=my.at)
