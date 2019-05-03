#==================================================================================================
# R:\ModflowBinary\P80headDifference.R
#==================================================================================================
#==================================================================================================
# Beginning of P80 Modflow Heads
#
# Created by Kevin A. Rodberg - February 2019
#
# Purpose: Create difference matrix of p80 Reference Condition Heads  
#          minus another P80 Simulation Heads from Layer 1 and  
#          the specified stress periods (POR) from Modflow runs using
# makes use of by R tools for Modflow:
#      rModflow::readHeadsbinByLay 
#         from install_github("KevinRodberg/rModflow")
#==================================================================================================
# source: 
# //ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ModflowBinary/P80heads.R
#==================================================================================================

list.of.packages <-c( "data.table","devtools","utils",
                      "tcltk2","rModflow",
                      "future.apply","future","listenv",
                      
                      "rasterVis","sp","maptools","rgeos","raster",
                      "ggplot2","RColorBrewer",
                      "tictoc",
                      "polynom")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (!'githubinstall' %in% installed.packages()[,"Package"]){
  install.packages('githubinstall')
}
if(length(new.packages)) install.packages(new.packages)
library(devtools)
if ("rModflow" %in% new.packages) devtools::install_github("KevinRodberg/rModflow")
lapply(list.of.packages,require, character.only=TRUE)

source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ReusableFunctions/tclFuncs.R")
skip <- FALSE 
if(exists('RCheadsFile') & exists('SIMheadsFile') ){
  if(file.exists(RCheadsFile) & file.exists(SIMheadsFile) ){
    if (utils::askYesNo(paste("Do you want to use the same binary heads selections?\n",
                               RCheadsFile,'\n',SIMheadsFile,'\n'),
                        prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))){
      cat('Bypassing data selections \n')
      skip <- TRUE
    } else { 
      skip <- FALSE 
    }
  } else {
    skip <- FALSE
  }
} 
if (!skip){
  
  #=================================================================
  # Choose Modflow Model to be processed via GUI
  # such as ECFTX, NPALM, LWCSIM, etc
  #=================================================================
  MFmodel.Params <- defineMFmodel()
  model <- chooseModel()
  M <- as.data.frame(MFmodel.Params[model,])
  
  #=================================================================
  # Select first Modflow Binary Heads file to process
  #=================================================================
  winA <- tktoplevel()
  msg = paste('Identify Binary Heads file for :', model)
  lbl.message <- tk2label(winA, text = msg, font = fontHeading)
  tkgrid(lbl.message, padx = 30)
  tkraise(winA)
  MFmodel.Params[model,]$mpath <- 
    '\\\\ad.sfwmd.gov\\dfsroot\\data\\wsd\\SUP\\proj\\CFWI_WetlandStress\\Update2018\\ModelRuns\\*.*'
  mpath <- toString(MFmodel.Params[model,]$mpath)
  RCheadsFile<-choose.files(default=mpath)
  tkdestroy(winA)
  
  if (length(RCheadsFile) == 0) {
    exit("User cancelled HeadsFile choice")
  }
  #=================================================================
  # Select second Modflow Binary Heads file to process
  #=================================================================
  winA <- tktoplevel()
  msg = paste('Identify Binary Heads file for :', model)
  lbl.message <- tk2label(winA, text = msg, font = fontHeading)
  tkgrid(lbl.message, padx = 30)
  tkraise(winA)
  MFmodel.Params[model,]$mpath <- 
    '\\\\ad.sfwmd.gov\\dfsroot\\data\\wsd\\SUP\\proj\\CFWI_WetlandStress\\Update2018\\ModelRuns\\*.*'
  mpath <- toString(MFmodel.Params[model,]$mpath)
  SIMheadsFile<-choose.files(default=mpath)
  tkdestroy(winA)
  
  if (length(RCheadsFile) == 0 || length(SIMheadsFile) ==0) {
    exit("User cancelled HeadsFile choices")
  }
  
  #=================================================================
  # Estimate number of stress periods in Heads file
  #=================================================================
  fileSz1 <- file.info(RCheadsFile)$size
  fileSz2 <- file.info(SIMheadsFile)$size
  fileSz <- min(fileSz1,fileSz2)
  TtlStrPd = fileSz / ( M$nlays * ((M$ncols * M$nrows * 4) + 44))
  
  #=================================================================
  # Define range of Stress Periods to read
  #=================================================================
  SP_rng <- readRange()
  if (max(SP_rng) > TtlStrPd || min(SP_rng) < 1) {
    exit('Out of Range')
  }                     
}

#=================================================================
# Retrieve Heads for Layer1 from 
# Reference Condition (RC) and Simulation (SIM) Runs
# 2 files are read in asyncronously using the future package
#=================================================================
if (is.null(MFLay) ){
  MFLay <-1
}
maxSP <- as.integer(TtlStrPd)
plan(multiprocess)
processed= listenv(NULL)

#====================================================================
# tic() Initiates stacked timers and and toc() echos elapsed time
#====================================================================

tic("Modflow Binary Heads Data Processing")
tic("Heads Retrieval")
cat(paste("Initiating call to readHeadsbinByLay for Layer ", 
          MFLay, " as Reference ",
          "Condition  [+]\nwith input from ", RCheadsFile, '\n'))
processed[[1]] <- future({readHeadsbinByLay(RCheadsFile, 
                                            MFLay, maxSP)})

cat(paste("Initiating call to readHeadsbinByLay for Layer ", 
          MFLay, " as Model Simulation ",
          "of Interest  [:]\nwith input from ", SIMheadsFile, '\n'))
processed[[2]] <- future({readHeadsbinByLay(SIMheadsFile, 
                                            MFLay, maxSP)})

#====================================================================
# Wait for values from future with progress indicators
#====================================================================
cat(paste('Waiting for background processing to complete','\n'))

while (!resolved(processed[[1]])){
  if (!resolved(processed[[2]])){
    cat("+")
  }
  cat(":")
}
cat("\n")

#====================================================================
# Reformat Layer1 as 3D array using col, row, StressPeriod dimensions
#====================================================================

Layer1RC2d <- array(future::value(processed[[1]]),c(M$ncols,M$nrows,maxSP))
Layer1SIM2d<- array(future::value(processed[[2]]),c(M$ncols,M$nrows,maxSP))
toc()

#====================================================================
# Process P80 calculations for each model cell in parallel
#====================================================================
tic("P80 Calculations")
cat(paste('Initiating Percentile rank calculations','\n'))

qRC <- future_apply (Layer1RC2d,MARGIN=c(1,2),
                    FUN=stats::quantile,probs=c(.2),na.rm=T)
qSIM <- future_apply (Layer1SIM2d,MARGIN=c(1,2),
                     FUN=stats::quantile,probs=c(.2),na.rm=T)

toc()
toc()
