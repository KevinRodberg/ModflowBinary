library (data.table)
library(tcltk2)

source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ResuableFunctions/tclFuncs.R")
source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ResuableFunctions/modflowDataFuncs.R")

#=================================================================
# Beginning of Script
#
# Created by Kevin A. Rodberg - February 2018
#
# Purpose: Create series of raster figures to review Modflow Binary data
#          and large ET or Recharge datasets/

#=================================================================
# Choose Modflow Binary Heads file
#=================================================================
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

wellfile = "//whqhpc01p/hpcc_shared/jgidding/LECSR/LOX18/ALT2/wellinfo.prn.stage"
PointNames <- read.table(wellfile, header = TRUE, sep = "")
rownames(PointNames)<-PointNames$Station
PointNames$Station <- NULL
Kpnts<-nrow(PointNames)
PointNames<-PointNames[,c(2,1,3)]
names(PointNames)
meltPoints<-as.vector(melt(PointNames)$value)
PointVector<-array(meltPoints,c(Kpnts,3))

#===============================================
# Retrieve Heads by Layer
#===============================================
to.read <- file(headsFile, "rb")
pointValues <- readHeadsbinAtPnts(to.read, SP_rng, PointVector)

close(to.read)

