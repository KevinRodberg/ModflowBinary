library (data.table)
library(tcltk2)
library(rModflow)

setwd('R:/ModflowBinary/sharable')
source ("./tclFuncs.R")
source('R:/rModflow/R/modflowDataFuncs.R')

#=================================================================
# Beginning of Script
#
# Created by Kevin A. Rodberg - February 2018
#
# Purpose: Create series of raster figures to review Modflow Binary data
#          and large ET or Recharge datasets/

#=================================================================
# Choose Modflow Cell by Cell Budget Term file
#=================================================================
MFmodel.Params <- defineMFmodel()
model <- chooseModel()

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

options <- chooseBudgetTerms(CBCterms)
n1=options$n1
n2=options$n2

#===============================================
# Retrieve first Budget Term Data Set Layer = 1
#===============================================
to.read <- file(cbbFile, "rb")
Layer <- 1
print (paste("Retrieving", trimws(CBCterms[[n1]])))
CBCdata1 <- readCBCbinByTerm(to.read, CBCterms[[n1]], SP_rng, Layer)
close(to.read)

#===============================================
# Retrieve second Budget Term Data Set Layer = 1
#===============================================
to.read <- file(cbbFile, "rb")
Layer <- 1
print (paste("Retrieving", trimws(CBCterms[[n2]])))
CBCdata2 <- readCBCbinByTerm(to.read, CBCterms[[n2]], SP_rng, Layer)
close(to.read)

#===============================================
# Add Budget terms n1 + n2
#===============================================
DiffVector <- CBCdata1 + CBCdata2
print(paste("Adding: ", trimws(CBCterms[[n1]]), " + ", trimws(CBCterms[[n2]])))
nsp = length(SP_rng)
