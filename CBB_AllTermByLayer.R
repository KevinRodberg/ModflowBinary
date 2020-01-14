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

plan(multiprocess, .skip = TRUE)

#=================================================================
# R:\ModflowBinary\CBB_AllTermByLayer.R
#=================================================================
# Beginning of Sum Modflow CBC
#
# Created by Kevin A. Rodberg - October 2018
# Revised January 2020
#
# Purpose: Create matrix of mean All CBC terms and each layer for all simulation 
#          stress periods from Modflow Binary data

#=================================================================
# Choose Modflow Binary CBC file
#=================================================================
source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ReusableFunctions/tclFuncs.R")

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
basePath=choose.dir(default = "basePath", caption = "Select folder")

makeSumLayer<- function(cbbFile,selectLayer,SP_rng,term,basePath){
  to.read <- file(cbbFile, "rb")
  filPtr<- to.read 
  print (paste("Retrieving", trimws(term)))
  term2Read <-trimws(term)
  
  maxSP <- as.integer(TtlStrPd)
  CBCdata1 <- readCBCbinByTerm(to.read, term, SP_rng, selectLayer)

  #toc()
  yourTheme = rasterTheme(region = brewer.pal('BrBG', n = 9))
  # Reformat CBCdata1 as 3D array using col, row, StressPeriod dimensions
  # create dataframe of Head values for this Layer 
  CBCsMatrix<- array(CBCdata1,c(M$ncols,M$nrows,maxSP))
  negData <-CBCdata1
  posData <-CBCdata1
  negData[ negData >0 ] <- 0
  negCBCsMatrix<- array(negData,c(M$ncols,M$nrows,maxSP))
  posData[ posData <0 ] <- 0
  posCBCsMatrix<- array(posData,c(M$ncols,M$nrows,maxSP))
  outBuds<-apply(negCBCsMatrix, c(3),sum)
  InBuds<-apply(posCBCsMatrix, c(3),sum)
  #tic(paste('Creating sum raster for layer ',selectLayer))
  #xf <-future_apply (CBCsMatrix,MARGIN=c(1,2),FUN=stats::quantile,probs=c(.1),na.rm=T)
  xf <-future_apply (CBCsMatrix,MARGIN=c(1,2),FUN=sum,na.rm=T)
  xfIn <-sum(InBuds)
  xfOut <-sum(outBuds)
  xfNet <-xfIn+xfOut

  # my.at = seq(-80,100,15)
  # levelplot(raster(t(xf[,])),par.settings = yourTheme,at=my.at)
  
  SumSim <- t(xf[,])

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
  
  SumRasL1<-raster::raster(SumSim,xmin,xmax,ymin,ymax, crs=HARNSP17ft)
  my.at = c(minValue(SumRasL1),-1,1,maxValue(SumRasL1))
  levelplot(SumRasL1,par.settings = yourTheme,at=my.at)
  Sys.sleep(0)
  filename = paste0(basePath,'/Sum_',term2Read,'_Lay',selectLayer,'.tif')
  raster::writeRaster(SumRasL1, filename, format="GTiff", overwrite=TRUE)
  budgets <-list("term"=term,"Layer"=selectLayer,"In"=xfIn,"Out"=xfOut,"Net"=xfNet,"outBuds"=outBuds,"InBuds"=InBuds)
  return(budgets)
}

ix=0
tic("Budget Creation")
allBudgets <- list()
 for (term in CBCterms){
   data <- listenv()
   for (selectLayer in seq(1,M$nlays)){
     #===============================================
     # Retrieve CBC by Layer
     #===============================================
     cat(paste("Start making rasters for",term,"Layer",selectLayer, '\n'))
     ix = ix + 1
     data[[ix]] %<-%  makeSumLayer(cbbFile,selectLayer,SP_rng,term,basePath)
     #data[[ix]] <-  makeSumLayer(cbbFile,selectLayer,SP_rng,term,basePath)
   }
   cat(paste("Rasters will be created in ",basePath,"in a couple minutes \n"))
    lastOne <- futureOf(data[[ix]])
    while (!resolved(lastOne)) {
      cat ('-')
      Sys.sleep(0.5)    
    }

   #-------------------------------------------------
   # unlist results returned from FUNCTION "biasByYearMon"
   #-------------------------------------------------
   budTerm = trim(unlist(lapply(data,"[[",1)))
   Layer = unlist(lapply(data,"[[",2))
   In = unlist(lapply(data,"[[",3))
   Out = unlist(lapply(data,"[[",4))
   Net = unlist(lapply(data,"[[",5))
   outBuds =unlist(lapply(data,"[[",6))
   InBuds = unlist(lapply(data,"[[",7))
   cols <- c(paste0('SP',SP_rng))
   Budget<-as.data.frame(cbind(budTerm,Layer,In,Out,Net))
   In_SP_Cols <- as.data.frame(matrix(unlist(InBuds),nrow=length(Layer),byrow=T))
   Out_SP_Cols <- as.data.frame(matrix(unlist(outBuds),nrow=length(Layer),byrow=T))
   Net_SP_Cols <- matrix(unlist(InBuds),nrow=length(Layer),byrow=T) + 
     matrix(unlist(outBuds),nrow=length(Layer),byrow=T)
   names(Out_SP_Cols)<-cols
   inBudget <-cbind(budTerm,Layer,In,Out,Net,In_SP_Cols)
   outBudget <-cbind(budTerm,Layer,In,Out,Net,Out_SP_Cols)
   netBudget <-cbind(budTerm,Layer,In,Out,Net,Net_SP_Cols)
   allBudgets[[ix]] <-Budget
   write.csv(file=paste0(basePath,'\\In_', trimws(term),'.csv'),inBudget)
   write.csv(file=paste0(basePath,'\\Out_', trimws(term),'.csv'),outBudget)
   write.csv(file=paste0(basePath,'\\Net_', trimws(term),'.csv'),netBudget)
   print(Budget[Budget$budTerm==trimws(term),])
 }
exportBudget<-do.call(rbind,allBudgets)
cat (paste('Source Data:',cbbFile,'\n'))
write.csv(file=paste0(basePath,'\\WaterBudgetSummary.csv'),exportBudget)
toc()


