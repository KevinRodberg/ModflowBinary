#===================================================================================
#   Script Name:  PltMFHeadsVsHOB.R
#
#       Purpose:  Plot hydrographs comparing Modflow Binary Heads values 
#                 with Observed values
#===================================================================================
#  Developed by: Kevin A. Rodberg, Science Supervisor 
#                 Resource Evaluation Section, Water Supply Bureau, SFWMD
#                 (561) 682-6702
#
#--
#   package management:
#     provide automated means for first time use of script to automatically 
#	  install any new packages required for this code, with library calls 
#	  wrapped in a for loop.
#--

pkgChecker <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}

list.of.packages <-c( "githubinstall","data.table","tcltk2","rModflow","devtools",
                      "tictoc","tidyr","dplyr","plyr","future", 
                      "listenv","ggplot2")

suppressWarnings(pkgChecker(list.of.packages))
tic("All Processes:")
if (!"rModflow" %in% installed.packages()[,"Package"]) {
  devtools::install_github("KevinRodberg/rModflow")
}
source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ReusableFunctions/tclFuncs.R")

getBinaryData<- function(headsFile, SP_rng, PointVector,M) {
  #
  # Retrieves a timeseries of Modflow head values for specific model cells
  #
  # Args:
  #   headsFile:    Modflow binary heads file name
  #   SP_rng:       User selected Stress Period Range defined in tclFuncs()
  #   PointVector:  Array of model cell indices organized by (Lay,Row,Col) to select cells rom 3D Heads matrix
  #   M:            M <- as.data.frame(MFmodel.Params[model,])
  #     
  # Returns:
  #   A timeseries series of head values at specified model cell locations
  #     
  
  to.read = file(headsFile, "rb")
  pointValues <- rModflow::readHeadsbinAtPnts(to.read, SP_rng, PointVector)
  close(to.read)
  return(pointValues)
}

#===================================================================================
# Choose SFWMD Modflow Model from predined data in rModflow::defineMFmodel()
# and Modflow Binary Heads file using tk GUI interface with utils::choose.files
#===================================================================================
MFmodel.Params <- rModflow::defineMFmodel()
model <- rModflow::chooseModel()
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


#===================================================================================
# Estimate number of stress periods in Heads file
# and select range of Stress Periods to read
#===================================================================================
fileSz <- file.info(headsFile)$size
TtlStrPd = fileSz / ( M$nlays * ((M$ncols * M$nrows * 4) + 44))
SP_rng <- readRange()
if (max(SP_rng) > TtlStrPd || min(SP_rng) < 1) {
  exit('Out of Range')
}

#===================================================================================
# choose & read wellfile formated as "HOB_ID" , "LAYER"  , "Row"  ,   "Column_"
# with utils::choose.files
#===================================================================================
wellfilePath = paste0('//whqhpc01p/hpcc_shared/dbandara/CFWI/ECFTX/Model/Transient/',
                      'TR11a/Postprocess/HuaR_HOB/*.*')
wellfile<-choose.files(default=wellfilePath)

PointNames <- read.table(wellfile, header = TRUE, sep = ",")
colNames<- c('Row','Col','Lay','HOB_ID','Lat','Lon')
indx<-lapply(sapply(colNames,toupper),grep,sapply(names(PointNames),toupper))
Kpnts<-nrow(PointNames)

#===================================================================================
# Provide crosswalk to HOB_id for figure titles
#===================================================================================
nameLkup<- PointNames[,c('HOB_ID','commonName')]

#===================================================================================
# meltPoints vector should be organized as Columns, Rows, Layers order
# to match format of Modflow binary heads array to be read in by 
# rModflow::readHeadsbinAtPnts.
#===================================================================================
meltPoints<-as.vector(melt(PointNames[,c(indx$Col[1],indx$Row[1],indx$Lay[1])])$value)

#===================================================================================
# PointNames should be organized Name, Layer, Row, Column to match
# organization of pointValues returned from rModflow::readHeadsbinAtPnts
#===================================================================================
PointNames<-PointNames[,c(indx$HOB_ID[1],indx$Lat[1],indx$Lon, 
                          indx$Lay[1],indx$Row[1],indx$Col[1])]
PointVector<-array(meltPoints,c(Kpnts,3))

#===================================================================================
# Retrieve Heads selected by PointVector (Lay,Row,Col) in a background process)
#   non-"future" method: 
        # to.read = file(headsFile, "rb")
        # pointValues <- rModflow::readHeadsbinAtPnts(to.read, SP_rng, PointVector)
        # close(to.read)
#===================================================================================
tic("Binary Heads file processing")

#--
#  Set environment for mutliprocessing
#--
plan(multiprocess,workers=availableCores())

processed= listenv(NULL) # for reading binary heads at points
results <- listenv(NULL) # for plot routine

processed[[1]] <- future({getBinaryData(headsFile, SP_rng, PointVector,M)})

tic("Heads OBS file processing")
hobsfilePath = '//whqhpc01p/hpcc_shared/dbandara/CFWI/ECFTX/Model/Transient/TR11a/*.hobpr'
hobsfile<-choose.files(default=hobsfilePath)
hobpr<-read.table(hobsfile,header = TRUE)
OBSHeadWide<-dcast(hobpr,OBSERVATION.NAME~rowid(OBSERVATION.NAME),value.var="OBSERVED.VALUE")
OBSheads<-melt(OBSHeadWide,"OBSERVATION.NAME")
names(OBSheads)<-c('Station','SP','Head')
OBSheads[OBSheads$Head < -999.0,]<-NA
toc()
cat('\nChoose output directory for figures...\n')
path <- tk_choose.dir(default='//ad.sfwmd.gov/dfsroot/data/wsd/SUP/proj')
#====================================================================
# Wait for values from future with progress indicators
#====================================================================
cat(paste('\nWaiting for background Binary Heads processing to complete','\n'))

while (!resolved(processed[[1]])){
  cat("+")
}
cat("\n")
toc()

pointValues <- future::value(processed[[1]])
HeadRecs<-cbind(do.call(rbind,replicate(length(SP_rng),PointNames, simplify=FALSE)),pointValues)
names(HeadRecs)<-c("Station","LAY","ROW","COL","Head","SP")


# SimHeadsWide <- dcast(HeadRecs,Station+ROW+COL+LAY ~SP,value.var="Head")

SIMheads <- HeadRecs[,c('Station','SP','Head')]
OBSvSIMheads<- merge(OBSheads,SIMheads, 
                     by.x=c('Station','SP'),
                     by.y=c('Station','SP'),
                     all.y=TRUE)

names(OBSvSIMheads)<-c('Station','SP','OBS','SIM')
OBSvSIMheads <- melt(OBSvSIMheads,id.vars=c('Station','SP'))
OBSvSIMheads$SP<- as.numeric(OBSvSIMheads$SP)
OBSvSIMheads<-OBSvSIMheads[order(OBSvSIMheads$Station,OBSvSIMheads$SP,OBSvSIMheads$variable),]

plotLines <- function(fileName,prmt,commonName,data,msg){
  cat(prmt,'\n')
  p<-ggplot(data) +
    geom_path(aes(SP,value, col=variable))  +
    scale_x_continuous(name= "Stress Periods (monthly)", 
                       limits=c(0,133),
                       breaks = seq(1,133,12)) +
    labs(title=commonName,y = "Heads\n(feet NAVD 88)") +
    guides(color=guide_legend(title="Water Level Heads"))
  
  ggsave(filename=fileName,width=10,height=6.66,units="in",dpi=300)
  graphics.off()
}



tic("Figures processing")
cat(paste('\nCreating hydrographs....\n'))

graphics.off()

x=0
stationList = unique(OBSvSIMheads$Station)
ilen = length(unlist(stationList))
for (station in stationList[1:25]){
# for (station in stationList[1:25]){
  x= x + 1
  cat(paste('\r', x, format(x / ilen * 100, digits = 2, nsmall = 2), "%"))
  
  fileName= paste0(path,'/',station,'.png')
  commonName = as.character(nameLkup[!is.na(nameLkup$HOB_ID) & nameLkup$HOB_ID==station,]$commonName)
  results[[x]] <- future({plotLines(fileName,station,commonName,
                                    OBSvSIMheads[OBSvSIMheads$Station==station,])})
}


#====================================================================
# Wait for values from future with progress indicators
#====================================================================
cat(paste('Waiting for background figure processing to complete','\n'))

while (!resolved(results[[x]])){
  cat("+")
}
cat("\n")
toc()
toc()
# 
# compareHeads<-merge(SimHeadsWide,OBSHeadWide,by.x="Station", by.y="OBSERVATION.NAME")
# test=melt(compareHeads,c("Station","ROW","COL","LAY"))
# names(OBSHeadWide)
