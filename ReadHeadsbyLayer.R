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

if (!"rModflow" %in% installed.packages()[,"Package"]) {
  devtools::install_github("KevinRodberg/rModflow")
}
source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ReusableFunctions/tclFuncs.R")

#=================================================================
# Beginning of Script
#
# Created by Kevin A. Rodberg - February 2018
#
# Purpose: 
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
wellfilePath = '//whqhpc01p/hpcc_shared/dbandara/CFWI/ECFTX/Model/Transient/TR11a/Postprocess/HuaR_HOB/*.*'
wellfile<-choose.files(default=wellfilePath)
path <- choose.dir(default='Y:/')
PointNames <- read.table(wellfile, header = TRUE, sep = ",")
colNames<- c('Row','Col','Lay','HOB_ID')
indx<-lapply(sapply(colNames,toupper),grep,sapply(names(PointNames),toupper))
Kpnts<-nrow(PointNames)
# meltPoints vector should be organized as Columns, Rows, Layers order
# to match format of Modflow binary heads array when read in.
meltPoints<-as.vector(melt(PointNames[,c(indx$Col[1],indx$Row[1],indx$Lay[1])])$value)
nameLkup<- PointNames[,c('HOB_ID','commonName')]

# PointNames should be organized Name, Layer, Row, Column
PointNames<-PointNames[,c(indx$HOB_ID[1],indx$Lay[1],indx$Row[1],indx$Col[1])]
PointVector<-array(meltPoints,c(Kpnts,3))

#===============================================
# Retrieve Heads by Layer
#===============================================
tic("Heads file processing")
pointValues <- readHeadsbinAtPnts(to.read, SP_rng, PointVector)
HeadRecs<-cbind(do.call(rbind,replicate(length(SP_rng),PointNames, simplify=FALSE)),pointValues)
names(HeadRecs)<-c("Station","LAY","ROW","COL","Head","SP")
SimHeadsWide <- dcast(HeadRecs,Station+ROW+COL+LAY ~SP,value.var="Head")
close(to.read)

toc()

hobpr<-read.table('//whqhpc01p/hpcc_shared/dbandara/CFWI/ECFTX/Model/Transient/TR11a/Final_Calib_20190208/TR11.hobpr',header = TRUE)
OBSHeadWide<-dcast(hobpr,OBSERVATION.NAME~rowid(OBSERVATION.NAME),value.var="OBSERVED.VALUE")
OBSheads<-melt(OBSHeadWide,"OBSERVATION.NAME")
names(OBSheads)<-c('Station','SP','Head')
OBSheads[OBSheads$Head < -999.0,]<-NA

SIMheads <- HeadRecs[,c('Station','SP','Head')]
OBSvSIMheads<- merge(OBSheads,SIMheads, 
                     by.x=c('Station','SP'),by.y=c('Station','SP'),
                     all.y=TRUE)

names(OBSvSIMheads)<-c('Station','SP','OBS','SIM')
OBSvSIMheads <- melt(OBSvSIMheads,id.vars=c('Station','SP'))
OBSvSIMheads$SP<- as.numeric(OBSvSIMheads$SP)
OBSvSIMheads<-OBSvSIMheads[order(OBSvSIMheads$Station,OBSvSIMheads$SP,OBSvSIMheads$variable),]

plotLines <- function(fileName,prmt,commonName,data,msg){
  cat(prmt,'\n')
  graphics.off()
  # p<-ggplot(subset(data,Station %in% prmt)) +
  p<-ggplot(data) +
    geom_path(aes(SP,value, col=variable))  +
    scale_x_continuous(name= "Stress Periods (monthly)", limits=c(0,133),
                       breaks = seq(1,133,12)) +
    labs(title=commonName,y = "Heads\n(feet NAVD 88)") +
    guides(color=guide_legend(title="Water Level Heads"))
  
  ggsave(filename=fileName,width=10,height=6.66,units="in",dpi=300)
}

#--
#  Set environment for mutliprocessing
#--
plan(multiprocess)
results <- listenv()

tic("Figures processing")
x=0
stationList = unique(OBSvSIMheads$Station)
ilen = length(unlist(stationList))
  for (station in stationList[1:25]){
  x= x + 1
  cat(paste('\r', x, format(x / ilen * 100, digits = 2, nsmall = 2), "%"))
    
  fileName= paste0(path,'/',station,'.png')
  commonName = nameLkup[nameLkup$HOB_ID==station,]$commonName
  results[[x]] <- future({plotLines(fileName,station,commonName,
                                    OBSvSIMheads[OBSvSIMheads$Station==station,])})
  
}
graphics.off()
toc()
# 
# compareHeads<-merge(SimHeadsWide,OBSHeadWide,by.x="Station", by.y="OBSERVATION.NAME")
# test=melt(compareHeads,c("Station","ROW","COL","LAY"))
# names(OBSHeadWide)
