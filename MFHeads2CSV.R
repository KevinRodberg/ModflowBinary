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

suppressWarnings(pkgChecker(c("optparse")))


makeFIGS <- FALSE
model <- "LWCSIM" 
defOutput <- getwd()
headsFile <- paste0("//whqhpc01p/hpcc_shared/dbandara/LWCSIM/Model/Transient/99_14_ManualCalib/",
                    "99_14_Final_manual_042519.hds")
#===================================================================================
# Non-interactive Alternative to batch file (or test command line arguments)
#   with execution from the console.
#
#   system("rscript r:/ModflowBinary/MFHeads2CSV.R -h")
# or
#   system("rscript r:/ModflowBinary/MFHeads2CSV.R --mod LWCSIM --out R:/ModflowBinary/figureOutput/.")
#===================================================================================
if (!interactive()) {
  parser = OptionParser();
  parser<-add_option(parser,c("-m", "--mod"), 
                     action = "store",
                     dest =   'model',
                     type =   "character", 
                     default ='LWCSIM', 
                     help =   "Model name", 
                     metavar= "character")
  
  parser<-add_option(parser,c("-f", "--figs"), 
                     action = "store_true",
                     type="logical", 
                     dest = 'makeFIGS',
                     default=FALSE, 
                     help="Create figures",
                     metavar="character")
  
  parser<-add_option(parser,c("-o", "--out"), 
                     type="character", 
                     default=".", 
                     help="output file name [default= %default]", 
                     metavar="character")
  
  parser<-add_option(parser,c("-b", "--bin"), 
                     type="character", 
                     default=headsFile, 
                     help=paste("Binary heads file [default= ",headsFile,"]"), 
                     metavar="character")
  
  opt = parse_args(parser);
  
  model     <- opt$model
  headsFile <- opt$bin
  defOutput <- opt$out
}

interactive <- interactive()


list.of.packages <-c( "githubinstall","data.table","tcltk2","rModflow","devtools",
                      "tictoc","tidyr","dplyr","plyr","future", 
                      "listenv","ggplot2")

suppressWarnings(pkgChecker(list.of.packages))
tic.clear()
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

SP_rng <- seq(1,192)
#===================================================================================
# Choose SFWMD Modflow Model from predined data in rModflow::defineMFmodel()
# and Modflow Binary Heads file using tk GUI interface with utils::choose.files
#===================================================================================
MFmodel.Params <- rModflow::defineMFmodel()

if (interactive()){
  model <- rModflow::chooseModel()
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
}
M <- as.data.frame(MFmodel.Params[model,])

#===================================================================================
# Estimate number of stress periods in Heads file
# and select range of Stress Periods to read
#===================================================================================
fileSz <- file.info(headsFile)$size
TtlStrPd = fileSz / ( M$nlays * ((M$ncols * M$nrows * 4) + 44))
if (interactive){
  SP_rng <- readRange()
  if (max(SP_rng) > TtlStrPd || min(SP_rng) < 1) {
    exit('Out of Range')
  }
}

#===================================================================================
# choose & read wellfile formated as  "NAME","ROW","COL","LAY"
# with utils::choose.files
#===================================================================================
# wellfilePath = paste0('//whqhpc01p/hpcc_shared/dbandara/LWCSIM/Model/Postproces/',
#                       'ReadHeads/readHeads_sw2k/*.*')
# wellfile<-choose.files(default=wellfilePath)

PointNames<-read.csv('//ad.sfwmd.gov/dfsroot/data/wsd/MOD/LWCSIM/LWCSIM_HOB_info4Shiny.csv')
colNames<- c('HOB_ID','Row','Col','Lay')

indx<-lapply(sapply(colNames,toupper),grep,sapply(names(PointNames),toupper))
Kpnts<-nrow(PointNames)

#===================================================================================
# meltPoints vector should be organized as Columns, Rows, Layers order
# to match format of Modflow binary heads array to be read in by 
# rModflow::readHeadsbinAtPnts.
#===================================================================================
cat ('Initializing datasets \n')
meltPoints<-as.vector(melt(PointNames[,c(indx$Col[1],indx$Row[1],indx$Lay[1])])$value)

#===================================================================================
# PointNames should be organized Name, Layer, Row, Column to match
# organization of pointValues returned from rModflow::readHeadsbinAtPnts
#===================================================================================
PointNames<-PointNames[,c(indx$HOB_ID[1],indx$Lay[1],indx$Row[1],indx$Col[1])]
PointVector<-array(meltPoints,c(Kpnts,3))

#===================================================================================
# Retrieve Heads selected by PointVector (Lay,Row,Col) in a background process)
#   non-"future" method: 
        # to.read = file(headsFile, "rb")
        # pointValues <- rModflow::readHeadsbinAtPnts(to.read, SP_rng, PointVector)
        # close(to.read)
#===================================================================================

#--
#  Set environment for mutliprocessing
#--
plan(multiprocess,workers=availableCores())

processed= listenv(NULL) # for reading binary heads at points
results <- listenv(NULL) # for plot routine

tic("Read Binary Heads file:")

processed[[1]] <- future({getBinaryData(headsFile, SP_rng, PointVector,M)})

hobsfile<-paste0('//whqhpc01p/hpcc_shared/dbandara/LWCSIM/Model/Postproces/ReadHeads/',
                 'readHeads_sw2k/readHeads_sw2k/LWCSIM_HOB_info_99_14_R.csv')
OBSHeadWide<-read.csv(hobsfile)

OBSheads<-melt(OBSHeadWide,id.vars=c("HOB_ID","Row","Column_","lyr201808"))

OBSheads$variable <- gsub('STP', '', OBSheads$variable)
OBSheads$variable<- as.numeric(OBSheads$variable)

names(OBSheads)<-c('Station','Row','Col','Lay','SP','Head')
OBSheads[OBSheads$Head == 0.0,]$Head<-NA

#====================================================================
# Wait for values from future with progress indicators
#====================================================================
cat(paste('\nWaiting for background Binary Heads processing to complete','\n'))

while (!resolved(processed[[1]])){
  cat("+")
}
cat("\n")
toc() # Timing for Binary Heads

pointValues <- future::value(processed[[1]])
HeadRecs<-cbind(do.call(rbind,replicate(length(SP_rng),PointNames, simplify=FALSE)),pointValues)

names(HeadRecs)<-c("Station","LAY","ROW","COL","Head","SP")

SIMheads <- HeadRecs[,c('Station','SP','Head')]

OBSvSIMheads<- merge(OBSheads,SIMheads, 
                      by.x=c('Station','SP'),
                      by.y=c('Station','SP'),
                      all.y=TRUE)
 
names(OBSvSIMheads)<-c('Station','SP', "Row","Col","Lay",'OBS','SIM')
OBSvSIMheads <- melt(OBSvSIMheads,id.vars=c('Station','SP', "Row","Col","Lay"))
OBSvSIMheads$SP<- as.numeric(OBSvSIMheads$SP)
OBSvSIMheads<-OBSvSIMheads[order(OBSvSIMheads$Station,OBSvSIMheads$SP,OBSvSIMheads$variable),]


hob<-dcast(OBSvSIMheads,Station+SP+Row+Col+Lay~variable)
hob$RESID<-hob$SIM-hob$OBS

hobinfo<-read.csv('//ad.sfwmd.gov/dfsroot/data/wsd/MOD/LWCSIM/LWCSIM_HOB_info4Shiny.csv')

newColumns <-c('nobs', 'mean_error', 'MAE', 'nMAE', 'rmse', 'mean_obs', 'max_obs', 
               'min_obs', 'sd_obs', 'range_obs', 'mean_sim', 'max_sim', 'min_sim', 
               'sd_sim', 'range_sim', 'R2', 'NS')

hobinfo[newColumns]<-NA
tic("Stats")
cat("Calc statistics\n")
Ntarget <-nrow(hobinfo)
xx= 0
nstations <- length(unique(hob$Station))
for (s in unique(hob$Station)){
  #cat('+')
  xx = xx + 1
  cat(paste('\r', xx, format(xx / nstations * 100, digits = 2, nsmall = 2), "%"))
  
  hobinfo[hobinfo$HOB_ID==s,]$nobs <- sum(!is.na(hob[hob$Station==s,]$RESID))
  hobinfo[hobinfo$HOB_ID==s,]$mean_error <- mean(hob[hob$Station==s,]$RESID,na.rm=T)
  hobinfo[hobinfo$HOB_ID==s,]$MAE <-  abs(hobinfo[hobinfo$HOB_ID==s,]$mean_error)       
  hobinfo[hobinfo$HOB_ID==s,]$rmse <- sqrt(mean(hob[hob$Station==s,]$RESID^2,na.rm = TRUE))
  hobinfo[hobinfo$HOB_ID==s,]$mean_obs <- mean(hob[hob$Station==s,]$OBS,na.rm=T)
  if (sum(!is.na(hob[hob$Station==s,]$OBS))>0){
    hobinfo[hobinfo$HOB_ID==s,]$min_obs <- min(hob[hob$Station==s ,]$OBS,na.rm=T)
    hobinfo[hobinfo$HOB_ID==s,]$max_obs <- max(hob[hob$Station==s,]$OBS,na.rm=T)
    hobinfo[hobinfo$HOB_ID==s,]$range_obs <- 
      max(hob[hob$Station==s,]$OBS,na.rm=T) - min(hob[hob$Station==s,]$OBS,na.rm=T)
  }
  if (hobinfo[hobinfo$HOB_ID==s,]$nobs>0 ){
    hobinfo[hobinfo$HOB_ID==s,]$nMAE <- 
      sum(abs(hob[hob$Station==s,]$RESID - mean(hob[hob$Station==s,]$RESID,na.rm = T)),
          na.rm = T) / hobinfo[hobinfo$HOB_ID==s,]$nobs
  }
  hobinfo[hobinfo$HOB_ID==s,]$mean_sim <- mean(hob[hob$Station==s,]$SIM,na.rm=T)
  hobinfo[hobinfo$HOB_ID==s,]$min_sim <- min(hob[hob$Station==s,]$SIM,na.rm=T)
  hobinfo[hobinfo$HOB_ID==s,]$max_sim <- max(hob[hob$Station==s,]$SIM,na.rm=T)   
  hobinfo[hobinfo$HOB_ID==s,]$sd_obs <- sd(hob[hob$Station==s,]$OBS,na.rm=T)
  hobinfo[hobinfo$HOB_ID==s,]$sd_sim <- sd(hob[hob$Station==s,]$SIM,na.rm=T)
  hobinfo[hobinfo$HOB_ID==s,]$range_sim <- 
    max(hob[hob$Station==s,]$SIM,na.rm=T) - min(hob[hob$Station==s,]$SIM,na.rm=T)
  ## Calculate R2 and Nash-Sutcliffe coefficient
  if (hobinfo[hobinfo$HOB_ID==s,]$nobs>3 && hobinfo[hobinfo$HOB_ID==s,]$sd_obs>0 ){
    hobinfo[hobinfo$HOB_ID==s,]$R2 <- 
      cor(hob[hob$Station==s,]$SIM, hob[hob$Station==s,]$OBS,use="na.or.complete")^2
    hobinfo[hobinfo$HOB_ID==s,]$NS<- 
      1 - sum(hob[hob$Station==s,]$RESID^2, na.rm = TRUE)/
      sum((hob[hob$Station==s,]$OBS - hobinfo[hobinfo$HOB_ID==s,]$mean_obs)^2, na.rm = TRUE)
  }
}
cat('\n')
toc() # Timing for Stats

plotLines <- function(fileName,prmt,commonName,statMsg,data){
  cat(prmt,'\n')
  sp<-length(data[data$variable=='SIM',]$value)
  p<-ggplot(data) +
    geom_path(aes(SP,value, col=variable))  +
    # geom_path(aes(SP,Head))  +
    scale_x_continuous(name= "Stress Periods (monthly)", 
                       limits=c(0,sp+1),
                       breaks = seq(0,sp,12)) +
    labs(title=commonName,subtitle=statMsg,y = "Heads\n(feet NAVD 88)") +
    guides(color=guide_legend(title="Water Level Heads")) +
    theme(plot.subtitle = element_text(family="mono",size=9))
  
  ggsave(filename=fileName,width=10,height=6.66,units="in",dpi=300)
  graphics.off()
}

path <- defOutput

if (makeFIGS == TRUE) {
  if (interactive){
    cat('\nChoose output directory for figures...\n')
    path <- tk_choose.dir(default=defOutput)
  }

  tic("Figures processing")
  cat(paste('\nCreating hydrographs....\n'))
  graphics.off()
  x=0
  stationList = unique(OBSvSIMheads$Station)
  #stationList = unique(SIMheads$Station)
  ilen = length(unlist(stationList))
  for (station in stationList){
    # for (station in stationList[1:25]){
    x= x + 1
    cat(paste('\r', x, format(x / ilen * 100, digits = 2, nsmall = 2), "%"))
    statMsg<-paste0(paste(format(names(hobinfo)[23:32],width=10,justify="right"),
          collapse="|"),'\n',
    paste(format(hobinfo[hobinfo$HOB_ID == station,c(23:32) ],
                 digits = 4, nsmall = 2,width=10,justify="right"), 
          collapse="|"),'\n',
    paste(format(names(hobinfo)[33:39],width=10,justify="right"),
          collapse="|"),'\n',
    paste(format(hobinfo[hobinfo$HOB_ID == station,c(33:39) ],
                 digits = 4, nsmall = 2,width=10,justify="right"), 
          collapse="|"))
    fileName= paste0(path,'/',station,'.png')
    # commonName = as.character(nameLkup[!is.na(nameLkup$HOB_ID) & nameLkup$HOB_ID==station,]$commonName)
    commonName = station
    results[[x]] <- future({
      plotLines(fileName,station,commonName,statMsg,
                OBSvSIMheads[OBSvSIMheads$Station==station,])
      })
  }
  
  #====================================================================
  # Wait for values from future with progress indicators
  #====================================================================
  cat(paste('Waiting for background figure processing to complete','\n'))
  
  while (!resolved(results[[x]])){
    cat("+")
  }
  cat("\n")
  toc() # Timing for figures
}

hobwide<-merge(merge(dcast(hob,Station~SP,value.var='SIM'),
               dcast(hob,Station~SP,value.var='OBS'), by='Station'), 
               dcast(hob,Station~SP,value.var='RESID'), by = 'Station')

NSP = length(SP_rng)
names(hobwide)<- c('Station',paste0("SIM_SP",1:NSP), 
                   paste0("OBS_SP",1:NSP),
                   paste0("RESID_SP",1:NSP))

hobpr<-merge(hobwide,hobinfo,by.x ="Station",by.y = "HOB_ID")

hobpr$FID <- NULL
names(hobpr)[names(hobpr) == 'Station'] <- 'HOB_ID'

# write.csv(hobpr,"//whqhpc01p/hpcc_shared/dbandara/LWCSIM/Model/Postproces/ShinyApp_1999_2014/Data/output/LWCSIM_HOB_out_99_14_051119.csv",quote = FALSE, row.names = FALSE)
write.csv(hobpr,"//ad.sfwmd.gov/dfsroot/data/wsd/MOD/LWCSIM/LWCSIM_HOB_out_99_14_051719.csv",quote = FALSE, row.names = FALSE)
write.csv(hobpr,paste0(path,"/LWCSIM_HOB_out_99_14_051719.csv"),quote = FALSE, row.names = FALSE)

toc() #Timing for ALL Processes

