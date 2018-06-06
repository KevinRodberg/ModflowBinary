library (data.table)
library(rModflow)
#===============================================
# Function to exit a little more nicely
#===============================================
exit <- function(msg){
  print(paste0("ERROR:", msg))
  closeAllConnections()
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}

#===============================================
# Prompt for user to enter an integer
#===============================================
readinteger <- function(){
  n <- readline(prompt = "Enter an integer: ")
  return(as.integer(n))
}

#===============================================
# Accept string defining range of integers vals
# reform as a unique sequence
#===============================================
readRange <- function(){
  rngStr <- readline(prompt = paste("Enter range >"))
  print (rngStr)
  rngStr <- gsub(" ", ",", rngStr)
  rngStr <- gsub(",,", ",", rngStr)
  print (rngStr)
  
  df <- as.vector(rngStr)
  rng <-
    sapply(df, function(x)
      dget(textConnection(paste('c(', x, ')'))))
  rng <- unique(sort(rng))
  return(rng)
}

#===============================================
# Define Modflow Cell by Cell Budget Term file
#===============================================

model = 'NPALM'
model <- readline(prompt = paste("Enter mode choice [NPALM,ECFTX,LWCSIM] >"))
if (model == 'NPALM') {
  cbbFile <-
    choose.files(default = "\\\\whqhpc01p\\hpcc_shared\\krodberg\\NPB\\*.*")
} else if (model == 'ECFTX')   {
  cbbFile <-
    choose.files(default = "\\\\whqhpc01p\\hpcc_shared\\dbandara\\CFWI\\ECFTX\\Model\\Transient\\*.*")
} else if (model == 'LWCSIM') {
  cbbFile <-
    choose.files(default = "\\\\ad.sfwmd.gov\\dfsroot\\data\\wsd\\MOD\\LWCSIM\\*.*")
  
} else {exit('Incorrect model choice')}

#cbbFile <-'//whqhpc01p/hpcc_shared/krodberg/NPB/FWO/fort.700'
#cbbFile <-'//ad.sfwmd.gov/dfsroot/data/wsd/MOD/RCE/101816_Test33_PESTRUN/LWCSIM_MF_Test33.cbb'
#cbbFile <-'//whqhpc01p/hpcc_shared/dbandara/CFWI/ECFTX/Model/Transient/Run_UB_update/ecftx.cbb'
to.read = file(cbbFile, "rb")

#===============================================
# Scan CBC file for avaialble Budget Terms
# and model characteristics
#===============================================
CBCTermSet = listBinHeaders(to.read)
ncols <- as.numeric(CBCTermSet[4])
nrows <- as.numeric(CBCTermSet[3])
nlays <- as.numeric(CBCTermSet[2])
CBCterms <- unlist(CBCTermSet[[5]])

#===============================================
# Calculate number of stress periods in CBC file
#===============================================
fileSz <- file.info(cbbFile)$size
TtlStrPd = fileSz / ((length(CBCterms) * ((
  ncols * nrows * nlays * 4
) + 36)))

#===============================================
#
#===============================================
cat(paste0("Total Stress Periods =", TtlStrPd,'\n\n\n'))
cat("Choose Periods of interest \n")
cat("i.e.: 1:3,5,7:100,200 \n")
SP_rng <- readRange()
if (max(SP_rng) > TtlStrPd || min(SP_rng) < 1) {
  exit('Out of Range')
}

#===============================================
# Display List of available Budget Terms in CBC file
#===============================================
for (i in seq(1, length(CBCterms))) {
  cat (paste(i, ':', trimws(CBCterms[[i]]), '\n'))
}
close(to.read)

#===============================================
# Select budget Term to extract from CBC file
#===============================================
print(n1 <- readinteger())
if(n1>length(CBCterms)){exit('Invalid Budget Term Number')}

#===============================================
# Select budget Term to extract from CBC file
#===============================================
print(n2 <- readinteger())
if(n2>length(CBCterms)){exit('Invalid Budget Term Number')}

#===============================================
# Retrieve first Budget Term Data Set
#===============================================
to.read <- file(cbbFile, "rb")

print (paste("Retrieving", trimws(CBCterms[[n1]])))
CBCdata1 <-
  readCBCbinByTerm(to.read, CBCterms[[n1]], SP_rng)
close(to.read)

#===============================================
# Retrieve second Budget Term Data Set
#===============================================
to.read <- file(cbbFile, "rb")

print (paste("Retrieving", trimws(CBCterms[[n2]])))
CBCdata2 <-
  readCBCbinByTerm(to.read, CBCterms[[n2]], SP_rng)
close(to.read)

#===============================================
# Subtract Budget term n1 - n2
#===============================================

DiffVector <- CBCdata1 - CBCdata2
print(paste("Subtracting: ",trimws(CBCterms[[n1]]), " - ", trimws(CBCterms[[n2]])))
nsp = length(SP_rng)
DiffVectorBySP <-  array(DiffVector, dim = c(ncols, nrows, nlays, nsp))
rchVectorBySP <- array(CBCdata1, dim=c(ncols,nrows,nlays,nsp))
ETVectorBySP <- array(CBCdata2, dim=c(ncols,nrows,nlays,nsp))

print(paste("max value:",max(ETVectorBySP)))

source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ModflowBinary/netRCH_anim.R")