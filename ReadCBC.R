library (data.table)
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
# Read just the Header record from a Cell by Cell file
#===============================================
readCBCHeader <- function(filPtr) {
  while (length(record <- readBin(filPtr, raw(), 36)) > 0)
  {
    ints <- readBin(record, integer(), 9)
    txt <- intToUtf8(record[8L + seq(16)])
    KSTP <- ints[1]
    KPER <- ints[2]
    NC <- ints[7]
    NR <- ints[8]
    K <- ints[9]
    header <- list(
      KSTP = KSTP,
      KPER = KPER,
      TEXT = txt,
      NC = NC,
      NR = NR,
      K = K
    )
    return(header)
  }
  return(list())
}
#===============================================
#  Create list of Budget Term Headers available in 
# Modflow Cell by Cell binary file
#===============================================
listBinHeaders <- function(filPtr) {
  #  CBCterms <- vector("list", 100)
  CBCterms <- list()
  firstHeader <- readCBCHeader(filPtr)
  kntFloats <- firstHeader$K * firstHeader$NR * firstHeader$NC
  cbcBlock <- readBin(filPtr, double(), n = kntFloats, size = 4)
  iknt <- 1
  CBCterms[[iknt]] <- firstHeader$TEXT
  repeat {
    iknt <- iknt + 1
    thisHeader <- readCBCHeader(filPtr)
    # Don't read past EOF
    if (length(thisHeader) > 0) {
      if (thisHeader$TEXT == firstHeader$TEXT) {
        cbcBlock <- readBin(filPtr, double(), n = kntFloats, size = 4)
        break
      } else {
        CBCterms[[iknt]] <- thisHeader$TEXT
        cbcBlock <-
          readBin(filPtr, double(), n = kntFloats, size = 4)
      }
    }
    # Prevent Runaway
    if (iknt > 20) {
      break
    }
  }
  CBCTermSet <-
    list(firstHeader$TEXT,
         firstHeader$K,
         firstHeader$NR,
         firstHeader$NC,
         CBCterms)
  
  return(CBCTermSet)
}

#===============================================
# Function search for a defined budget term
# and returns a vector of values 
# by Stress Periods idenified in range of values
#===============================================
readCBCbinByTerm <- function(filPtr, term, SP_rng) {
  bigVector <- vector('numeric')
  HeaderRead <- readCBCHeader(filPtr)
  kntFloats <- HeaderRead$K * HeaderRead$NR * HeaderRead$NC
  cbcBlock <- readBin(filPtr, double(), n = kntFloats, size = 4)
  i <- 1
  cat(paste("0%.."))
  if (HeaderRead$TEXT == term  &&
      is.element(HeaderRead$KPER, SP_rng &&
                 thisHeader$KPER <= max(SP_rng))) {
#    print ("found on first try")
    bigVector <- c(bigVector , cbcBlock)
    i <- i + 1
  } else{
    repeat {
      thisHeader <- readCBCHeader(filPtr)
      # Don't read past EOF
      if (length(thisHeader) > 0) {
        if (thisHeader$TEXT == term  &&
            is.element(thisHeader$KPER, SP_rng) &&
            thisHeader$KPER <= max(SP_rng)) {
          i <- i + 1
#          print(paste(thisHeader$TEXT, thisHeader$KSTP, thisHeader$KPER))
          cbcBlock <-
            readBin(filPtr, double(), n = kntFloats, size = 4)
          bigVector <- c(bigVector, cbcBlock)
        } else {
          cbcBlock <- readBin(filPtr, double(), n = kntFloats, size = 4)
        }
      }
      # don't read everything unless necessary
      if (length(thisHeader) == 0) {
        cat('\n')
        break
      }
      
      if ( thisHeader$KPER > max(SP_rng)) {
        cat('\n')
        break
      }
      cat(paste('\r',format(as.numeric(thisHeader$KPER)/max(SP_rng)*100,digits=2,nsmall=2),"%"))
      
    }
  }
  
  return(bigVector)
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
CBCTermSet <- listBinHeaders(to.read)
ncols <- as.numeric(CBCTermSet[5])
nrows <- as.numeric(CBCTermSet[4])
nlays <- as.numeric(CBCTermSet[3])
CBCterms <- unlist(CBCTermSet[6])

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

print(paste("max value:",max(NetRchBySP)))

source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ModflowBinary/netRCH_anim.R")