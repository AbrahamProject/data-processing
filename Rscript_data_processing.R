################################################################################
#
# Abraham Moyal
# 16.02.23
# Data Processing
#
################################################################################

# Libraries

install.packages("devtools")
library(devtools)

################################################################################
# source of all used functions

source_url( "https://raw.githubusercontent.com/AbrahamProject/data-processing/main/data_processing.R" )


################################################################################
#-------------------------------------------------------------------------------
################################################################################

# CID data processing

# Input for all data processing


Directory <- "W:/users/Abraham/Exp 004/MS data processing/MS dial output/CID Alignment"
file <-"Area_CID15min.txt"
FC<-5
n_CID <- c(1,2,5,6,8,10)
fix <- 4
################################################################################

# Data import 

dataframe <- data_import(Directory,file)

data <-dataframe$data

standard <- dataframe$standard[n_CID,]

################################################################################

# Filtering for ms2 identified lipids

data <- filtering(data)

################################################################################

# Remove features based on blank samples

type="EP"

dataframe_EP <- blanking(data,type,FC)

type <- "SP"

dataframe_SP <- blanking(data,type,FC)

################################################################################

# Normalizaion using internal standard

type="EP"

dataframe_EP <- normalization(dataframe_EP,type,standard,fix)

type="SP"

dataframe_SP <- normalization(dataframe_SP,type,standard,fix)

################################################################################

#   OD correction using OD data
################################################################################

dataframe_EP <- OD_correction(dataframe_EP)

dataframe_SP <- OD_correction(dataframe_SP)

################################################################################

# CID data storage

CID <- list(dataframe_EP,
            dataframe_SP
)

names(CID) <-c("CID-EP-data","CID-SP-data")
################################################################################
#-------------------------------------------------------------------------------
################################################################################

# EAD data processing
#--------------------
# Input for all data processing


Directory <- "W:/users/Abraham/Exp 004/MS data processing/MS dial output/EAD  Alignment"
file <-"EAD_15min_all.txt"
FC<-5
n_EAD <- c(1,3)
fix <- 1
################################################################################

# Data import 

dataframe <- data_import(Directory,file)

data <-dataframe$data

standard <- dataframe$standard[n_EAD,]


################################################################################

# Filtering for ms2 identified lipids

data <- filtering(data)

################################################################################

# Remove features based on blank samples

type="EP"

dataframe_EP <- blanking(data,type,FC)

type <- "SP"

dataframe_SP <- blanking(data,type,FC)

################################################################################

# Normalizaion using internal standard

type="EP"

dataframe_EP <- normalization(dataframe_EP,type,standard,fix)

type="SP"

dataframe_SP <- normalization(dataframe_SP,type,standard,fix)

################################################################################

#   OD correction using OD data
################################################################################

dataframe_EP <- OD_correction(dataframe_EP)

dataframe_SP <- OD_correction(dataframe_SP)

################################################################################

# EAD data storage

EAD <- list(dataframe_EP,
            dataframe_SP
)

names(EAD) <-c("EAD-EP-data","EAD-SP-data")


################################################################################

# Data storage


DATA <- list(CID=CID,EAD=EAD)

################################################################################
