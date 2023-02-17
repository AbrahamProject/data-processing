################################################################################
#
# Abraham Moyal
# 16.02.23
# Data Processing
# 1) Blanking
# 2) Filtering for ms2 identified lipids
# 2) Normalization
# 3) OD correction
################################################################################
#
# Load libraries

library(factoextra)
library(ggplot2)
library(dplyr)
library(ggVennDiagram)
library(tidyverse)
library(rcdk)
library(cluster)
library(svDialogs)
library(plotly)
library(rsconnect)
library(stringi)
library(insight)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")


################################################################################
#
# Print the functions and the parameters


title <- "This script contains 5 functions: "
function_0 <- "     1) Data import data_import(Directory,file)"
function_1 <- "     2) Blanking: blanking(dataframe,type,FC)"
function_2 <- "     3) Filtering: filtering(dataframe)"
function_3 <- "     4) Normalization: normalization(dataframe,type,standard)"
function_4 <- "     5) OD correction: OD_correction(dataframe)"

print_color(paste0(title,"\n",function_1,"\n",function_2,"\n",function_3,"\n",function_4),col="red")


title <-"     1) Data import data_import(Directory,file)"
aim <- "Aim: Import MS-dial alignment results"
function1 <-"Parameters: The function has 2 parameters a. Directory b. file"
a<-"  a. Directory: Location of the MS-dial output"
b<-"  b.file: file name of MS-dial output"
print_color(paste0(" ","\n",title,"\n",aim,"\n",function1,"\n",a,"\n",b),col="blue")


title <-"     2) Blanking: blanking(dataframe,type)"
aim <- "Aim: Remove sample features which are below 5 fold change of blanks"
function1 <-"Parameters: The function has 2 parameters a. dataframe b. type c. FC"
a<-"  a. dataframe: Alignment result from MS-dial full data frame"
b<-"  b. type: EP or SP depending on the phase"
c<-"  c. FC: Value of FC cuttoff"
print_color(paste0(" ","\n",title,"\n",aim,"\n",function1,"\n",a,"\n",b),col="blue")


title <-"     3) Filtering: filtering(dataframe)"
aim <- "Aim: Filter for features which correspond to ms2 identified lipids"
function1 <-"Parameters: The function has 2 parameters dataframe"
a<-"  dataframe: Alignment result from MS-dial full data frame"
print_color(paste0(" ","\n",title,"\n",aim,"\n",function1,"\n",a),col="blue")


title <-"     4) Normalization: normalization(dataframe,type,standard)"
aim <- "Aim: Normalize lipid classes with there corresponding internal standard"
function1 <-"Parameters: The function has 3 parameters a. dataframe b. type c. standard d. fix"
a<-"  a. dataframe: Alignment result from MS-dial full data frame"
b<-"  b.type: EP or SP depending on the phase"
c<-"  c. standard: data of the internal standards of all samples"
d<-"  d. fix: standard used for non present classes"
print_color(paste0(" ","\n",title,"\n",aim,"\n",function1,"\n",a,"\n",b,"\n",c,"\n",d),col="blue")



title <-"     5) OD correction: OD_correction(dataframe)"
aim <- "Aim: Normalize each sample by it's measured OD"
function1 <-"Parameters: The function has 1 parameters: dataframe "
a<-"  dataframe: Alignment result from MS-dial full data frame"
print_color(paste0(" ","\n",title,"\n",aim,"\n",function1,"\n",a),col="blue")




################################################################################
#
# Functions

# Data Import
#-------------------------------------------------------------------------------

data_import<-function(Directory,file){
  setwd(Directory)
  data<-read.delim(file)
  
  colname <- data[4,]
  
  data<- data[5:nrow(data),]
  colnames(data) <- colname
  
  standard_Index <- grep("(d7)",data$`Metabolite name`)
  standard_Index <- data$`Alignment ID`[standard_Index]
  
  index_d7<-match(standard_Index,data$`Alignment ID`)
  
  standard <- data[index_d7,]
  
  data <- data[-index_d7,]
  
  
  
  list<-list(data,
             standard)
  names(list)<-c("data","standard")
  return(list)
}


# Blanking
#-------------------------------------------------------------------------------

blanking <-function(dataframe,type,FC){
  seqEP<-c(138:142)
  seqSP<-c(144:147)
  if(type=="EP"){
    seq <- seqEP
  }else if(type=="SP"){
    seq <- seqSP
  }else{
    print("type, takes value 'EP' or 'SP' ")
  }
  
  data<-dataframe
  sequence <- seq(36,137,3)
  for(j in 1:nrow(dataframe)){
    mean_int <- mean(as.numeric(dataframe[j,seq]),na.rm=T)
    for (i in sequence){
      sample <- mean(as.numeric(dataframe[j,(i:(i+2))]),na.rm=T)
      
      if(sample<FC*mean_int){
        data[j,c(i,(i+1),(i+2))] <- c(0,0,0)
      }
      
    }
  }
  return(data)
}


# Filtering
#-------------------------------------------------------------------------------

filtering<-function(dataframe){
  data<-dataframe
  data_1 <- data[which(data$Ontology !="Unknown"),]
  data_2 <- data_1[which(data_1$Ontology !="null"),]
  data_3 <- data_2[which(data_2$Ontology !="Others"),]
 
  
  data_ms2 <- data_3[-grep("w/o",data_3$`Metabolite name`),]
  
  return(data_ms2)
}


# Normalization
#-------------------------------------------------------------------------------

normalization <- function(dataframe, type, standard, fix){
  if(type=="EP"){
   sequence<-c(36:81,85:86)
   x <- 0.04/0.12
  }else if(type=="SP"){
    sequence<-c(87:131,135:137) 
    x <- 0.04/0.48
  }else{
    print("type, takes value 'EP' or 'SP' ")
  }
  
  
  dataframe$standard<-NA
  standard$standard <- 1
  
  for(i in 1:nrow(standard)){
    name <- standard$Ontology[i]
    index <- which(dataframe$Ontology == name)
    
    if(length(index)==0){
            dataframe$standard[index] <- 1
            
            for(j in sequence){
              dataframe[index,j]<- (as.numeric(dataframe[index,j])/as.numeric(standard[i,j]) )*x
      }
    }
  }
  
  
  
  for(i in 1:nrow(dataframe)){
    if(is.na(dataframe$standard[i])==T){
      for(j in sequence ){
        dataframe[i,j] <- (as.numeric(dataframe[i,j])/as.numeric(standard[fix,j]) )*x
        
      }
    }
  }
  
 return(dataframe)  
}


# OD correction
#-------------------------------------------------------------------------------
OD_correction <- function(dataframe){
  setwd("W:/users/Abraham/Exp 004/od measurement")
  EP_OD <- read.csv("Exp_004_EP_harvest_summary_for_R.csv",header = F)
  EP<-cbind(EP_OD[,c(2,11,12,13,14,15,16,3,4,5,6,7,8,9,10)],c(1,1,1),EP_OD[,1])
  
  SP_OD <- read.csv("Exp_004_SP_2_harvest_summary_for_R.csv",header = F)/4
  SP<-cbind(SP_OD[,c(2,11,12,13,14,15,16,3,4,5,6,7,8,9,10)],c(1,1,1),SP_OD[,1])
  
  
  OD<-c(unlist(EP),unlist(SP))
  
  
  
  setwd(Directory)
  
  
  data <- dataframe
  index<-c(1:102)
  for(i in index){
    j <- i+35
    data[,j] <- as.numeric(data[,j])/as.numeric(OD[i])
    
  }      
  
  return(data)
  
}

