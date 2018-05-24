#####################################################################
# Step3_Temporal_Series_Generation
#
# Christopher Ruybal and Karl Oetjen 
# Created: 11/10/2017
# Last modified: 11/10/2017 
#
####################################################################

# Load Libraries
library(openxlsx)
library(gtools)
library(data.table)
library(zoo)
library(parallel)
library(snow)

# **** NOTE: All Excel Files to be used must be CLOSED and NOT OPENED ELSEWHERE **** ######

# Location were input files are located
indir = 'C:/Users/kaoetjen/Desktop/PhD/Paper 2 Homologous Series/Homologous Series Screening R/Filtered Homo Series/'
list.name = list.files('C:/Users/kaoetjen/Desktop/PhD/Paper 2 Homologous Series/Homologous Series Screening R/Filtered Homo Series/',pattern = 'Day ')
file.list = paste(indir,list.name,sep="")

# Sort file list as increasing days
file.list = mixedsort(file.list)

# Indicate order of "days" in file list
days <- c("1","2","6","8","10","12","14","16","18","21","25","31","45","64","87")

# Read each Excel file and name list according to day
df.list <- lapply(file.list, openxlsx::read.xlsx)
names(df.list) <- days

# Combine all tables into single table 
df <- rbindlist(df.list, idcol = "id")

# Set as data.frame
df <- data.frame(df)

# Remove some columns
df$X1 <- NULL
df[,13:(dim(df)[2])] <- NULL

# Set names 
names(df) <- c("Day","series","ID","Y/N","C1","C2","C3","C4","C5","C6","C7","C8")

# Populate series with correct values so that every row has a series value 
df$series <- na.locf(df$series)

data <- df 

dataK <- df
  
  # Run as cluster to decrease computational times 
  # START clustering. Set to 1 core less than computer has available. 
  nclus = detectCores() - 1
  clus <- c(rep("localhost", nclus))
  
  # set up cluster and data
  cl <- parallel::makeCluster(clus, type = "SOCK")
  parallel::clusterEvalQ(cl, c(library(openxlsx),library(gtools),library(data.table),library(zoo)))
  parallel::clusterExport(cl, list("dataK","data"))
  
  # split rows:
  Vec_row <- dim(dataK)[1]
  splt = rep(1:nclus, each = ceiling(Vec_row/nclus), length.out = Vec_row)
  newdlst = lapply(as.list(1:nclus), function(w) dataK[splt == w,])
  
  
  system.time(out.clus <- parLapply(cl, newdlst, function(lst){
    
    # get table dimensions 
    row = dim(lst)[1]
    col = dim(lst)[2]
  
    # preallocate table to store information 
    Unique_days <- data.frame(matrix(NA,nrow = dim(lst)[1],ncol = dim(lst)[2]+3))
    names(Unique_days) <- c("Days","series","ID","Y/N","C1","C2","C3","C4","C5","C6","C7","C8","MaxLength","Similiar to Day","ID of Day")
  
  for(i in 1:row){
    
    # create index for day of interest
    idx1 <- which(data$Day==lst$Day[i])
    
    # comparison of each row 
    res <- data.frame(matrix(0,nrow = dim(data)[1]-length(idx1),ncol=col))
    
    # For columsn 5 to end 
    for(j in 5:col){
      
      # Get appropriate row
      A <- lst[i,]
      # Get value to check
      B <- A[j]
      
      # create index for day of interest
      idx <- which(data$Day==A$Day)
      
      # Get remaining data table not equal to the day that is being checked 
      C <- data[-idx,]
      
      # temp matrix
      temp <- data.frame(matrix(NA,nrow = dim(C)[1],ncol=col))
      
      
      # calculate ppm for each row of C
      for (k in 1:(dim(C)[1])){
        for (l in 5:col){
          
          temp[k,l] <- abs(((B - C[k,l])/B)*1000000)}}
      
      # get indx of those lesss than 10 (ppm error)
      idx <- temp <= 10
      # Store
      res[idx] <- 1
    }
    
    # Now that we have all those similiar to series 1, let's find the larger of the 2
    # length original 
    Aclip <- A[,5:length(A)]
    lenO <- length(Aclip[!is.na(Aclip)])
    # length new
    lenN <- rowSums(res)
    # identify rows with more than 3 (length of series)  #**************************************######
    more <- which(lenN >= 3)
    # get rows in C
    lrows <- C[more,]
    lrows2 <- lrows[,5:length(lrows)]
    # get row names
    #lrowName <- row.names(lrows)
    
    # get row day and ID
    lrowDday <- lrows$Day
    lrowID <- lrows$ID
    

    temp2 <- data.frame(matrix(NA,nrow = length(more),ncol=1))
    # determine max length of new
    for (xx in 1:(dim(lrows2)[1])){
      
      temp2[xx,1] <- length(lrows2[xx,!is.na(lrows2[xx,])])
    }
    
    # which is larger row
    larg <- which.max(temp2[,1])
    # now combine orignal and largest of news
    new <- rbind(A,lrows[larg,])
    new2 <- new[,5:length(new)]
    
    temp2 <- data.frame(matrix(NA,nrow = dim(new2)[1],ncol=1))
    # determine max length of new
    for (xx in 1:(dim(new)[1])){
      
      temp2[xx,1] <- length(new2[xx,!is.na(new2[xx,])])
    }
    # which is larger between new and old 
    larg <- which.max(temp2[,1])
    
    ### IF rows are identical length, then choose row with smaller day
    
    if (anyDuplicated(temp2)>0){
      
      # See which day came first 
      smaller_day <- min(as.numeric(new[,1]))
      # Index for smaller day
      idxL <- which(new[,1]==as.character(smaller_day))  
      # store larger in unique
      Unique_days[i,1:col] <- new[idxL,]
      
    }else{
      
      # store larger in unique
      Unique_days[i,1:col] <- new[larg,]
      
    }
    
    if (is.na(new[larg,2])){
      uu <- unique(new[,2])
      if (all(is.na(unique(new[,2])))){
        'skip'
      }else{
        Unique_days[i,2] <- uu[!is.na(unique(new[,2]))]
      }
    }
    
    Unique_days[i,col+1] <- temp2[larg,]
    Unique_days[i,col+2] <- paste(c(A$Day,lrowDday), collapse=",") 
    Unique_days[i,col+3] <- paste(c(A$ID,lrowID), collapse=",") 
    
  }
    return(Unique_days)
    
    }
  ))
  
  # Stop cluster on cores 
  stopCluster(cl)
  
  # Combine predictions and variance from each core used 
  combined <- rbind(out.clus[[1]],out.clus[[2]],out.clus[[3]])
  
  # Finally get unique 
  col <- 12
  final_unique <- combined[!duplicated(combined[,1:col]),]
  

library(xlsx)
write.xlsx(final_unique,"C:/Users/kaoetjen/Desktop/PhD/Paper 2 Homologous Series/Homologous Series Screening R/Filtered Homo Series/Well S Similar Series Final Results.xlsx")


