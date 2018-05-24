#####################################################################
# Step2_Filtering_out_crossmeshed_and_duplicate_series 
#
# Christopher Ruybal and Karl Oetjen 
# Created: 11/10/2017
# Last modified: 11/10/2017 
#
####################################################################


###################################################
#############DUPLICATE SCREENING###################
###################################################

# Load Library
library(openxlsx)

# Set directory
setwd('C:/Users/kaoetjen/Desktop/PhD/Paper 2 Homologous Series/Homologous Series Screening R')

# Load excel file
day1 <- openxlsx::read.xlsx('C:/Users/kaoetjen/Desktop/PhD/Paper 2 Homologous Series/Homologous Series Screening R/Day 1 Homo Series from PK.xlsx',sheet='Sheet1')

# Define a function 
# vec is a row (series)
fun <- function(vec)
{
  row = dim(vec)[1]
  col = dim(vec)[2]
  
  # place to store unique series (adjust the name "C7" to the correct vector length)
  series_unique <- data.frame(matrix(NA,nrow = dim(vec)[1],ncol = dim(vec)[2]+2))
  names(series_unique) <- c("series","ID","Y/N","C1","C2","C3","C4","C5","C6","C7","C8","MaxLength","Similiar to original rows...")
  
  
  for(i in 1:row){
    # comparison of each row temp matrix
    res <- data.frame(matrix(0,nrow = row-1,ncol=col))
    
    # 4:col is determined by data position
    for(j in 4:col){
    
    A <- vec[i,]
    B <- A[j]
    
    C <- vec[-i,]
    
    # temp matrix
    temp <- data.frame(matrix(NA,nrow = row-1,ncol=col))
   
    # 4:col is determined by data position
    # calculate ppm for each row of C
    for (k in 1:(row-1)){
      for (l in 4:col){
        
      temp[k,l] <- abs(((B - C[k,l])/B)*1000000)}}
      
    # get indx of those lesss than 10 (ppm error)
    idx <- temp <= 10
    # Store
    res[idx] <- 1
    }
    
    # Now that we have all those similiar to series 1, let's find the larger of the 2
    # length original 
    Aclip <- A[,4:length(A)]
    lenO <- length(Aclip[!is.na(Aclip)])
    # length new
    lenN <- rowSums(res)
    # identify rows with more than 3 (length of series)
    more <- which(lenN >= 3)
    # get rows in C
    lrows <- C[more,]
    lrows2 <- lrows[,4:length(lrows)]
    # get row names
    lrowName <- row.names(lrows)
    
    temp2 <- data.frame(matrix(NA,nrow = length(more),ncol=1))
    # determine max length of new
    for (xx in 1:(dim(lrows2)[1])){
      
      temp2[xx,1] <- length(lrows2[xx,!is.na(lrows2[xx,])])
    }
    
    # which is larger row
    larg <- which.max(temp2[,1])
    # now combine orignal and largest of news
    new <- rbind(A,lrows[larg,])
    new2 <- new[,4:length(new)]
    
    temp2 <- data.frame(matrix(NA,nrow = dim(new2)[1],ncol=1))
    
    # determine max length of new
    for (xx in 1:(dim(new)[1])){
      
      temp2[xx,1] <- length(new2[xx,!is.na(new2[xx,])])
    }
    
    # which is larger between new and old 
    larg <- which.max(temp2[,1])
    
    ### IF rows are identical length, then choose row with larger ID
    
    if (anyDuplicated(temp2)>0){
      
      # See which ID is larger
      larger_ID <- max(new[,2])
      # Index for larger ID
      idxL <- which(new[,2]==larger_ID)  
      # store larger in unique
      series_unique[i,1:col] <- new[idxL,]
      
    }else{
      
      # store larger in unique
      series_unique[i,1:col] <- new[larg,]
      
    }
  
    if (is.na(new[larg,1])){
      uu <- unique(new[,1])
      if (all(is.na(unique(new[,1])))){
        'skip'
      }else{
        series_unique[i,1] <- uu[!is.na(unique(new[,1]))]
      }
    }
    
    
    series_unique[i,col+1] <- temp2[larg,]
    series_unique[i,col+2] <- paste(c(toString(i),lrowName), collapse=",") 
    
  }

  # Finally get unique 
  final_unique <- series_unique[!duplicated(series_unique[,1:col]),]
  
  # identify if Y/N is Yes. User must specify this in the original file
  idx_yn <- which(final_unique$`Y/N` == "Yes")
  
  final_unique <- final_unique[idx_yn,]
  
}


# Run function created above 
dup_results <- fun(day1)



################################################
#############ADDUCT SCREENING###################
################################################

# Load libraries
library(nontarget)
library(nontargetData)
library(calibrate)
library(plotly)

# Source other scripts
source('C:/Nontarget_Data/Scripts/plothomol3.R')
source('C:/Nontarget_Data/Scripts/Interactive.R')

#Load adduct list
data("adducts")

#Load isotope list
data("isotopes")


## ***********Directory containing data to plot**********
indir <- 'C:/Nontarget_Data/Well S/Positive/'

##*********** Directory that you want to put your plots in*****************
outdir <- 'C:/Nontarget_Data/Well S/Positive/'


#*************load file you want to look at name it the sample list name*******
fname <- paste(indir,'Well S ACN Day1 DUP short.csv',sep='')

#************Name Run Varible - change this to sample name **************** 
peaklist<-read.csv(fname)

#### Load file with known peak names ####
namelist<-read.csv('C:/Nontarget_Data/SETAC/dummy.csv')


#Run adduct search 
adducts <- nontarget::adduct.search(peaklist, adducts, rttol=0.01, mztol = 1, ppm = TRUE, use_adducts = c("M+H","M+K","M+Na", "M+NH4"), ion_mode = "positive")


adducts2 <- as.data.frame(adducts$adducts)

adducts2 <- adducts2[adducts2$`group ID`!=0,]

adducts2$`group ID` <- as.integer(adducts2$`group ID`)

sortedadducts <- adducts2[order(adducts2$`group ID`),]

#create unique IDS for group#
ID_group <- unique(sortedadducts$`group ID`)


#for each unique group#

#for (ID in ID_group) {
#  
#  idx = which(adducts2$`group ID` == ID)
#  
#  #get group data
#  A = adducts2[idx,]
#  
#  #condense dup_results
#  dup_results2 <- dup_results[,4:10]
#  
#  for (dd in A$`m/z`){
#    
#    print(any(dup_results2==dd))
#    
#    print(which(dup_results2==dd))
#    
#  }
#}


# vec is going to be adducts2 and dup_results 
fun2 <- function(add1,dup1)
{
  row = dim(dup1)[1]
  col = dim(dup1)[2]+2
  
  # place to store unique series (adjust the name "C7" to the correct vector length)
  series_unique <- data.frame(matrix(NA,nrow = dim(dup1)[1],ncol = dim(dup1)[2]+2))
  names(series_unique) <- c("series","ID","Y/N","C1","C2","C3","C4","C5","C6","C7","C8","MaxLength","Duplicate rows","Adduct Rows","Adduct Type")
  

  # temp matrix to hold group ID's
  temp <- data.frame(matrix(NA,nrow = row,ncol=col-7))
  A <- dup1
  
  # matrix to hold adducts
  temp_adducts <- data.frame(matrix(NA,nrow = row,ncol=1))
  colnames(temp_adducts) <- "adducts(s)"
  
  # matrix to hold adducts
  temp_name <- data.frame(matrix(NA,nrow = row,ncol=1))
  colnames(temp_name) <- "ID"
  
  
  
  for(aa in 1:dim(add1)[1]){
    # each m/z
    C <- add1$`m/z`[aa]
    
    #creating temp matrix with group ID
    for(i in 1:row){
    for(j in 4:(col-4)){
      
      
      B <- A[i,j]
      
      
      # does adduct m/z exists in dup_results
      if (is.element(C,B)){
        
        temp[i,j-3] <- add1$`group ID`[aa]
        temp_adducts[i,1] <- add1$`adduct(s)`[aa]
        temp_name[i,1] <- A$ID[i]
        
      }else{
        
        'skip'
      }
      
      }
    }
  }
    
    
  test<- temp
  test2 <- temp
  test2$sameAs <- NA
  test2$longest <- NA
  test2$length <- NA
  
  for (i in 1:dim(test)[1]){
    
    A <- test[i,]
    
    for (j in 1:dim(test)[1]){
      
      B <- test[j,]
      
      if (i==j){
        'skip'
      }else{
        
        tst <- c(unique(A[!is.na(A)]),unique(B[!is.na(B)]))
        tst <- tst[duplicated(tst)]
        l <- length(tst)
        
        #### set  number of dupclicates allowed (l<2) #############################################
        if (l<2){
          'skip'
        }else{
          
          if (is.na(test2$sameAs[i])){
            
            test2$sameAs[i] <- paste(i,j,sep=',')
            
          }else{
            
            test2$sameAs[i] <- paste(test2$sameAs[i],j,sep=',')
            
          }
        }
      }
    }
    
    if (is.na(test2$sameAs[i])){
      'skip'
      
    }else{
    idx <- as.numeric(strsplit(test2$sameAs[i],split=",",fixed=TRUE)[[1]])
    idx <- sort(idx)
    
    temp2 <- matrix(NA,nrow=1,ncol=2)
    for (n in idx){
      
      D <- test[n,]
      f <- length(unique(D[!is.na(D)]))
      
      if (is.na(temp2[1])){
        
        temp2[1,1] <- f
        temp2[1,2] <- n
        
      }else if(temp2[1,1]==f){
        'skip'
        
      }else{
        
        val <- temp2[1,1]
        
        if (val>f){
          'skip'
        }else{
          
          temp2[1,1] <- f
          temp2[1,2] <- n
          
        }
      }
    }
    }
    
    if (is.na(test2$sameAs[i])){
      'skip'
      
    }else{
    
    test2$longest[i] <- temp2[1,2]
    
    test2$length[i] <- temp2[1,1]
    
    rm(temp2)
    
    }
  }
  
  # Combine test2 with the identified adducts
  test2 <- cbind(temp_name,test2,temp_adducts)
  
  # for each row of test2, check if "sameAs" is NA. IF it is, then rename "adduct(s)" to NA
  for (ii in 1:dim(test2)[1]){
    
    # get value
    same <- test2$sameAs[ii]
    
    if (is.na(same)){
      
      test2$`adducts(s)`[ii] <- NA
      
    }else{
      'skip'
    }
  }
  

  #create column with IDs
  test2$AdductOfID <- NA
  
 # for each row in test2
  for (rr in 1:dim(test2)[1]){
    
    # get sameAs
    ss <- test2$sameAs[rr]
    
    if (is.na(ss)){
      'skip'
    }else{
      
      idx <- as.numeric(strsplit(ss,split=",",fixed=TRUE)[[1]])
      
      test2$AdductOfID[rr] <- paste(dup1$ID[idx],collapse=',')
      
    }
  }
 

test2 <- cbind(test2[,1:12],test2[,14],test2[,13])
colnames(test2) <- c("ID","X1","X2","X3","X4","X5","X6","X7","X8","sameAs","longest","length","AdductOfID","adducts(s)")

#Combine dup_results with test2
test3 <- cbind(dup_results,test2[,10:14])

  
for(jj in 1:dim(test3)[1]){
  
  if(is.na(test3$longest[jj])){
    'skip'
  }else{
  
  #get longest row
  longR <- test3$longest[jj]
  
  #get ID of that row
  ID_row <- test3$ID[longR]
  
  #reference of AdductIDs
  ID_add <- test3$AdductOfID[longR]
  
  #get ids not equal to ID of longest row
  ID_add <- as.character(ID_add)
  
  #split string
  idx <- as.numeric(strsplit(ID_add,split=",",fixed=TRUE)[[1]])
  
  # get ID not equal to longest row
  idN <- which(ID_row != idx)
  
  # get id(s) to delete
  ID_delete <- idx[idN]
  
  #now, delete corresponding ID
  #index IDS
  idx <- which(is.element(test3$ID,ID_delete))
    
  #delete
  #test3 <- test3[-idx,]
    
  test3$ID[idx] <- NA
  
  rm(idx)
  }
}


# find IDs with NA

idx <- which(is.na(test3$ID))

if(length(idx)==0){
  'skip'
}else{
test3 <- test3[-idx,]
}

return(test3)

}


final_res <- fun2(adducts2,dup_results)

library(xlsx)
write.xlsx(final_res,"C:/Users/kaoetjen/Desktop/PhD/Paper 2 Homologous Series/Homologous Series Screening R/Filtered Homo Series/Well S Homo Filtered Day 1.xlsx")
