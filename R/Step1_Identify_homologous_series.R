##################################################################################
# Step1_Identify_homologous_series
#
# Karl Oetjen and Christopher Ruybal
# Last modified: 5/10/2018
#
#################################################################################

# Non-Target Screening Tool
#Package Source "nontarget" developed by Martin Loos <Martin.Loos@eawag.ch>

# Load libraries 
library(nontarget)
library(nontargetData)
library(calibrate)
library(openxlsx)
library(plotly)

# Source additional scripts 
source('C:/Nontarget_Data/Scripts/plothomol3.R')
source('C:/Nontarget_Data/Scripts/Interactive Update 51517.R')


###################################
##### CREATING DATA PATH ##########
###################################

#Load adduct list
data("adducts")

#Load isotope list
data("isotopes")



# ***********Directory containing data to plot**********
indir <- 'C:/Nontarget_Data/Fluro/Simon/Krista POS GW/'

#*********** Directory that you want to put your plots in*****************
outdir <- 'C:/Nontarget_Data/Fluro/Simon/Krista POS GW/Results/'


#*************load file you want to look at name it the sample list name*******
fname <- paste(indir,'070715 Pos Jacksonville  Short.csv',sep='')

#************Name Run Varible - change this to sample name **************** 
peaklist<-read.csv(fname)

#### Load file with known peak names ####
namelist<-read.csv('C:/Nontarget_Data/SETAC/dummy.csv')




####################################
##### 1 MAKE ISOTOPE LIST ##########
####################################

#make isotope list that you want to look at you change canhe this to look at the isotope list [ View(isotopes) ]
iso<-make.isos(isotopes, use_isotopes=c("13C","15N","34S","37Cl","81Br","41K","13C","15N","34S","37Cl","81Br","41K","18O","2H"), use_charges=c(1,1,1,1,1,1,2,2,2,2,2,2,2,2))


###################################
##### 2 RUN PATTERN SEARCH ########
###################################

#Run Pattern Search
pattern<-pattern.search(peaklist,iso,cutint=10000,rttol=c(-0.05,0.05),mztol=5,mzfrac=0.1,ppm=TRUE,inttol=0.2,rules=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),deter=FALSE,entry=50)

# Plot isotope and masss defect 
plotisotopes(pattern);
plotdefect(pattern,elements=c("C"))


###########################################
####### 3 RUNNING ADDUCT SEARCH#############
###########################################

#Run adduct search 
adducts <- nontarget::adduct.search(peaklist, adducts, rttol=0.5, mztol = 5, ppm = TRUE, use_adducts = c("M+H","M+K","M+Na", "M+NH4"), ion_mode = "positive")

#Plot adduct
plotadduct(adducts)

#To show single pattern group and its relation to adduct groups######

plotall(pattern,adducts);
plotgroup(pattern,adducts,groupID=1,massrange=10,allmass=FALSE);


############################################
########### 4 HOMOLOGUE SEARCH #############
############################################

# Screen for homologue series
homol <- homol.search(peaklist,isotopes,elements=c("C","H","O"), use_C=TRUE, minmz=5, maxmz=1200, minrt=0.5, maxrt=45, ppm=TRUE, mztol=5, rttol=0.05, minlength=3, mzfilter=FALSE,vec_size=1E6)

#Plot results 
plothomol(homol,xlim=FALSE,ylim=FALSE,plotlegend=TRUE);

# plot series, m/z / mass defect
plothomol(homol,plotdefect=TRUE);

# If Peak Name is Known (PeakID,Name,Group ID)
plothomol3(homol,xlim=FALSE,ylim=FALSE,plotlegend=TRUE);


datAll <- data.frame(mz = homol[[1]][,1],
                     RT = homol[[1]][,3],
                     PeakID = homol[[1]][,4],
                     GroupID = homol[[1]][,5])

m = 0
KNO <- matrix(NA, nrow = length(namelist$Name[!is.na(namelist$Name)]), ncol = 3)

cnam <- c('PeakID','Name','GroupID')
colnames(KNO) = cnam

for (rrr in 1:length(datAll$mz)){
  if (!is.na(namelist$Name[rrr])) {
    m = m + 1
    KNO[m,1] = toString(datAll$PeakID[rrr])
    KNO[m,2] = toString(namelist$Name[rrr])
    KNO[m,3] = toString(datAll$GroupID[rrr])
    
  }
}



## setup a workbook with 3 worksheets

wb <- createWorkbook("Karl")

addWorksheet(wb = wb, sheetName = "Single Peaks", gridLines = FALSE)
writeDataTable(wb = wb, sheet = 1, x = homol[[1]][,1:5])

addWorksheet(wb = wb, sheetName = "Series", gridLines = TRUE)
writeData(wb = wb, sheet = 2, x = homol[[3]][,c(1,3,4,2)])

addWorksheet(wb = wb, sheetName = "Known Peaks", gridLines = TRUE)
writeData(wb = wb, sheet = 3, x = KNO)

saveWorkbook(wb,'C:/Nontarget_Data/Fluro/Simon/Krista POS GW/Results/070715 Pos Jacksonville  Short.xlsx',overwrite = TRUE)

# 1. spreadsheet defines membership of single peaks in a series
# write.csv(homol[[1]][,1:5],file.path(outdir,"results_George_CHO.csv"), row.names=FALSE)
# columns 1-4 = peaklist & peak IDs
# column 5 = ID of (homologues) series

# 2. spreadsheet defines series and the peaks therein
# write.csv(homol[[3]][,c(1,3,4,2)],file.path(outdir,"results2_Geogre_CHO.csv"), row.names=FALSE)
# column 1 = series ID
# column 2 = mean m/z increment in a series
# column 3 = mean RT increment in a series
# column 4 = IDs of peaks in the series
# all other columns: series properties


#############################################
######## 5 Combine Results ##################
#############################################

# Combine grouping results to components 

comp<-combine(pattern,adducts,homol,dont=FALSE, rules=c(TRUE,FALSE,FALSE));comp[[7]];

plotisotopes(comp);

plotcomp(comp,compoID=1,peakID=FALSE);


################################################
##### Plot interactive Version (html) ##########
################################################

plot_interact(homol)
