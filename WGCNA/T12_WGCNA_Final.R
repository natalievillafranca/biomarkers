#### data input and cleaning ####
#isogroupDN59835_c1_g1 needed to be removed from T0 to match T1

rm(list=ls()) 
setwd('/Users/natalievillafranca/Desktop/biomarkers')
options(stringsAsFactors=FALSE)
library("conflicted")
library(tidyverse)
library(readxl)
library(dplyr)

#Reading in the table of counts per isogroup by sample
CountsHost=read.table("AllCountsNH_T12.tab",header=TRUE,row.names=1)

# figure out which samples didn't map to any isogroups
CountsHost_totals<-CountsHost

# summarise_all(sum)
iso.zero<-apply(CountsHost_totals,2,function(x) all(x==0))
iso.zero 
#there were none 

head(CountsHost)
length(CountsHost[,1]) # 36650 isogroups 

names(CountsHost)

# Read the data from CSV files
meta_data <- read.csv("t12_metadata.csv")
names(meta_data)

growth_data <- read.csv("AcerMorphologyData_Imputed2.csv")
names(growth_data)
IDKeyT0 <- read.csv("IDKeyT0.csv")

#creating a new column named sample similar to the one in meta_data
growth_data$sample <- paste0("G", growth_data$Genotype, "r", growth_data$FragID)
print(growth_data)
growth_data <- growth_data[, c("sample", names(growth_data)[-which(names(growth_data) == "sample")])]

#reading in normalized relative abundances from microbe data 
asv <- read.csv("TaxonRelAbunByGenusSample_ps_rare_9_5_24.csv")
asv <- asv %>% select(-c(1, 2))
asv <- t(asv)
new_colnames <- asv[1, ]
asv <- asv[-1, ]
colnames(asv) <- new_colnames
asv <- as.data.frame(asv)
asv$sample <- rownames(asv) 
asv <- asv %>%
  select(sample, everything())
rownames(asv) <- NULL

alpha <- read.csv("alpha_diversity_ps_rare.csv")
colnames(alpha)[colnames(alpha) == "X"] <- "sample"
alpha <- as.data.frame(alpha)

#merging data
merged_data <- merge(meta_data, growth_data, by = "sample", all = TRUE)
merged_data <- merge(merged_data, asv, by = "sample", all = TRUE)
merged_data <- merge(merged_data, alpha, by = "sample", all = TRUE)
merged_data <- merged_data[!is.na(merged_data$SRR), ]
names(merged_data)
names(CountsHost) <- gsub(".*_(SRR\\d+)_.*", "\\1", names(CountsHost))

#create a pairing of sample names from counts file with metadata from experimental design
merged_data <- merged_data[order(merged_data$SRR), ]
idx<-which(merged_data$SRR %in% colnames(CountsHost))
merged_data$SRR[idx]==names(CountsHost)
sampleID<-merged_data$SRR[idx]

####creating DESeq2 dataset - to log transform counts data for WGCNA based pipeline ####
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)

#Remove isogroups with low counts from dataset - count less than 10 in more than 90% of samples
CountsHost$low = apply(CountsHost[,1:187],1,function(x){sum(x<=10)})  #making new column counting number of samples with counts <=10 within each isogroup (host) 
colnames(CountsHost)
CountsHost$low
counts<-CountsHost[-which(CountsHost$low>168.3),] #168.3 is 90% of 187 samples - get rid of count less than 10 in more than 90% of samples

nrow(counts) #18774 genes pass filter => high expression isogroups 

countstrim<-counts[,c(1:187)]

#making sure all isogroups from t0 are in t12 before rlogging 
conflicts_prefer(base::setdiff)

countstrimT0 <-read.table("countstrim_T0.tab",header=TRUE,row.names=1) 

setdiff(rownames(countstrimT0), rownames(countstrim))

# Identify missing columns in datExprT12 that exist in datExpr
missing_rows <- setdiff(rownames(countstrimT0), rownames(countstrim))

missing_df <- as.data.frame(matrix(0, ncol = ncol(countstrim), nrow = length(missing_rows), 
                                   dimnames = list(missing_rows, colnames(countstrim))))

# Combine with the original dataframe
countstrim <- bind_rows(countstrim, missing_df)

#going to create a model matrix of identifying info -- in this case it would be genotype and site, this turns them
#into new columns and are identified as either 1 or 0 as either presence or absence 
genotype<-merged_data$genotype
site<-merged_data$Site
conditions=data.frame(cbind(genotype))
conditions$genotype<-as.factor(conditions$genotype)
head(conditions)
summary(conditions)

#now transform counts data using DESeq2, makes structure of data and conditions in a readable format for Deseq2
ddsCOUNTS<-DESeqDataSetFromMatrix(countData=countstrim,colData=conditions,design = ~genotype)

#save(ddsCOUNTS, file="ddsCOUNTS_T12.RData")
lnames=load(file="ddsCOUNTS_T12.RData")

#rlog transform, this takes a very long time, do in HPC with at least 248 gb and 64 cpus so it won't take too long. once you've done, you should
#be able to save as an RData file to just lnames into R, so you don't have to run it every single time. 
#rlogCOUNTS_T12<-rlog(ddsCOUNTS,blind=TRUE) #use blind=TRUE to not account for experimental design
#save(rlogCOUNTS_T12, file="rlogCOUNTS_T12.RData")

lnames = load(file="rlogCOUNTS_T12.RData") 


head(assay(rlogCOUNTS))
##Table of rlog transformed values for each gene for each sample

#make rlogCOUNTS a dataframe to be easily worked with 
dat=as.data.frame(assay(rlogCOUNTS))
colnames(dat)<-names(counts[1:187])
dat <- dat%>%select(sort(names(.)))
boxplot(dat)
#how do expression plots look overall? most genes are from 7-13 for most samples. 
#Host: about 8 have very high outliers (>20)

##### WGCNA Installation
#BiocManager::install("WGCNA", force = TRUE)
#install.packages("WGCNA")
#install.packages("Hmisc")
library(Hmisc)

#### Data input, cleaning and pre-processing
library(WGCNA)
disableWGCNAThreads()
options(stringsAsFactors=F)
head(dat)
dim(dat)

#create CSV file of dat sample names to compare to sample names to see where mismatches are
datnames <- names(dat)

#write.csv(datnames, file = "datnames.csv", row.names = F)

#transpose expression data
datExpr0=as.data.frame(t(dat[,1:187]))
#save(datExpr0, file="datExpr0.RData")

#check for genes with too many missing values
#badGenes = which(!gsg$goodGenes)
#badSamples = which(!gsg$goodSamples)

# Show bad genes and samples by name
#badGeneNames = colnames(datExpr0)[badGenes]
#badSampleNames = rownames(datExpr0)[badSamples]

# View them
#print(badGeneNames)

#removing isogroupDN59835_c1_g1 from further analysis
#datExpr0 = datExpr0[, colnames(datExpr0) != "isogroupDN59835_c1_g1"]
#left with 18774 isogroups now 


#Coding in the colony data
allTraits<-merged_data[idx,]
Test1<-allTraits
names(Test1)

for (i in 1:11){ # a loop to store all of the treatments as individual columns denoted 1/0 Y/N
  goi<-unique(Test1$Site)[i]
  tmp<-data.frame(ifelse(Test1$Site==goi,1,0))
  colnames(tmp)<-goi
  Test1<-cbind(Test1,tmp)
}

for (i in 1:10){ # a loop to store all of the genos as individual columns denoted 1/0 Y/N
  goi<-unique(Test1$genotype)[i]
  tmp<-data.frame(ifelse(Test1$genotype==goi,1,0))
  colnames(tmp)<-goi
  Test1<-cbind(Test1,tmp)
}

#Test1$GxT<-paste(Test1$genotype,"x",Test1$Site) #adding in all of the possible interactions terms
#for (i in 1:100){
#  goi<-unique(Test1$GxT)[i]
#  tmp<-data.frame(ifelse(Test1$GxT==goi,1,0))
#  colnames(tmp)<-goi
#  Test1<-cbind(Test1,tmp)
#}

library(dplyr)
names(Test1)

names(Test1)

#all traits without interactions or site since there was not significant site correlation with modules 
Test=Test1[,c(2,4, 5,13, 15:51, 58:107, 110:119)]

names(Test)


library(naniar)
Test$T3_Status <- replace(Test$T3_Status, Test$T3_Status == "A", 0)
Test$T3_Status <- replace(Test$T3_Status, Test$T3_Status == "D", 1)
Test$T3_Status <- replace(Test$T3_Status, Test$T3_Status == "M", NA)
Test$T3_Status<- as.numeric(Test$T3_Status)

Test$T6_Status <- replace(Test$T6_Status, Test$T6_Status == "A", 0)
Test$T6_Status <- replace(Test$T6_Status, Test$T6_Status == "D", 1)
Test$T6_Status <- replace(Test$T6_Status, Test$T6_Status == "M", NA)
Test$T6_Status<- as.numeric(Test$T6_Status)

Test$T9_Status <- replace(Test$T9_Status, Test$T9_Status == "A", 0)
Test$T9_Status <- replace(Test$T9_Status, Test$T9_Status == "D", 1)
Test$T9_Status <- replace(Test$T9_Status, Test$T9_Status == "M", NA)
Test$T9_Status<- as.numeric(Test$T9_Status)

Test$T12_Status <- replace(Test$T12_Status, Test$T12_Status == "A", 0)
Test$T12_Status <- replace(Test$T12_Status, Test$T12_Status == "D", 1)
Test$T12_Status <- replace(Test$T12_Status, Test$T12_Status == "M", NA)
Test$T12_Status<- as.numeric(Test$T12_Status)

summary(Test) #all columns are numeric, moving on 

allTraits<-Test
#allInteract<-Test2

# Form a data frame analogous to expression data that will hold the clinical traits
Samples = rownames(datExpr0);
rownames(allTraits) = allTraits$SRR;
traitRows = match(Samples, allTraits$SRR);
datTraits = allTraits[traitRows, -1];

#rownames(allInteract) = allInteract$SRR
#datInteract=allInteract[,-1]

table(rownames(datTraits)==rownames(datExpr0)) 
#should return TRUE if datasets align correctly, otherwise your names are out of order

#sample dendrogram and trait heat map showing outliers
#unsigned - values are either 0 or 1 
#signed - values are between -1 to 1
library(flashClust)
conflicts_prefer(WGCNA::cor)
A=adjacency(t(datExpr0),type="signed")                 #SELECT SIGNED OR UNSIGNED HERE
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")

library(dplyr)
# Convert traits to a color representation where red indicates high values
str(datTraits)
names(datTraits)
#making sure asvs are numeric not character 
datTraits[[2]] <- as.numeric(datTraits[[2]])
datTraits <- datTraits %>%
  mutate(across(46:75, as.numeric))
str(datTraits)
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
print(datTraits)
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlierC=outlierColor,traitColors)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")

# Remove outlying samples from expression and trait data 
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datExprT12=datExpr0[!remove.samples,]
datTraitsT12=datTraits[!remove.samples,]
#file


####FOR MODULE PRESERVATION, GETTING RID OF ISOGROUPS THAT DO NOT MATCH T0 DATTRAITS 
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Load the expression and trait data saved in the first part
#this lnames file was also run in july but forgot to change the name 
lnames = load(file = "BiomarkersTraitsAndInteractionsAssigned_Sep24.RData") 
#The variable lnames contains the names of loaded variables

# Load network data saved in the second part. Change this based on the level of merge you want to work with
lnames = load(file = "BiomarkersNetworkSymb_rlog_signed_unmerged_sft9_Sep24ASV.RData")

#lnames = load(file = "BiomarkersTraitsAndInteractionsAssigned_T12_March.RData")

#Changing 34 to 36 in genotype column 
datTraitsT12 <- datTraitsT12 %>%
  mutate(genotype = ifelse(genotype == 34, 36, genotype))
#changing 34 to 36 in column 36 
names(datTraitsT12)[names(datTraitsT12) == "34"] <- "36"

names(datTraitsT12)
# rename the traits to look nice for the heatmap, and reorder them so they are grouped by symb etc
datTraitsT12=datTraitsT12%>%relocate("50", .after = "T12_Vinter")
datTraitsT12=datTraitsT12%>%relocate("3", .after ="50")
datTraitsT12=datTraitsT12%>%relocate("1", .after ="3")
datTraitsT12=datTraitsT12%>%relocate("31", .after ="1")
datTraitsT12=datTraitsT12%>%relocate("44", .after ="31")
datTraitsT12=datTraitsT12%>%relocate("36", .after ="44")
datTraitsT12=datTraitsT12%>%relocate("7", .after ="36")
datTraitsT12=datTraitsT12%>%relocate("41", .after ="7")
datTraitsT12=datTraitsT12%>%relocate("62", .after ="41")
datTraitsT12=datTraitsT12%>%relocate("13", .after ="62")

resval_map <- c("50" = 1, "3" = 2, "1" = 3, "31" = 4, "44" = 5, 
                "36" = 6, "7" = 7, "41" = 8, "62" = 9, "13" = 10)

# Ensure genotype is character for matching
datTraitsT12$genotype <- as.character(datTraitsT12$genotype)

# Add new column using the mapping
datTraitsT12$ResValueRank <- resval_map[datTraitsT12$genotype]



resvalsum_map <- c("50" = 6197.946, "3" = 3669.834, "1" = 3078.816, "31" = 2944.329,
                   "7" = 2895.230, "36" = 2697.192, "44" = 2682.318, "62" = 2253.453,
                   "41" = 1623.347, "13" = 796.824)

# Assign values to new column based on genotype
datTraitsT12$ResValSum <- resvalsum_map[as.character(datTraitsT12$genotype)]

###added may 20 
###SETTING UP FOR MODULE PRESERVATION
#only keeping datTraitsT12 columns that match with datTraits columns (isogroups)
#have to combine and then separate because T12 has more isogroups (columns) and the number 
#of isogroups have to match those of datTraits
setdiff(colnames(datExpr), colnames(datExprT12))

datExprT12 <- datExprT12 %>% select(any_of(colnames(datExpr)))

setdiff(colnames(datExpr), colnames(datExprT12))


#save(datExprT12, datTraitsT12, file="BiomarkersTraitsAndInteractionsAssigned_T12_March.RData")
#save(datTraitsT12, file="T12_BiomarkersTraits.RData")

#move into module pres