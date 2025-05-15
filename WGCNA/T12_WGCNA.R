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
#lnames=load(file="ddsCOUNTS_T12.RData")

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

