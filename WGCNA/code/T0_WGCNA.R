#### data input and cleaning ####
rm(list=ls()) 
options(stringsAsFactors=FALSE)
library("conflicted")
library(tidyverse)
library(readxl)
library(dplyr)

#Reading in the table of counts per isogroup by sample
CountsHost=read.table("AllCountsNH.tab",header=TRUE,row.names=1)

# figure out which samples didn't map to any isogroups
CountsHost_totals<-CountsHost

# summarise_all(sum)
iso.zero<-apply(CountsHost_totals,2,function(x) all(x==0))
iso.zero 
#there were none 

head(CountsHost)
length(CountsHost[,1]) # 33650 isogroups 
names(CountsHost)

# Read the data from CSV files
meta_data <- read.csv("meta.csv")
meta_data <- mutate(meta_data, sample = gsub("G34", "G36", sample))
names(meta_data)
meta_data <- meta_data[-c(276:278), ]

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

#create a pairing of sample names from counts file with metadata from experimental design
idx<-which(merged_data$SRR %in% colnames(CountsHost))
merged_data$SRR[idx]==names(CountsHost)
sampleID<-merged_data$SRR[idx]

####creating DESeq2 dataset - to log transform counts data for WGCNA based pipeline ####
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)

#Remove isogroups with low counts from dataset - count less than 10 in more than 90% of samples
CountsHost$low = apply(CountsHost[,1:275],1,function(x){sum(x<=10)})  #making new column counting number of samples with counts <=10 within each isogroup (host) 
colnames(CountsHost)
CountsHost$low
counts<-CountsHost[-which(CountsHost$low>247.5),] #247.5 is 90% of 275 samples - get rid of count less than 10 in more than 90% of samples

nrow(counts) #10785 genes pass filter => high expression isogroups ## new ref: 10785 genes pass filter

countstrim<-counts[,c(1:275)]

#write.table(countstrim, file = "countstrim_T0.tab")

#going to create a model matrix of identifying info -- in this case it would be genotype and site, this turns them
#into new columns and are identified as either 1 or 0 as either presence or absence 
genotype<-merged_data$genotype
site<-merged_data$Site
conditions=data.frame(cbind(genotype))
conditions$genotype<-as.factor(conditions$genotype)
head(conditions)
summary(conditions)

#now transform counts data using DESeq2, makes structure of data and conditions in a readable format for Deseq2
#ddsCOUNTS<-DESeqDataSetFromMatrix(countData=countstrim,colData=conditions,design = ~genotype)

#save(ddsCOUNTS, file="ddsCOUNTS.RData")

#rlog transform, this takes a very long time, do in HPC with at least 248 gb and 64 cpus so it won't take too long. once you've done, you should
#be able to save as an RData file to just lnames into R, so you don't have to run it every single time. 
#rlogCOUNTS<-rlog(ddsCOUNTS,blind=TRUE) #use blind=TRUE to not account for experimental design

lnames = load(file="rlogCOUNTS.RData") 

head(assay(rlogCOUNTS))
##Table of rlog transformed values for each gene for each sample

#make rlogCOUNTS a dataframe to be easily worked with 
dat=as.data.frame(assay(rlogCOUNTS))
colnames(dat)<-names(counts[1:275])
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
datExpr0=as.data.frame(t(dat[,1:275]))
#save(datExpr0, file="datExpr0.RData")

#check for genes with too many missing values
gsg=goodSamplesGenes(datExpr0, verbose=3)
gsg$allOK #TRUE

#Coding in the colony data
allTraits<-merged_data[idx,]
Test1<-allTraits
names(Test1)

for (i in 1:10){ # a loop to store all of the treatments as individual columns denoted 1/0 Y/N
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
Test=Test1[,c(2,4, 5,13, 15:51, 58:102, 104:118)]

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
datExpr=datExpr0[!remove.samples,]
datTraits=datTraits[!remove.samples,]
#this brings datTraits down to 268 from 275
#file
save(datExpr, datTraits, file="BiomarkersTraitsAndInteractionsAssigned_Sep24.RData")
save(datTraits, file="BiomarkersTraits_March24.RData")


################Moving on!  Network construction and module detection (tutorial PDF Step 2)
#####################################################################

library(WGCNA)
options(stringsAsFactors = FALSE)
library(flashClust)
disableWGCNAThreads()

#host
lnames = load(file="BiomarkersTraitsAndInteractionsAssigned_Sep24.RData") 
lnames

#Figure out proper SFT 
# # Choose a set of soft-thresholding powers
powers = c(seq(1,18,by=1)); #may need to adjust these power values to hone in on proper sft value
# # Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType="signed", verbose = 2) 
#want smallest value, to plateau closest to 0.9 (but still under)

# # Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of .9
abline(h=0.9,col="red")
# # Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

########################## going with 9, right below the line 
###Calculate the adjacencies with the soft thresholding power
softPower=9
adjacency=adjacency(datExpr, power=softPower,type="signed") #must change method type here too!!

#translate the adjacency into topological overlap matrix (TOM) and calculate the corresponding dissimilarity: 
#(to minimize effects of noise and spurious associations)
##Topological Overlap Matrix (TOM): TOM is a measure used in network analysis to assess the 
##similarity between nodes (genes) based on their shared neighbors. It quantifies the extent to which two nodes 
##have similar network connections. Higher TOM values indicate that two nodes are more closely connected in the network.
TOM=TOMsimilarity(adjacency,TOMType = "signed")

#Dissimilarity TOM (dissTOM): This is the complement of the TOM, representing the dissimilarity between nodes. 
#It is calculated as 1 - TOM.
##used in clustering and module detection processes within WGCNA. 
##Since higher TOM values indicate more similarity, dissTOM provides a measure of dissimilarity that 
##can be used with clustering algorithms (such as hierarchical clustering) to group genes into modules 
##based on their network dissimilarity.
dissTOM= 1-TOM

save(adjacency, TOM, dissTOM, file="BiomarkersSamplesAndTraitsTOM_signed_Sep24ASV.RData")

######################open this file instead if adjacency and TOM run on supercomputer, otherwise skip to geneTree:
lnames = load(file="BiomarkersSamplesAndTraitsTOM_signed_Sep24ASV.RData")
lnames
#######################################

#####Use hierarchical clustering to produce a tree of genes (based on what exactly? adjacency and similarity/topology)
geneTree= flashClust(as.dist(dissTOM), method="average")
par(mar=c(1,1,1,1))
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
###each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes

###Cutting up the tree into highly co-expressed gene modules (by branches)
minModuleSize=35 #we only want large modules, this is considered relatively high
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods) #lists the modules and how many genes are in each one

### HOST
#1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   
#1983 1544 1450 1358 1351  934  897  447  183  173  125  115  112   66   47 

##Plot module assignments under the gene tree
dynamicColors= labels2colors(dynamicMods)
#plot dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors sft=9, min Mod size=35")

#Merge modules whose expression profiles are very similar or choose not to merge
#calculate eigengenes
MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = 9)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigengenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")
#plot
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#start with a 95% similarity merge to get an initial sense of module trait correlations
#MEDissThres= 0 # 
#MEDissThres= 0.10
MEDissThres=0.18
#MEDissThres=0.2

##the lower the height cut, the more the modules have to be similar (more conservative)
#plot the cut line into tree
abline(h=MEDissThres, col="red")

merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
#merge module colors and find new eigengenes on the new merged modules
mergedColors= merge$colors
mergedMEs= merge$newMEs
length(unique(mergedColors)) #shows how many modules there are now, 
#there are 13 

##plot new module colors on gene tree under previous colors to see how they change
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors= mergedColors
#create numerical lables corresponding to the colors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs #a dataframe of module Eigengenes for each module for each sample

#save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file= "BiomarkersNetworkSymb_rlog_signed_unmerged_sft9_Sep24ASV.RData")


#################Relating modules to traits and finding important genes
#######################################################################
#can start here if you restarted R session

library(WGCNA)
library(stringr)
library(flashClust)
library(dplyr)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Load the expression and trait data saved in the first part
#this lnames file was also run in july but forgot to change the name 
lnames = load(file = "BiomarkersTraitsAndInteractionsAssigned_Sep24.RData") 
#The variable lnames contains the names of loaded variables

# Load network data saved in the second part. Change this based on the level of merge you want to work with
lnames = load(file = "BiomarkersNetworkSymb_rlog_signed_unmerged_sft9_Sep24ASV.RData")

lnames = load(file = "BiomarkersTraitsAndInteractionsAssigned_T12_March.RData")


#######################Replot module dendrogram

MEList= moduleEigengenes(datExpr, moduleColors, softPower=9)$eigengenes
MEs<- MEList

#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
abline(h=.20,col="blue")
abline(h=.15,col="green")

#now for module trait heatmap
#correlate eigengenes with external clinical traits to look for most significant associations
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

names(datTraits)

#Changing 34 to 36 in genotype column 
datTraits <- datTraits %>%
  mutate(genotype = ifelse(genotype == 34, 36, genotype))
#changing 34 to 36 in column 36 
names(datTraits)[names(datTraits) == "34"] <- "36"

names(datTraits)
# rename the traits to look nice for the heatmap, and reorder them so they are grouped by symb etc
datTraits=datTraits%>%relocate("50", .after = "T12_Vinter")
datTraits=datTraits%>%relocate("3", .after ="50")
datTraits=datTraits%>%relocate("1", .after ="3")
datTraits=datTraits%>%relocate("31", .after ="1")
datTraits=datTraits%>%relocate("44", .after ="31")
datTraits=datTraits%>%relocate("36", .after ="44")
datTraits=datTraits%>%relocate("7", .after ="36")
datTraits=datTraits%>%relocate("41", .after ="7")
datTraits=datTraits%>%relocate("62", .after ="41")
datTraits=datTraits%>%relocate("13", .after ="62")

resval_map <- c("50" = 1, "3" = 2, "1" = 3, "31" = 4, "44" = 5, 
                "36" = 6, "7" = 7, "41" = 8, "62" = 9, "13" = 10)

# Ensure genotype is character for matching
datTraits$genotype <- as.character(datTraits$genotype)

# Add new column using the mapping
datTraits$ResValueRank <- resval_map[datTraits$genotype]


resvalsum_map <- c("50" = 6197.946, "3" = 3669.834, "1" = 3078.816, "31" = 2944.329,
                   "7" = 2895.230, "36" = 2697.192, "44" = 2682.318, "62" = 2253.453,
                   "41" = 1623.347, "13" = 796.824)

# Assign values to new column based on genotype
datTraits$ResValSum <- resvalsum_map[as.character(datTraits$genotype)]
                                     
  
                                     
#Status, TLE, and SA for all time points 
#datTraits1=datTraits[,c(7:9,11,12,16:18,20,21,25:27,29,30,34:36,38,39, 44:53, 87, 59, 76, 77)]

datTraits2=datTraits[,c(8, 9,15, 16, 22, 23, 29, 30, 42:55, 61, 67, 70, 78, 80, 89, 92:100,102)]
#head(datTraits1)



# Recalculate MEs with color labels (MEs=module eigengenes)
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=9)$eigengenes 
MEs = orderMEs(MEs0)

#only keeping purple and tan
#MEs <- MEs[, c("MEpurple", "MEtan", "MEgreenyellow", "MEpink")]

#correlations of traits with eigengenes
# for overall heatmap
moduleTraitCor2 = cor(MEs, datTraits2, use = "p"); #p=pearsons
moduleTraitPvalue = corPvalueStudent(moduleTraitCor2, nSamples);
Colors=sub("ME","",names(MEs)) #just takes off "MEs" in front of the color names

#correlations of genes with eigengenes 
##used in next step of correlation process, will get done later in tutorial (can be down now or later)
moduleGeneCor=cor(MEs,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#represent module trait correlations as a heatmap
# module-trait correlations
library(RColorBrewer)
modLabels=sub("ME","",names(MEs))

ps=signif(moduleTraitPvalue,1)
cors=signif(moduleTraitCor2,2)
textMatrix = ps;

#displays only significant p values
textMatrix[ps>0.05]="-"
dim(textMatrix) = dim(moduleTraitCor2)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor2, 2),  "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "")

# Display the correlation values within a heatmap plot - all modules and all traits
par(mar = c(4, 5.5, 1, 1));
labeledHeatmap(Matrix = moduleTraitCor2,
               xLabels = names(datTraits2),
               ySymbols = modLabels,
               yLabels = modLabels,
               colorLabels = FALSE,
               colors = colorRampPalette(c("blue","lightblue","white","coral","red"))(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               cex.lab = 0.7,
               cex.lab.x = 0.5,  # Adjusts x-axis text size
               zlim = c(-0.7,0.7))
               
#purple module seems the most interesting, especially with regard to genotype and survival 
# module size barplot
labelShift=300 # increase to move module size labels to the right
par(mar = c(6, 8.5, 3, 3));
mct=table(moduleColors)
mct[modLabels]
x=barplot(mct[rev(modLabels)],horiz=T,las=1,xlim=c(0,4500),col=rev(modLabels))
text(mct[rev(modLabels)]+labelShift,y=x,mct[rev(modLabels)],cex=0.9) 



#selecting purple and tan
MEs2 = MEs%>%select("MEpurple", "MEtan","MEblack","MEbrown","MEturquoise","MEsalmon","MEgreenyellow","MEpink",
                    "MEcyan","MEmagenta","MEmidnightblue","MEblue","MEgreen")
# for pretty fig 
moduleTraitCor = cor(MEs2, datTraits2, use = "p"); #p=pearsons
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs2)) #just takes off "MEs" in front of the color names

#correlations of genes with eigengenes 
##used in next step of correlation process, will get done later in tutorial (can be down now or later)
moduleGeneCor2=cor(MEs2,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor2, nSamples)


#module membership vs gene significance plot -- part of module validation process
#how strong are the correlations? are there specific isogroups in here that might be worth looking 
#at in GO analysis? 
#Gene relationship to trait and important modules:
#############Treatment
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
custom_order <- c(50, 3, 1, 31, 44, 36, 7, 41, 62, 13)
datTraits$genotype <- factor(datTraits$genotype, levels = custom_order)
datTraits <- datTraits[order(datTraits$genotype), ]

weight = as.data.frame(datTraits2$"ResValueRank") 

#change to your trait name of interest
names(weight) = "ResValRank"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p")); 
#finds pearson correlations of GE and module eigengenes
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
#making dataframe of pvalues for gene module membership values in each module

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p")); 
# data frame of correlations between expression and trait of interest using pearson correlations
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)); 
#p-values for correlations between expression and trait of interest
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
#geneTraitSignificance is a dataframe of correlation/covariance of GE data and stage value
#GSPvalue is a list of correlation of gTS pvalue

#modules of interest for trait
moduleCols=c( "purple")

#plot scatter plots of gene significance vs module membership for all of these modules of interest
#add correlation and p-value, use this to look at how strong the modules are.
par(mfrow=c(1,1))
par(mar = c(2, 2, 2, 2));
par(bg = "white");
for (module in moduleCols) {
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("ModMem", module),
                     ylab = "Gene Sig for Trait",
                     main = paste("MM vs. GS\n"),
                     cex.main = 1, cex.lab = 1, cex.axis = 1.2, col = module)
} 

#this is the same as above but adds the isogroup name to each dot so you can pull out individual ones 
#that look interesting 
for (module in moduleCols) {
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("ModMem", module),
                     ylab = "Gene Sig for Trait",
                     main = paste("MM vs. GS\n"),
                     cex.main = 1, cex.lab = 1, cex.axis = 1.2, col = module)
  text(abs(geneModuleMembership[moduleGenes, column]),
       abs(geneTraitSignificance[moduleGenes, 1]),
       labels = rownames(geneModuleMembership)[moduleGenes], pos = 3, cex = 0.6, col = "black")
}
library(ggplot2)

###Visualization of Gene Networks
#figure 1 
############################################
#heatmaps of module expression with bar plot of eigengene
#use this to look at the different samples and make sure the modules make sense
#we want blocks by treatment, no single sample driving the differences.

#start with modules that looked strongest: #purple looks strongest for T3_Status, 
#green, pink, and black are best for T3_TLE
library(tibble)
#sort ME by traits 
MEs$sort=rownames(MEs)
datTraits <- datTraits %>%
  rownames_to_column("SRR")

# Joining using full_join
MEs <- full_join(MEs, datTraits, by = c("sort" = "SRR"))
print(MEs)

MEs <- MEs %>%
  column_to_rownames("sort")

#if you want to order by genotype ranking 
survival_genotype_order <- c(36, 1, 50, 3, 44, 7, 31, 13, 62, 41)
custom_order <- c(50, 3, 1, 31, 44, 36, 7, 41, 62, 13)
MEs$genotype <- factor(MEs$genotype, levels = custom_order)
MEs.sorted <- MEs[order(MEs$genotype), ]

#MEs.sorted <- na.omit(MEs.sorted[, "MD3-55", drop = FALSE])

#for microbes
#scatterplot_data <- data.frame(MEtan = MEs.sorted$MEtan, HIMB11 = MEs.sorted$HIMB11)
#ggplot(scatterplot_data, aes(x=HIMB11, y=MEtan)) +
#  geom_point()

#visualizing module eigengene  expression with boxplots 
# Create a data frame for boxplot
library(ggplot2)
library(tidyverse)
library(egg)
library(multcompView)
#create boxplot using ggplot2

##stats for purple module and genotype 
anova <- aov(MEpurple ~ genotype, data = MEs.sorted)
summary(anova)
#Df Sum Sq Mean Sq F value Pr(>F)    
#genotype      9 0.7718 0.08576   96.97 <2e-16 ***
#  Residuals   258 0.2282 0.00088                   
tukey <- TukeyHSD(anova)
tukey_pvals <- tukey$genotype[, "p adj"]

# Generate letters
letters <- multcompLetters(tukey_pvals)$Letters

# Convert to dataframe
letters.df <- tibble(genotype = names(letters), letter = letters)

letters.df <- letters.df %>%
  mutate(genotype = factor(genotype, levels = custom_order)) %>%
  arrange(genotype)

custom_order <- c(50, 3, 1, 31, 44, 36, 7, 41, 62, 13)
MEs$genotype <- factor(MEs$genotype, levels = custom_order)
MEs.sorted <- MEs[order(MEs$genotype), ]

which.module= "purple" #pick module of interest

ME=MEs.sorted[, paste("ME",which.module, sep="")]
datExpr_reordered <- datExpr[match(rownames(MEs.sorted), rownames(datExpr)), ]
genes=datExpr_reordered[,moduleColors==which.module ]

#plotting purple module and genotype 
PurpleME <- ggplot(MEs.sorted, aes(x = factor(genotype), y = MEpurple)) +  
  geom_boxplot(fill = "#cbabd1", outlier.shape = NA, size = 0.1) +  # Reduce line thickness
  geom_jitter(width = 0.2, size = .25, color = "#3d2443", alpha = 0.7) +  # Add jittered dots for each sample
  labs(x = "Genotype", y = "MEpurple") +  
  theme_minimal() +  # Apply minimal theme
  theme(
    panel.grid = element_blank(),  # Remove background grids
    axis.text.x = element_text(color = "black", size = 6),  # Reduce x-axis text size
    axis.text.y = element_text(color = "black", size = 6),  # Reduce y-axis text size
    axis.title.x = element_text(size = 7),  # Adjust x-axis title size
    axis.title.y = element_text(size = 7)   # Adjust y-axis title size
  )
PurpleME

survival_genotype_order <- c(36, 1, 50, 3, 44, 7, 31, 13, 62, 41)
MEs$genotype <- factor(MEs$genotype, levels = survival_genotype_order)
MEs.sorted <- MEs[order(MEs$genotype), ]

which.module= "tan" #pick module of interest

ME=MEs.sorted[, paste("ME",which.module, sep="")]
datExpr_reordered <- datExpr[match(rownames(MEs.sorted), rownames(datExpr)), ]
genes=datExpr_reordered[,moduleColors==which.module ]

##stats for purple module and genotype 
anova <- aov(MEtan ~ genotype, data = MEs.sorted)
summary(anova)
#Df Sum Sq Mean Sq F value Pr(>F)    
#genotype      9 0.7718 0.08576   96.97 <2e-16 ***
#  Residuals   258 0.2282 0.00088                   
tukey <- TukeyHSD(anova)
tukey_tan <- as.data.frame(tukey$genotype)

tukey_pvals <- tukey$genotype[, "p adj"]

# Generate letters
letters <- multcompLetters(tukey_pvals)$Letters

# Convert to dataframe
letters.df <- tibble(genotype = names(letters), letter = letters)

letters.df <- letters.df %>%
  mutate(genotype = factor(genotype, levels = survival_genotype_order)) %>%
  arrange(genotype)

TanME <- ggplot(MEs.sorted, aes(x = factor(genotype), y = MEtan)) +  
  geom_boxplot(fill = "#f1ddbf", outlier.shape = NA, size = 0.1) +  # Reduce line thickness
  geom_jitter(width = 0.2, size = .25, color = "#AB7E4C", alpha = 0.7) +  # Add jittered dots for each sample
  labs(x = "Genotype", y = "MEtan") +  
  theme_minimal() +  # Apply minimal theme
  theme(
    panel.grid = element_blank(),  # Remove background grids
    axis.text.x = element_text(color = "black", size = 6),  # Reduce x-axis text size
    axis.text.y = element_text(color = "black", size = 6),  # Reduce y-axis text size
    axis.title.x = element_text(size = 7),  # Adjust x-axis title size
    axis.title.y = element_text(size = 7)   # Adjust y-axis title size
  )
TanME

custom_order <- c(50, 3, 1, 31, 44, 36, 7, 41, 62, 13)
MEs$genotype <- factor(MEs$genotype, levels = custom_order)
MEs.sorted <- MEs[order(MEs$genotype), ]

which.module= "pink" #pick module of interest

ME=MEs.sorted[, paste("ME",which.module, sep="")]
datExpr_reordered <- datExpr[match(rownames(MEs.sorted), rownames(datExpr)), ]
genes=datExpr_reordered[,moduleColors==which.module ]

##stats for purple module and genotype 
anova <- aov(MEpink ~ genotype, data = MEs.sorted)
summary(anova)
#Df Sum Sq Mean Sq F value Pr(>F)    
#genotype      9 0.4768 0.05297   26.12 <2e-16 ***
#Residuals   258 0.5232 0.00203                     
tukey <- TukeyHSD(anova)
tukey_pink <- as.data.frame(tukey$genotype)
tukey_pvals <- tukey$genotype[, "p adj"]

# Generate letters
letters <- multcompLetters(tukey_pvals)$Letters

# Convert to dataframe
letters.df <- tibble(genotype = names(letters), letter = letters)

letters.df <- letters.df %>%
  mutate(genotype = factor(genotype, levels = custom_order)) %>%
  arrange(genotype)

pinkME <- ggplot(MEs.sorted, aes(x = factor(genotype), y = MEpink)) +  
  geom_boxplot(fill = "pink", outlier.shape = NA, size = 0.1) +  # Reduce line thickness
  geom_jitter(width = 0.2, size = .25, color = "#E75480", alpha = 0.7) +  # Add jittered dots for each sample
  labs(x = "Genotype", y = "MEpink") +  
  theme_minimal() +  # Apply minimal theme
  theme(
    panel.grid = element_blank(),  # Remove background grids
    axis.text.x = element_text(color = "black", size = 6),  # Reduce x-axis text size
    axis.text.y = element_text(color = "black", size = 6),  # Reduce y-axis text size
    axis.title.x = element_text(size = 7),  # Adjust x-axis title size
    axis.title.y = element_text(size = 7)   # Adjust y-axis title size
  )

pinkME

custom_order <- c(50, 3, 1, 31, 44, 36, 7, 41, 62, 13)
MEs$genotype <- factor(MEs$genotype, levels = custom_order)
MEs.sorted <- MEs[order(MEs$genotype), ]
which.module= "greenyellow" #pick module of interest

ME=MEs.sorted[, paste("ME",which.module, sep="")]
datExpr_reordered <- datExpr[match(rownames(MEs.sorted), rownames(datExpr)), ]
genes=datExpr_reordered[,moduleColors==which.module ]

anova <- aov(MEgreenyellow ~ genotype, data = MEs.sorted)
summary(anova)
#Df Sum Sq Mean Sq F value Pr(>F)    
#genotype      9  0.367 0.04078   16.62 <2e-16 ***
#Residuals   258  0.633 0.00245                     
tukey <- TukeyHSD(anova)
tukey_greenyellow <- as.data.frame(tukey$genotype)

tukey_pvals <- tukey$genotype[, "p adj"]

# Generate letters
letters <- multcompLetters(tukey_pvals)$Letters

# Convert to dataframe
letters.df <- tibble(genotype = names(letters), letter = letters)

letters.df <- letters.df %>%
  mutate(genotype = factor(genotype, levels = custom_order)) %>%
  arrange(genotype)

greenyellowME <- ggplot(MEs.sorted, aes(x = factor(genotype), y = MEgreenyellow)) +  
  geom_boxplot(fill = "greenyellow", outlier.shape = NA, size = 0.1) +  # Reduce line thickness
  geom_jitter(width = 0.2, size = .25, color = "forestgreen", alpha = 0.7) +  # Add jittered dots for each sample
  labs(x = "Genotype", y = "MEgreenyellow") +  
  theme_minimal() +  # Apply minimal theme
  theme(
    panel.grid = element_blank(),  # Remove background grids
    axis.text.x = element_text(color = "black", size = 6),  # Reduce x-axis text size
    axis.text.y = element_text(color = "black", size = 6),  # Reduce y-axis text size
    axis.title.x = element_text(size = 7),  # Adjust x-axis title size
    axis.title.y = element_text(size = 7)   # Adjust y-axis title size
  )
greenyellowME
library(patchwork)

# Compute y-axis limits
y_min <- min(MEs.sorted$MEpurple, MEs.sorted$MEtan, MEs.sorted$MEpink, MEs.sorted$MEgreenyellow, na.rm = TRUE)
y_max <- 0.25  # Set maximum limit to 0.3

# Apply consistent y-axis limits
PurpleME <- PurpleME + ylim(y_min, y_max)
TanME <- TanME + ylim(y_min, y_max)
pinkME <- pinkME + ylim(y_min, y_max)
greenyellowME <- greenyellowME + ylim(y_min, y_max)

# Combine plots horizontally
combined_plot <- PurpleME + TanME + pinkME + greenyellowME + 
  plot_layout(nrow = 1)

# Display the plot
combined_plot


########finding hub genes########
hub <- chooseTopHubInEachModule(datExpr, moduleColors)
#purple isogroupDN102558_c0_g1
#tan isogroupDN53976_c1_g1 Guanylate-binding protein 7 OS=Mus musculus OX=10090 GN=Gbp7 PE=1 SV=2 E(blastx)=8e-46
#pink isogroupDN4562_c0_g1
#greenyellow isogroupDN68057_c0_g1 

#############Network Analysis to functional annotation and gene ontology###############
##read in table with the annotated genes for the isogroups we mapped our reads to
# this is from the refs folder on the HPC - saved to personal computer using scp
annot=read.table("Acer_iso2go.tab",header=FALSE,sep="\t") #not all isogroups are annotated 
probes = names(datExpr)
probes2annot = match(probes,annot$V1) 
sum(is.na(probes2annot)) # 3278 NAs
datGS.Traits=data.frame(cor(datExpr,datTraits,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(datExpr,moduleColors)$eigengenes
datKME=signedKME(datExpr, datME, outputColumnName="MM.")
datOutput=data.frame(ProbeID=names(datExpr),annot[probes2annot,],moduleColors,datKME,datGS.Traits)
#datOutput=datOutput[order((datOutput$MM.black),decreasing=T),]
#write.table(datOutput,"AnnotatedNetworkAnalysisResultsSymb_rlog_signed_sft8_merge90_Jun22.csv",row.names=F,sep=",")

#Creating files for GO analysis of modules
##############Execute this entire block of code through creating VSD files 
#Host: purple, green, pink, black

## can also to rank-based GO with 0 for absent gene and then kME for the genes present in the module.
#below: creating output files for that.

### need to change color name here to module of interest.
purple=datOutput%>% 
  select("ProbeID","moduleColors","MM.purple")%>%
  mutate(kME = case_when(moduleColors == "purple" ~ MM.purple))%>%
  mutate(kME= ifelse(is.na(kME), 0, kME))%>%
  select("ProbeID", "kME")

#output here would be suitable for fisher's exact test in 
#makes values either 0 or 1 
purple = datOutput %>%
  select("ProbeID", "moduleColors", "MM.purple") %>%
  mutate(kME = case_when(moduleColors == "purple" ~ MM.purple)) %>%
  mutate(kME = ifelse(is.na(kME), 0, kME)) %>%
  mutate(kME = ifelse(kME > 0, 1, kME)) %>%
  select("ProbeID", "kME")

write.csv(purple, file = "unmerged_hostGO_MM_purple_numericASV.csv",row.names = F, quote=F)  

### need to change color name here to module of interest.
tan=datOutput%>% 
  select("ProbeID","moduleColors","MM.tan")%>%
  mutate(kME = case_when(moduleColors == "tan" ~ MM.tan))%>%
  mutate(kME= ifelse(is.na(kME), 0, kME))%>%
  select("ProbeID", "kME")
write.csv(tan, file = "unmerged_hostGO_MM_tan_numericASV.csv",row.names = F, quote=F)  

tan = datOutput %>%
  select("ProbeID", "moduleColors", "MM.tan") %>%
  mutate(kME = case_when(moduleColors == "tan" ~ MM.tan)) %>%
  mutate(kME = ifelse(is.na(kME), 0, kME)) %>%
  mutate(kME = ifelse(kME > 0, 1, kME)) %>%
  select("ProbeID", "kME") 
write.csv(tan, file = "fishers_unmerged_hostGO_MM_tan_numericASV.csv",row.names = F, quote=F)  


### need to change color name here to module of interest.
pink=datOutput%>% 
  select("ProbeID","moduleColors","MM.pink")%>%
  mutate(kME = case_when(moduleColors == "pink" ~ MM.pink))%>%
  mutate(kME= ifelse(is.na(kME), 0, kME))%>%
  select("ProbeID", "kME")
write.csv(pink, file = "unmerged_hostGO_MM_pink_numericASV.csv",row.names = F, quote=F)  

pink = datOutput %>%
  select("ProbeID", "moduleColors", "MM.pink") %>%
  mutate(kME = case_when(moduleColors == "pink" ~ MM.pink)) %>%
  mutate(kME = ifelse(is.na(kME), 0, kME)) %>%
  mutate(kME = ifelse(kME > 0, 1, kME)) %>%
  select("ProbeID", "kME") 
write.csv(pink, file = "fishers_unmerged_hostGO_MM_pink_numericASV.csv",row.names = F, quote=F)  



### need to change color name here to module of interest.
greenyellow=datOutput%>% 
  select("ProbeID","moduleColors","MM.greenyellow")%>%
  mutate(kME = case_when(moduleColors == "greenyellow" ~ MM.greenyellow))%>%
  mutate(kME= ifelse(is.na(kME), 0, kME))%>%
  select("ProbeID", "kME")
write.csv(greenyellow, file = "unmerged_hostGO_MM_greenyellow_numericASV.csv",row.names = F, quote=F)  

greenyellow = datOutput %>%
  select("ProbeID", "moduleColors", "MM.greenyellow") %>%
  mutate(kME = case_when(moduleColors == "greenyellow" ~ MM.greenyellow)) %>%
  mutate(kME = ifelse(is.na(kME), 0, kME)) %>%
  mutate(kME = ifelse(kME > 0, 1, kME)) %>%
  select("ProbeID", "kME") 
write.csv(greenyellow, file = "fishers_unmerged_hostGO_MM_greenyellow_numericASV.csv",row.names = F, quote=F) 
## above files go into the GO MWU analysis
##https://github.com/z0on/GO_MWU?tab=readme-ov-file 








#didnt use this 
#looking at most significant GO terms for each module - using GO MWU output
purple_terms <- read.csv("MWU_BP_unmerged_hostGO_MM_purple_numeric.csv",header=T, sep= " ")

#creating table with top 10 most significant GO terms for each module
#install.packages("flextable")
library(flextable)

purple_terms_paper<-purple_terms%>%
  arrange(pval)%>%
  select(pval, nseqs, term, name)

purple_terms_top10<-purple_terms_paper[c(1:10),]
purple_terms_top10$module="purple"



theme_table <- function(x) {
  border_remove(x) %>% 
    valign(valign = "center", part = "all") %>% 
    align(align = "center", part = "all") %>% 
    fontsize(part = "all", size = 10) %>% 
    bold(part = "header", bold = TRUE) %>%
    bold(part = "body", j = 1, bold = TRUE) %>% 
    color(color = "black", part = "header") %>% 
    bg(part = "header", bg = "transparent")
}


ft <- as_grouped_data(purple_terms_top10, groups = c("module"), expand_single = TRUE) %>% 
  as_flextable(hide_grouplabel = FALSE, col_keys = c("pval", "nseqs", "name")) %>% 
  set_header_labels(pval = "P-value", nseqs = "No. Sequences", name = "GO Term") %>%
  set_formatter(pval = function(x){
    formatC(x, format = "e", digits = 1)
  }) %>%
  theme_table() %>% 
  align(i = ~!is.na(name), align = "left", part = "body") %>% 
  bg(i = ~ module %in% "purple", bg = "#C3B1E1") %>%
  hline(i = rep(c(TRUE, TRUE, TRUE, TRUE), length = nrow_part(.)))
ft

width(ft, j = 3, width = 10, unit = "in")

save_as_image(width(ft, j = 3, width = 4, unit = "in"), "Final Figures/GO_table_red.pdf",  webshot = "webshot")


#Making VSD files by module for GO plot functions
Modcolors=c("purple")
for (col in Modcolors) {
  
  vs=t(datExpr)
  cands=names(datExpr[moduleColors==paste(col)])   
  
  c.vsd=vs[rownames(vs) %in% cands,]
  print(nrow(c.vsd)) #should correspond to module size
  write.csv(c.vsd,file=paste("host_vsd_MM_asv",col,".csv", sep=""),quote=F) #change to host or sym
}

write.csv(datTraits, file = "datTraits.csv", row.names = TRUE)

Modcolors=c("tan")
for (col in Modcolors) {
  
  vs=t(datExpr)
  cands=names(datExpr[moduleColors==paste(col)])   
  
  c.vsd=vs[rownames(vs) %in% cands,]
  print(nrow(c.vsd)) #should correspond to module size
  write.csv(c.vsd,file=paste("host_vsd_MM_asv",col,".csv", sep=""),quote=F) #change to host or sym
}

write.csv(datTraits, file = "datTraits.csv", row.names = TRUE)

Modcolors=c("pink")
for (col in Modcolors) {
  
  vs=t(datExpr)
  cands=names(datExpr[moduleColors==paste(col)])   
  
  c.vsd=vs[rownames(vs) %in% cands,]
  print(nrow(c.vsd)) #should correspond to module size
  write.csv(c.vsd,file=paste("host_vsd_MM_asv",col,".csv", sep=""),quote=F) #change to host or sym
}

write.csv(datTraits, file = "datTraits.csv", row.names = TRUE)

Modcolors=c("black")
for (col in Modcolors) {
  
  vs=t(datExpr)
  cands=names(datExpr[moduleColors==paste(col)])   
  
  c.vsd=vs[rownames(vs) %in% cands,]
  print(nrow(c.vsd)) #should correspond to module size
  write.csv(c.vsd,file=paste("host_vsd_MM_asv",col,".csv", sep=""),quote=F) #change to host or sym
}

write.csv(datTraits, file = "datTraits.csv", row.names = TRUE)

Modcolors=c("greenyellow")
for (col in Modcolors) {
  
  vs=t(datExpr)
  cands=names(datExpr[moduleColors==paste(col)])   
  
  c.vsd=vs[rownames(vs) %in% cands,]
  print(nrow(c.vsd)) #should correspond to module size
  write.csv(c.vsd,file=paste("host_vsd_MM_asv",col,".csv", sep=""),quote=F) #change to host or sym
}

write.csv(datTraits, file = "datTraits.csv", row.names = TRUE)
