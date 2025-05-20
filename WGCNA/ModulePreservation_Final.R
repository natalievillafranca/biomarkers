#module preservation
rm(list=ls()) 
setwd('/Users/natalievillafranca/Desktop/biomarkers')

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

lnames = load(file = "May_T12_Biomarkers_Expr_andTraits.RData")

lnames = load(file = "BiomarkersTraitsAndInteractionsAssigned_T12_March.RData")


###SETTING UP FOR MODULE PRESERVATION
#only keeping datTraitsT12 columns that match with datTraits columns (isogroups)
#have to combine and then separate because T12 has more isogroups (columns) and the number 
#of isogroups have to match those of datTraits
setdiff(colnames(datExpr), colnames(datExprT12))

datExprT12 <- datExprT12 %>% select(any_of(colnames(datExpr)))

setdiff(colnames(datExpr), colnames(datExprT12))


#fixing the genotype name error
colnames(datTraits)[colnames(datTraits) == "34"] <- "36"

#combining datTraits 
#datTraitsCombined <- bind_rows(datTraits, datTraitsT12)

which(colnames(datExpr) == "isogroupDN59835_c1_g1") #corresponds to column 1281 
moduleColors <- moduleColors[-1281]

#removing isogroupDN59835_c1_g1 from both T0 and T12 datExpr
datExpr <- datExpr[, !colnames(datExpr) %in% "isogroupDN59835_c1_g1"]
datExprT12 <- datExprT12[, !colnames(datExprT12) %in% "isogroupDN59835_c1_g1"]

###module preservation 

multiExpr = list(A1=list(data=datExpr),A2=list(data=datExprT12))
multiColor = list(A1 = moduleColors) 

gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
                      nPermutations=30,maxGoldModuleSize=100,maxModuleSize=500) 
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2 
stats[order(-stats[,2]),c(1:2)]

#maxModuleSize at 500 
#moduleSize Zsummary.pres
#pink                447     29.776564
#tan                 115     27.381062
#brown               500     27.132658
#green               500     25.711277
#greenyellow         125     23.472379
#black               500     22.598101
#blue                500     20.628364
#salmon              112     18.803863
#purple              173     15.414257
#cyan                500     15.276403
#midnightblue         47     15.238363
#magenta             500     12.685186
#turquoise           500      9.818397
#gold                100      6.463661


colors = names(table(moduleColors))

# First calculate MEs for 2021, since we haven't done that yet
PCs2_T12 = moduleEigengenes(datExprT12, colors=moduleColors)
ME2_T12 = PCs2_T12$eigengenes

geneModuleMembership2 = signedKME(datExprT12, ME2_T12)
colnames(geneModuleMembership2)=paste("PC",colors,".cor",sep="")
MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(t(datExprT12))[[2]])
colnames(MMPvalue2)=paste(colors,".pval",sep="")
kMEtable2 = cbind(Gene,Gene,moduleColors)

for (i in 1:length(colors)){
  kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i]) 
}
colnames(kMEtable2)=colnames(kMEtable1)

#write.csv(kMEtable2,"kMEtable2.csv",row.names=FALSE)


lnames = load(file = "datTraitsT12.RData") 


###now going to plot purple module eigengene expression by genotype for T12 
###Visualization of Gene Networks
############################################
#heatmaps of module expression with bar plot of eigengene
#use this to look at the different samples and make sure the modules make sense
#we want blocks by treatment, no single sample driving the differences.

#start with modules that looked strongest: #purple looks strongest for T3_Status, 
#green, pink, and black are best for T3_TLE
library(tibble)
#sort ME by traits 
ME2_T12$sort=rownames(ME2_T12)
datTraitsT12 <- datTraitsT12 %>%
  rownames_to_column("SRR")

# Joining using full_join
ME2_T12 <- full_join(ME2_T12, datTraitsT12, by = c("sort" = "SRR"))
print(ME2_T12)

ME2_T12 <- ME2_T12 %>%
  column_to_rownames("sort")

###############
#MEpurple <- MEs.sorted[, paste("ME", which.module, sep = "")]

################

#if you want to order by genotype ranking 
survival_genotype_order <- c(36, 1, 50, 3, 44, 7, 31, 13, 62, 41)

#visualizing module eigengene  expression with boxplots 
# Create a data frame for boxplot
library(ggplot2)
library(tidyverse)
library(egg)
library(multcompView)
#create boxplot using ggplot2

custom_order <- c(50, 3, 1, 31, 44, 36, 7, 41, 62, 13)
ME2_T12$genotype <- factor(ME2_T12$genotype, levels = custom_order)
MEs.sorted <- ME2_T12[order(ME2_T12$genotype), ]

anova <- aov(MEpurple ~ genotype, data = MEs.sorted)
summary(anova)
#Df Sum Sq Mean Sq F value Pr(>F)    
#genotype      9 0.8258 0.09175   92.68 <2e-16 ***
#Residuals   176 0.1742 0.00099                  
tukey <- TukeyHSD(anova)
tukey_df <- as.data.frame(tukey$genotype)
tukey_pvals <- tukey$genotype[, "p adj"]

# Generate letters
letters <- multcompLetters(tukey_pvals)$Letters

# Convert to dataframe
letters.df <- tibble(genotype = names(letters), letter = letters)

letters.df <- letters.df %>%
  mutate(genotype = factor(genotype, levels = custom_order)) %>%
  arrange(genotype)

#custom_order <- c(50, 3, 1, 31, 44, 36, 7, 41, 62, 13)
#MEs.sorted$genotype <- factor(MEs.sorted$genotype, levels = custom_order)
#MEs.sorted <- MEs[order(MEs.sorted$genotype), ]

which.module= "purple" #pick module of interest

ME=MEs.sorted[, paste("ME",which.module, sep="")]
datExpr_reordered <- datExprT12[match(rownames(MEs.sorted), rownames(datExprT12)), ]
genes=datExpr_reordered[,moduleColors==which.module ]

library(multcompView)
##stats for purple module and genotype 
anova <- aov(MEpurple ~ genotype, data = MEs.sorted)
summary(anova)
#Df Sum Sq Mean Sq F value Pr(>F)    
#genotype      9 0.8258 0.09175   92.68 <2e-16 ***
#Residuals   176 0.1742 0.00099                    
tukey <- TukeyHSD(anova)
tukey_purple <- as.data.frame(tukey$genotype)
letters.df <- data_frame(multcompLetters2(TukeyHSD(aov(MEpurple ~ genotype, data = MEs.sorted))$genotype[,4])$Letters)


#plotting purple module and genotype 
T12PurpleME <- ggplot(MEs.sorted, aes(x = factor(genotype), y = MEpurple)) +  
  geom_boxplot(fill = "#cbabd1", outlier.shape = NA, size = 0.1) +  # Reduce line thickness
  geom_jitter(width = 0.2, size = .25, color = "#3d2443", alpha = 0.7) +  # Add jittered dots for each sample
  labs(x = "Genotype", y = "T12 MEpurple") +  
  theme_minimal() +  # Apply minimal theme
  theme(
    panel.grid = element_blank(),  # Remove background grids
    axis.text.x = element_text(color = "black", size = 6),  # Reduce x-axis text size
    axis.text.y = element_text(color = "black", size = 6),  # Reduce y-axis text size
    axis.title.x = element_text(size = 7),  # Adjust x-axis title size
    axis.title.y = element_text(size = 7)   # Adjust y-axis title size
  )

survival_genotype_order <- c(36, 1, 50, 3, 44, 7, 31, 13, 62, 41)
ME2_T12$genotype <- factor(ME2_T12$genotype, levels = survival_genotype_order)
MEs.sorted <- ME2_T12[order(ME2_T12$genotype), ]

which.module= "tan" #pick module of interest

ME=MEs.sorted[, paste("ME",which.module, sep="")]
datExpr_reordered <- datExpr[match(rownames(MEs.sorted), rownames(datExprT12)), ]
genes=datExpr_reordered[,moduleColors==which.module ]

##stats for purple module and genotype 
anova <- aov(MEtan ~ genotype, data = MEs.sorted)
summary(anova)
#Df Sum Sq Mean Sq F value Pr(>F)    
#genotype      9 0.7718 0.08576   96.97 <2e-16 ***
#  Residuals   258 0.2282 0.00088                   
tukey <- TukeyHSD(anova)
tukey_df <- as.data.frame(tukey$genotype)
tukey_tan <- as.data.frame(tukey$genotype)

tukey_pvals <- tukey$genotype[, "p adj"]

# Generate letters
letters <- multcompLetters(tukey_pvals)$Letters

# Convert to dataframe
letters.df <- tibble(genotype = names(letters), letter = letters)

letters.df <- letters.df %>%
  mutate(genotype = factor(genotype, levels = survival_genotype_order)) %>%
  arrange(genotype)

T12TanME <- ggplot(MEs.sorted, aes(x = factor(genotype), y = MEtan)) +  
  geom_boxplot(fill = "#f1ddbf", outlier.shape = NA, size = 0.1) +  # Reduce line thickness
  geom_jitter(width = 0.2, size = .25, color = "#AB7E4C", alpha = 0.7) +  # Add jittered dots for each sample
  labs(x = "Genotype", y = "T12 MEtan") +  
  theme_minimal() +  # Apply minimal theme
  theme(
    panel.grid = element_blank(),  # Remove background grids
    axis.text.x = element_text(color = "black", size = 6),  # Reduce x-axis text size
    axis.text.y = element_text(color = "black", size = 6),  # Reduce y-axis text size
    axis.title.x = element_text(size = 7),  # Adjust x-axis title size
    axis.title.y = element_text(size = 7)   # Adjust y-axis title size
  )

custom_order <- c(50, 3, 1, 31, 44, 36, 7, 41, 62, 13)
ME2_T12$genotype <- factor(ME2_T12$genotype, levels = custom_order)
MEs.sorted <- ME2_T12[order(ME2_T12$genotype), ]

which.module= "pink" #pick module of interest

ME=MEs.sorted[, paste("ME",which.module, sep="")]
datExpr_reordered <- datExpr[match(rownames(MEs.sorted), rownames(datExprT12)), ]
genes=datExpr_reordered[,moduleColors==which.module ]

anova <- aov(MEpink ~ genotype, data = MEs.sorted)
summary(anova)


T12pinkME <- ggplot(MEs.sorted, aes(x = factor(genotype), y = MEpink)) +  
  geom_boxplot(fill = "pink", outlier.shape = NA, size = 0.1) +  # Reduce line thickness
  geom_jitter(width = 0.2, size = .25, color = "#E75480", alpha = 0.7) +  # Add jittered dots for each sample
  labs(x = "Genotype", y = "T12 MEpink") +  
  theme_minimal() +  # Apply minimal theme
  theme(
    panel.grid = element_blank(),  # Remove background grids
    axis.text.x = element_text(color = "black", size = 6),  # Reduce x-axis text size
    axis.text.y = element_text(color = "black", size = 6),  # Reduce y-axis text size
    axis.title.x = element_text(size = 7),  # Adjust x-axis title size
    axis.title.y = element_text(size = 7)   # Adjust y-axis title size
  )

custom_order <- c(50, 3, 1, 31, 44, 36, 7, 41, 62, 13)
ME2_T12$genotype <- factor(ME2_T12$genotype, levels = custom_order)
MEs.sorted <- ME2_T12[order(ME2_T12$genotype), ]

which.module= "greenyellow" #pick module of interest

ME=MEs.sorted[, paste("ME",which.module, sep="")]
datExpr_reordered <- datExpr[match(rownames(MEs.sorted), rownames(datExprT12)), ]
genes=datExpr_reordered[,moduleColors==which.module ]

anova <- aov(MEgreenyellow ~ genotype, data = MEs.sorted)
summary(anova)

T12greenyellowME <- ggplot(MEs.sorted, aes(x = factor(genotype), y = MEgreenyellow)) +  
  geom_boxplot(fill = "greenyellow", outlier.shape = NA, size = 0.1) +  # Reduce line thickness
  geom_jitter(width = 0.2, size = .25, color = "forestgreen", alpha = 0.7) +  # Add jittered dots for each sample
  labs(x = "Genotype", y = "T12 MEgreenyellow") +  
  theme_minimal() +  # Apply minimal theme
  theme(
    panel.grid = element_blank(),  # Remove background grids
    axis.text.x = element_text(color = "black", size = 6),  # Reduce x-axis text size
    axis.text.y = element_text(color = "black", size = 6),  # Reduce y-axis text size
    axis.title.x = element_text(size = 7),  # Adjust x-axis title size
    axis.title.y = element_text(size = 7)   # Adjust y-axis title size
  )

library(patchwork)

buffer <- 0.5  

# Compute y-axis limits
y_min <- min(MEs.sorted$MEpurple, MEs.sorted$MEtan, MEs.sorted$MEpink, MEs.sorted$MEgreenyellow, na.rm = TRUE)
y_max <- 0.3  # Set maximum limit to 0.3

# Apply consistent y-axis limits
T12PurpleME <- T12PurpleME + ylim(y_min, y_max)
T12TanME <- T12TanME + ylim(y_min, y_max)
T12pinkME <- T12pinkME + ylim(y_min, y_max)
T12greenyellowME <- T12greenyellowME + ylim(y_min, y_max)

# Combine plots horizontally
combined_plot <- T12PurpleME + T12TanME + T12pinkME + T12greenyellowME + 
  plot_layout(nrow = 1)

# Display the plot
combined_plot


#############Network Analysis to functional annotation and gene ontology###############
##read in table with the annotated genes for the isogroups we mapped our reads to
# this is from the refs folder on the HPC - saved to personal computer using scp
annot=read.table("Acer_iso2go.tab",header=FALSE,sep="\t") #not all isogroups are annotated 
probes = names(datExprT12)
probes2annot = match(probes,annot$V1) 
sum(is.na(probes2annot)) # 3277 NAs
datGS.Traits=data.frame(cor(datExprT12,datTraitsT12,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(datExprT12,moduleColors)$eigengenes
datKME=signedKME(datExprT12, datME, outputColumnName="MM.")
datOutput=data.frame(ProbeID=names(datExprT12),annot[probes2annot,],moduleColors,datKME,datGS.Traits)
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

write.csv(purple, file = "T12unmerged_hostGO_MM_purple_numericASV.csv",row.names = F, quote=F)  

Modcolors=c("purple")
for (col in Modcolors) {
  
  vs=t(datExprT12)
  cands=names(datExprT12[moduleColors==paste(col)])   
  
  c.vsd=vs[rownames(vs) %in% cands,]
  print(nrow(c.vsd)) #should correspond to module size
  write.csv(c.vsd,file=paste("T12host_vsd_MM_asv",col,".csv", sep=""),quote=F) #change to host or sym
}
