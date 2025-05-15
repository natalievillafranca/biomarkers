#######Survival Analysis#####
###code adapted from Million et al 2022 https://github.com/wyattmillion/Acer_Morphological_Plasticity
setwd('/Users/natalievillafranca/Desktop/biomarkers')
surv <- read.csv("AcerMorphologyData_Imputed2.csv")
#install.packages("corr")
library(survival) # for survival analysis 
library(RColorBrewer)
library(survminer)
str(surv)
library(dplyr)

surv$T12_TOD<-as.numeric(surv$T12_TOD)

####Making a survival object where TOD is time of death event, and Dead is whether or not a fragment was dead
score<-surv[complete.cases(surv$T12_TOD),] #removing ramets Missing from outplant 
score$Genotype <- factor(score$Genotype)
score$Genotype<-relevel(score$Genotype,ref="36")## picked genet with best survival as reference
died=Surv(score$T12_TOD, score$T12_Dead)
scurveG<- survfit(died~Genotype,data=score)
scurveG
#n events median 0.95LCL 0.95UCL
#Genotype=36 25      1     NA      NA      NA #96 
#Genotype=1  20      2     NA      NA      NA #90
#Genotype=50 26      3     NA      NA      NA #88.5
#Genotype=3  26      3     NA      NA      NA #88.5
#Genotype=44 25      3     NA      NA      NA #88
#Genotype=7  23      4     NA      NA      NA #82.6
#Genotype=31 26      5     NA      NA      NA #80.8 
#Genotype=13 25      8     NA      NA      NA #68
#Genotype=62 24      9     NA       9      NA #62.5
#Genotype=41 26     10     NA       9      NA #61.5

score$Genotype<-factor(score$Genotype,levels=c("36","1","50","3","44","7","31","13","62","41"))## picked genet with best survival as reference

## total restoration value rankings 
# Remove rows where T0_TLE is NA to remove samples where we have sequence data but not field data 
#this includes samples where there 0 outplant data as well as samples that have T0 data but then no T3 data (lost)
surv <- surv[!is.na(surv$T0_TLE), ]
surv <- surv[surv$T12_Status != "M", ]
#changing NAs (just 0) at T12 to values 0 
surv <- surv %>%
  mutate(T12_Vinter = ifelse(is.na(T12_Vinter), 0, T12_Vinter),
         T9_Vinter = ifelse(is.na(T9_Vinter), 0, T9_Vinter),
         T6_Vinter = ifelse(is.na(T6_Vinter), 0, T6_Vinter),
         T3_Vinter = ifelse(is.na(T3_Vinter), 0, T3_Vinter))


#creating restoration value metric --> sum of Vcx at T12 in surviving ramets by genotype 
#ordering by genotype from best to worst survivor 
Genotypes <- c(36, 1, 50, 3, 44, 7, 31, 13, 62, 41)
surv$Genotype <- factor(surv$Genotype, levels = Genotypes)
surv <- surv[order(surv$Genotype), ]

# Initialize an empty dataframe for results
TotResValue <- data.frame(Genotype = integer(), sum_T12Vinter = numeric())

# Loop through each genotype and calculate the sum of T12 Vinter
for (Genotype in Genotypes) {
  sum_T12_Vinter <- sum(surv$T12_Vinter[surv$Genotype == Genotype], na.rm = FALSE)
  TotResValue <- rbind(TotResValue, data.frame(Genotype = Genotype, sum_T12_Vinter = sum_T12_Vinter))
}
#best to worst: 50, 3, 1, 31, 44, 36, 7, 41, 62, 13

#spearman correlation -- comparing rank orders 
Ranks <- data.frame(
  genotypes = c(36, 1, 50, 3, 44, 7, 31, 13, 62, 41),
  surv_ranks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
  resval_ranks =c(6, 3, 1, 2, 5, 7, 4, 10, 9, 8)
)
cor.test(Ranks$surv_ranks, Ranks$resval_ranks, method = "spearman") # p-value = 0.03509, rho =0.6848485
#the survival rank and resvalue ranks numbers are significantly correlated 

JennaRanks <- data.frame(
  genotypes = c(36, 1, 50, 3, 44, 7, 31, 13, 62, 41),
  old_ranks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
  new_ranks =c(1, 10, 4, 3, 8, 5, 2, 7, 6, 9)
)
cor.test(JennaRanks$old_ranks, JennaRanks$new_ranks, method = "spearman") # p-value = 0.03509, rho =0.6848485

#pearson correlation -- comparing the values themselves 
death_v_vinter <- data.frame(
  genotypes = c(36, 1, 50, 3, 44, 7, 31, 13, 62, 41),
    deaths = c(1, 2, 3, 3, 3, 4, 5, 8, 9, 10), 
    vinter = c(2697.192, 3078.816, 6197.946, 3669.834, 2682.318, 2895.230, 2944.329, 796.824, 2253.453, 1623.347)
)
#checking for normality for pearson correlation
shapiro.test(death_v_vinter$deaths)
shapiro.test(death_v_vinter$vinter)
#both values are normally distributed 
cor.test(death_v_vinter$deaths, death_v_vinter$vinter, method = "pearson") #p = p-value = 0.08436, cor = -0.5715047 
#plotting 
death_v_vinter$genotypes <- factor(death_v_vinter$genotypes, levels = unique(death_v_vinter$genotypes))
ggplot(data = death_v_vinter, aes(x = deaths, y = vinter)) +
  geom_point(aes(fill = genotypes), shape = 21, color = "black", size = 3) +  # Use fill in geom_point
  geom_smooth(method = "lm", color = "blue") +  # Add linear regression without fill aesthetic
  scale_fill_brewer(palette = "RdBu") +
  theme_minimal() +
  labs(fill = "Genotypes")

