#######Survival Analysis#####
###code adapted from Million et al 2022 https://github.com/wyattmillion/Acer_Morphological_Plasticity and Dilworth in prep 
setwd('/Users/natalievillafranca/Desktop/biomarkers')
rm(list=ls())
wyattsurv <- read.csv("AcerMorphologyData_Imputed2.csv")

library(tidyverse)
library(readxl)
library(janitor)
library(ggsci)
library(gtools)
library(RColorBrewer)
library(survminer)
library(survival)
library(coxme)
library(lubridate)
library(boot)

####JENNA CODE 
# survivorship ####
################### Kaplan-Meier survival curves
surv<-read_xlsx("ORCC_bleaching_mortality_Jul23_QCed.xlsx", sheet = "Mortality_binary")%>%
  clean_names()%>%
  dplyr::select(site, genotype, tag_number, replicate, jun_23_surv, jul_23_surv, tod_jun, tod_jul)%>%
  dplyr::rename(jun_23=jun_23_surv)%>%
  dplyr::rename(jul_23=jul_23_surv)%>%
  mutate_at("site", as.factor)%>%
  mutate_at("genotype", as.factor)%>%
  mutate_at("tag_number", as.factor)%>%
  mutate_at("replicate", as.factor)%>%
  dplyr::mutate(across(ends_with("23"), as.numeric))%>%
  dplyr::mutate(across(starts_with("tod"), as.numeric))%>%
  na.omit()
str(surv)

#order of survivorship in Wyatt's paper
# c("36","1","50","3","44","7","31","13","62","41"))

# survival under ambient conditions - ie until June 
surv$genotype<-factor(surv$genotype,levels=c("36","41","3", "50", "13","44","1","31","7","62"))

#AMBIENT SURVIVAL
ambient.died=Surv(surv$tod_jun, surv$jun_23)

###### Plotting survival curve
ambient.scurveG <- survfit(ambient.died~genotype, data=surv) #function cannot handle mixed effects or interactions

colorG <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(10)
jennasurvplot <- ggsurvplot(
  ambient.scurveG, 
  data = surv,
  palette = colorG,
  fun = "pct", 
  conf.int = FALSE,
  legend = "right",
  legend.title = "Genotype",
  ggtheme = theme_bw() +
    theme(
      legend.text = element_text(size = 8, family = "Helvetica Neue", face = "bold"),   
      legend.title = element_text(size = 9, family = "Helvetica Neue", face = "bold"), 
      axis.title.x = element_text(size = 10, family = "Helvetica Neue", face = "bold"),
      axis.title.y = element_text(size = 10, family = "Helvetica Neue", face = "bold")
    ),
  ylim = c(40, 100), 
  break.x.by = 3, 
  xlab = "Months"
)
genotype_labels <- gsub("Genotype=", "", levels(surv$genotype))

jennasurvplot$plot <- jennasurvplot$plot + 
  scale_color_manual(values = colorG, labels = genotype_labels)

ambient.res<-data.frame(summary(ambient.scurveG)$table)
genet<-str_split(rownames(ambient.res),"=",simplify=TRUE)[,2]
df<-data.frame(ambient.res,genet,stringsAsFactors=TRUE)

df$n.surv<-df$n.start-df$events
df$PercSurv<-df$n.surv/df$n.start
rownames(df)<-seq_len(nrow(df))
df.sort<-df[order(df$genet),]

## Cox proportional hazards
ambient.genotype.cox<-coxph(ambient.died~genotype, data=surv)



ambient.cox.rank <- data.frame(ambient.genotype.cox$coefficients)%>%
  rownames_to_column()%>%
  dplyr::rename(genotype=rowname)%>%
  rename(ambient_cox = ambient.genotype.cox.coefficients)






##wyatts data 
#only keeping looe and dave's 
wyattsurv_twosites <- subset(wyattsurv, Site %in% c("Looe Key", "Dave's Ledge"))

wyattsurv_twosites <- wyattsurv_twosites %>%
  mutate(T9_Dead = ifelse(T9_Status == "A", 0, ifelse(T9_Status == "D", 1, NA)))

wyattsurv_twosites <- wyattsurv_twosites %>%
  mutate(T9_TOD = ifelse(T12_TOD %in% c(12, 13), 10, T12_TOD))


wyattsurv_twosites$Genotype<-factor(wyattsurv_twosites$Genotype,levels=c("36","1","3","7", "31","44","50","41", "13","62"))

million.deaths=Surv(wyattsurv_twosites$T9_TOD, wyattsurv_twosites$T9_Dead)

million.scurveG <- survfit(million.deaths~Genotype, data=wyattsurv_twosites) #function cannot handle mixed effects or interactions

colorG <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(10)
wyattsurvplot <- ggsurvplot(million.scurveG, data = wyattsurv_twosites,
                            palette=colorG, fun = "pct", conf.int = FALSE, legend="right", legend.title="Genotype",
                            ggtheme = theme_bw() +
                              theme(
                                legend.text = element_text(size = 8, family = "Helvetica Neue", face = "bold"),  # Legend text bold
                                legend.title = element_text(size = 9, family = "Helvetica Neue", face = "bold"), # Legend title bold
                                axis.title.x = element_text(size = 10, family = "Helvetica Neue", face = "bold"), # X-axis label bold
                                axis.title.y = element_text(size = 10, family = "Helvetica Neue", face = "bold")  # Y-axis label bold
                              ),
                            break.x.by = 3, xlab = "Months")

genotype_labels <- gsub("Genotype=", "", levels(wyattsurv_twosites$Genotype))

wyattsurvplot$plot <- wyattsurvplot$plot + 
  scale_color_manual(values = colorG, labels = genotype_labels)

million.res<-data.frame(summary(million.scurveG)$table)
genet<-str_split(rownames(million.res),"=",simplify=TRUE)[,2]
df<-data.frame(million.res,genet,stringsAsFactors=TRUE)

df$n.surv<-df$n.start-df$events
df$PercSurv<-df$n.surv/df$n.start
rownames(df)<-seq_len(nrow(df))
df.sort<-df[order(df$genet),]

wyatt.genotype.cox<-coxph(million.deaths~Genotype, data=wyattsurv_twosites)

wyatt.cox.rank <- data.frame(wyatt.genotype.cox$coefficients)%>%
  rownames_to_column()%>%
  dplyr::rename(genotype=rowname)%>%
  rename(million_cox = wyatt.genotype.cox.coefficients)

##COMPARING WYATT AND JENNAS USING A SPEARMAN 
cor.test(ambient.cox.rank$ambient_cox, wyatt.cox.rank$million_cox, method=c("spearman"))
#rho = 0.85, p = 0.006074 

library(patchwork)

##combo these guys up 

wyatt_plot <- wyattsurvplot$plot
jenna_plot <- jenna_plot + theme(axis.title.y = element_blank())

# Combine the plots and add labels "A" and "B" on the top left
combined_plot <- (wyatt_plot + jenna_plot) +
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(family = "Helvetica", face = "bold", size = 14))

# Display the plot
combined_plot


