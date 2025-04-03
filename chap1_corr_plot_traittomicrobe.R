setwd('/Users/natalievillafranca/Desktop/biomarkers')
options(stringsAsFactors=FALSE)
library("conflicted")
library(tidyverse)
library(readxl)
library(dplyr) 
library(WGCNA)
library(corrplot)

lnames = load(file="BiomarkersTraits_Nov24.RData") 

##corrplotting 
names(datTraits)
order <- c("T3_Status", "T12_Status", "T3_TLE", "T12_TLE", 
           "T3_Vinter", "T12_Vinter", "Synechococcus CC9902", "MD3-55",
           "NS5 marine group", "Candidatus Actinomarina", "HIMB11", "Shannon")

corrplot_datTraits <- datTraits[, order]
corrplot_datTraits <- corrplot_datTraits[!is.na(corrplot_datTraits$"MD3-55"), ]
corrplot_datTraits <- corrplot_datTraits[!is.na(corrplot_datTraits$T3_Status), ]

library(WGCNA)
### correlation stats ####
cor1 <- WGCNA::cor(corrplot_datTraits, method = c("pearson"), use = "pairwise.complete.obs")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(corrplot_datTraits)

library(pheatmap)
par(mar = c(1, 1, 4, 2))
####correlation heatmap####
library(corrplot)
col <- colorRampPalette(c("#4477AA","#77AADD","#FFFFFF","#EE9988", "#BB4444"))
heatmap <- corrplot(cor1,
                    tl.col = "black", 
                    tl.srt = 45, 
                    tl.cex = 0.5, 
                    p.mat = p.mat, 
                    insig = "label_sig", 
                    sig.level = c(.001, .01, .05), 
                    pch.cex = 0.8, 
                    pch.col = "white", 
                    type = "lower")


datTraits <- datTraits[!is.na(datTraits$`MD3-55`) & !is.na(datTraits$T3_Status), ]

hist(datTraits$`MD3-55`)
####subsequent stats####
result <- glm(T3_Status ~ `MD3-55`, family=binomial(link="logit"), data = datTraits)
summary(result) 




















#
# Create the bar plot for MD-55

ggplot(data_na, aes(x = genotype, y = `MD3-55`, fill = as.factor(T3_Status))) +
  geom_boxplot(position = "dodge") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black")+
  labs(fill = "T3 Status") + 
  theme_minimal()

data_no_na <- data %>%
  filter(!is.na(T12_Status))
data_no_na$genotype <- factor(data_no_na$genotype, levels = survival_order)
data_no_na <- data_no_na[order(data_no_na$genotype), ]
ggplot(data_no_na, aes(x = genotype, y = `MD3-55`, fill = as.factor(T12_Status))) +
  geom_boxplot(position = "dodge") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black")+
  labs(fill = "T12 Status") + 
  theme_minimal()

ggplot(data_na, aes(x = genotype, y = `Synechococcus CC9902`, fill = as.factor(T3_Status))) +
  geom_boxplot(position = "dodge") +
  labs(fill = "T3 Status") + 
  theme_minimal()

ggplot(data_na, aes(x = genotype, y = `Synechococcus CC9902`, fill = as.factor(T12_Status))) +
  geom_boxplot(position = "dodge") +
  labs(fill = "T12 Status") + 
  theme_minimal()















###T3 Volume
rosnerTest(data$T3_V, k = 5)
data_noout <- dplyr::filter(data, !T3_Vinter %in% c(31.79200))

r <- WGCNA::cor(data_noout$Candidatus.Puniceispirillum, data_noout$T3_V, method = "pearson", use = "pairwise.complete.obs")
p <- cor.test(data_noout$Candidatus.Puniceispirillum, data_noout$T3_V)$p.value
ggplot(data_noout, aes(T3_V, Candidatus.Puniceispirillum)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black", lty = 1) + 
  annotate("text", x=0.03, y=0.017, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=0.03, y=0.016, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic() 

#### T3 Vinter

rosnerTest(data$T3_Vinter, k = 5)
data_noout1 <- dplyr::filter(data, !T3_Vinter %in% c(105.44300))

r <- WGCNA::cor(data_noout1$Candidatus.Puniceispirillum, data_noout1$T3_Vinter, method = "pearson", use = "pairwise.complete.obs")
p <- cor.test(data_noout1$Candidatus.Puniceispirillum, data_noout1$T3_Vinter)$p.value
ggplot(data_noout1, aes(T3_Vinter, Candidatus.Puniceispirillum)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black", lty = 1) + 
  annotate("text", x=0.03, y=0.017, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=0.03, y=0.016, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic() 

###
#checking for outliers 
rosnerTest(data$T12_V, k = 5)
data_noout2 <- dplyr::filter(data, !T12_V %in% c(149.129,120.061, 114.954))

#correlation test 
r <- WGCNA::cor(data_noout2$Candidatus.Puniceispirillum, data_noout2$T12_V, method = "pearson", use = "pairwise.complete.obs")
p <- cor.test(data_noout2$Candidatus.Puniceispirillum, data_noout2$T12_V)$p.value
#plotting
ggplot(data_noout2, aes(T12_V, Candidatus.Puniceispirillum)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black", lty = 1) + 
  annotate("text", x=0.03, y=0.017, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=0.03, y=0.016, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic()

###
#checking for outliers 
rosnerTest(data$T12_Vinter, k = 5)
data_noout3 <- dplyr::filter(data, !T12_Vinter %in% c(1077.425,948.851, 796.728,434.731))
#correlation test 
r <- WGCNA::cor(data_noout3$Candidatus.Puniceispirillum, data_noout3$T12_Vinter, method = "pearson", use = "pairwise.complete.obs")
p <- cor.test(data_noout3$Candidatus.Puniceispirillum, data_noout3$T12_Vinter)$p.value
#plotting
ggplot(data_noout3, aes(T12_Vinter, Candidatus.Puniceispirillum)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black", lty = 1) + 
  annotate("text", x=0.03, y=0.017, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=0.03, y=0.016, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic()

###
#Synechococcus CC9902 with status 
rosnerTest(data$T3_Status, k = 5)#none
rosnerTest(data$T12_Status, k = 5)#none

#t3 status 
# t test
r <- WGCNA::cor(data$Synechococcus.CC9902, data$T3_Status, method = "pearson", use = "pairwise.complete.obs")
p <- cor.test(data$Synechococcus.CC9902, data$T3_Status)$p.value
ggplot(data %>% filter(!is.na(T3_Status)), aes(x = factor(T3_Status), y = Synechococcus.CC9902)) + 
  geom_boxplot(na.rm = TRUE) + 
  annotate("text", x=0.2, y=0.05, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=0.2, y=0.016, label=paste0("p = ", round(p, 3)), hjust=0) +  theme_classic()

#t12 status
r <- WGCNA::cor(data$Synechococcus.CC9902, data$T12_Status, method = "pearson", use = "pairwise.complete.obs")
p <- cor.test(data$Synechococcus.CC9902, data$T12_Status)$p.value
ggplot(data %>% filter(!is.na(T12_Status)), aes(x = factor(T12_Status), y = Synechococcus.CC9902)) + 
  geom_boxplot(na.rm = TRUE) + 
  annotate("text", x=0.2, y=0.05, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=0.2, y=0.016, label=paste0("p = ", round(p, 3)), hjust=0) + theme_classic()

##
#md355 with status
#t3
r <- WGCNA::cor(data$MD3.55, data$T3_Status, method = "pearson", use = "pairwise.complete.obs")
p <- cor.test(data$MD3.55, data$T3_Status)$p.value
ggplot(data %>% filter(!is.na(T3_Status)), aes(x = factor(T3_Status), y = MD3.55)) + 
  geom_boxplot(na.rm = TRUE) + 
  annotate("text", x=0.2, y=0.05, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=0.2, y=0.016, label=paste0("p = ", round(p, 3)), hjust=0) +  theme_classic()



#t12
r <- WGCNA::cor(data$MD3.55, data$T12_Status, method = "pearson", use = "pairwise.complete.obs")
p <- cor.test(data$MD3.55, data$T12_Status)$p.value

ggplot(data %>% filter(!is.na(T12_Status)), aes(x = factor(T12_Status), y = MD3.55)) + 
  geom_boxplot(na.rm = TRUE) + 
  annotate("text", x=0.2, y=0.05, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=0.2, y=0.016, label=paste0("p = ", round(p, 3)), hjust=0) +  theme_classic()

###what is best way to rn stats on this... ##shapiro test will tell you if something is normal, if not significant it is normal 
# use t test to compare 

#ALSO RUN STATS ON GENES IF POSSIBLE 

beta <- read.csv("beta_diversity_vectors_PCoA_wunifrac_filt_ps_rare_T0.csv")







