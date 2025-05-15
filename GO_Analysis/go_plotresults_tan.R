# this is the final script in GO analysis using MWU test, as in
# Voolstra et al PLoS ONE 2011, 6(5): e20392.

setwd('/Users/natalievillafranca/Desktop/biomarkers')

####setting data up####
#Note: must run this script while results of GO_MWU are loaded into R
#must run these to initiate, then can ignore...
gg=read.table("Acer_iso2gene.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)  #the iso2gene file for your species
dfull=read.csv("host_vsd_MM_asvtan.csv") #read in proper VSD file for genes in module - generated in wgcna script 
rownames(dfull)<-dfull$X #make gene names your rownames
efulldata=c(2:269) #the columns where your expression data are 2:59 host; 2:56 sym
library(pheatmap)

###############plotting what significant GO genes are
#"BP_unmerged_hostGO_MM_purple_numeric.csv" is without ASV 
#MWU_BP_unmerged_hostGO_MM_purple_numeric.csv is without ASV
golist=read.table("BP_unmerged_hostGO_MM_tan_numericASV.csv",sep="	",header=T) #read in proper GO list from gomwu output - BP/MF/CC and module name
in.mwu=read.table("MWU_BP_unmerged_hostGO_MM_tan_numericASV.csv",sep=" ",header=T)
d=read.csv("host_vsd_MM_asvtan.csv")#read in proper VSD file for genes in module - generated in wgcna script 
edata=c(2:269) #only columns with your expression data 2:56 host; 2:58 sym


gene=subset(golist,name=="negative regulation of receptor signaling pathway via JAK-STAT") #write GO term of interest here from your sig list#
t=rownames(gene)
t

is <- golist[rownames(golist) %in% t, c("seq", "term")]
length(is) #then pulls out the isogroups with this annotation

#### OR pulling from alternate lists (skip this if using GO term above)
#is<-c("isogroupDN102558_c0_g1", "isogroupDN72717_c4_g1", "isogroupDN74671_c2_g1", "isogroupDN57704_c0_g1", "isogroupDN74674_c1_g3") 
#in order: hub gene IRF, IRF T12vinter and prochlor, e. selectin, allograft inflammatory factor, tyrosine protein kinase receptor 


#####################first loop through genes matching GO term in module or in "good expression" subset
#if using specific isogroups (alternate lists)
#sel=c();gnms=c()
#for ( i in is){
#  if (i %in% d$X){
#    sel=rbind(sel,d[d$X==i,])
#    gnms=append(gnms,substr(paste(as.character(gg$V2[gg$V1==i]),collapse="."),1,50))
#  }
#}
#row.names(sel)=paste(gnms,sel$X,sep=".")
#nrow(sel)
#rownames(sel)

sel <- data.frame()  # Initialize as an empty data frame
gnms <- c()

for (i in is$seq) {  # Assuming 'seq' column in 'is' contains valid identifiers
  if (i %in% d$X) {
    sel <- rbind(sel, d[d$X == i, ])
    matched_values <- gg$V2[gg$V1 == i]
    gnms <- append(gnms, substr(paste(as.character(matched_values), collapse = "."), 1, 50))
  }
}

# Ensure rownames are assigned properly
if (nrow(sel) > 0) {
  rownames(sel) <- paste(gnms, sel$X, sep = ".")
}
print(nrow(sel))
print(rownames(sel))

exp=sel[,efulldata]
if (length(exp[,1])==1) { exp=rbind(exp,exp) }
nrow(exp)
#should match 'good genes' tabulation in figure
rownames(exp)

#------------------ heatmap of GOgenes stacked by symbiont density (or other trait...)
library(pheatmap)

#reading in traits 
lnames = load(file = "BiomarkersTraits_Nov24.RData") 
library(dplyr)

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



#36, 1, 50, 3, 44, 7, 31, 13, 62, 41
datTraits=datTraits%>%relocate("36", .after = "T12_Vinter")
datTraits=datTraits%>%relocate("1", .after ="36")
datTraits=datTraits%>%relocate("50", .after ="1")
datTraits=datTraits%>%relocate("3", .after ="50")
datTraits=datTraits%>%relocate("44", .after ="3")
datTraits=datTraits%>%relocate("7", .after ="44")
datTraits=datTraits%>%relocate("31", .after ="7")
datTraits=datTraits%>%relocate("13", .after ="31")
datTraits=datTraits%>%relocate("62", .after ="13")
datTraits=datTraits%>%relocate("41", .after ="62")


#genotype only 
datTraits1=datTraits[,c(1)]

#totresval_genotype_order <- c(50, 3, 1, 31, 44, 36, 7, 41, 62, 13)
survival_genotype_order <- c(36, 1, 50, 3, 44, 7, 31, 13, 62, 41)
datTraits$genotype <- factor(datTraits$genotype, levels = survival_genotype_order)
datTraits <- datTraits[order(datTraits$genotype),]

#traits <- traits[order(traits$Chao1, decreasing = TRUE),] #only use when using a microbe as a trait 
#traits_noNA <- traits[!is.na(traits$Chao1), ] #only use when using a microbe as a trait 
####
sorted=exp[,rownames(datTraits)] #sort rows by values of a trait *** if you don't want to sort, just say sorted=exp

#sorted<-exp #do this if you do not want to sort by predefined trait

############# Center scale
means=apply(sorted,1,mean) # means of rows
#hist(means)
expc=sorted-means #rescale expression data



library(RColorBrewer)
col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=1)(30) #adjust bias to center 0

#setting up annotation_col option to get legends saying what genotype and status each sample is 
#only if you are trying to see gene correlation with genotype as a trait 

Genotype <- datTraits$genotype
annotation_col <- data.frame(Genotype = Genotype)
rownames(annotation_col) <- colnames(expc)  # Ensure rownames match the column names of expc
annotation_colors=list(Genotype=c("50"="#053061", "3"="#2166ac", "1"="#4393c3", "31"="#92c5de", 
                                    "44"="#d1e5f0", "36"="#fddbc7", "7"="#f4a582", "41"="#d6604d",
                                     "62"="#b2182b", "13"="#67001f"))
genotype_breaks <- cumsum(table(Genotype))
pheatmap(
  expc, 
  color = col,
  cluster_cols = FALSE, 
  clustering_distance_rows = "binary",
  show_colnames = FALSE,
  border_color = "black",  # Adds black borders to all cells
  gaps_col = genotype_breaks # Adds lines between genotype groups
)

####heatmap of genes ####
#this will cluster in whatever order the dataframe is in with the legends provided above for annotation_col
pheatmap(expc, color = col, cluster_cols = FALSE,
         annotation_col = annotation_col, clustering_distance_rows = "binary", 
         show_colnames = FALSE)

gene=as.data.frame(t(sorted))
texp=data.frame(t(sorted))
sexp=stack(texp)
unique(sexp$ind) #unique returns a list of all the gene names in your GOterm for the module of choice

final <- bind_cols(gene,
                   genotype = datTraits$genotype, 
                   Prochlorococcus=datTraits$`Prochlorococcus MIT9313`,
                   Synechococcus=datTraits$`Synechococcus CC9902`,
                   MD355=datTraits$`MD3-55`, 
                   NS5 = datTraits$`NS5 marine group`,
                   HIMB11 = datTraits$HIMB11, 
                   CandidatusActinomarina = datTraits$`Candidatus Actinomarina`, 
                   Shannon = datTraits$Shannon)

#save(final, file="gene2microbe.RData")

####GE correlation with genotype ####
#running stats to see how correlated GE is with genotype and how different exp is between genotypes 
#with anova and post hoc tukey, not goint to check for normalization because we normalized counts in WGCNA with DESEQ2
#hub gene IRF2
AOV <- aov(final$`Guanylate-binding protein 7 OS=Mus musculus OX=100.isogroupDN53976_c1_g1` ~ final$genotype)
AIC(AOV) 
Tuk<- TukeyHSD(AOV)
tukey_results <- Tuk$`final$genotype`
pvals <- tukey_results[, "p adj"]
names(pvals) <- rownames(tukey_results)

cld_vec <- multcompLetters(pvals)$Letters
letter_df <- data.frame(
  genotype = names(cld_vec),
  cld = as.character(cld_vec),
  stringsAsFactors = FALSE
)

letter_df$genotype <- factor(letter_df$genotype, levels = totresval_genotype_order)
letter_df <- letter_df[order(letter_df$genotype), ]
print(letter_df)
#genotype cld
#       50   a
#        3   b
#        1   a
#       31   a
#       44   c
#       36   a
#        7   a
#       41   d
#       62   c
#       13   a

#

AOV <- aov(final$`Interferon regulatory factor 1 OS=Mus musculus OX=.isogroupDN73115_c4_g1` ~ final$genotype)
AIC(AOV) 
Tuk<- TukeyHSD(AOV)
tukey_results <- Tuk$`final$genotype`
pvals <- tukey_results[, "p adj"]
names(pvals) <- rownames(tukey_results)

cld_vec <- multcompLetters(pvals)$Letters
letter_df <- data.frame(
  genotype = names(cld_vec),
  cld = as.character(cld_vec),
  stringsAsFactors = FALSE
)

letter_df$genotype <- factor(letter_df$genotype, levels = totresval_genotype_order)
letter_df <- letter_df[order(letter_df$genotype), ]
print(letter_df)

#genotype cld
#       50  ab
#        3   c
#        1  ab
#       31  ab
#       44   d
#       36 bde
#        7  ac
#       41   f
#       62  de
#       13 abe

#
AOV <- aov(final$`Protein mono-ADP-ribosyltransferase PARP14 OS=Homo.isogroupDN72552_c1_g1` ~ final$genotype)
AIC(AOV) 
Tuk<- TukeyHSD(AOV)
tukey_results <- Tuk$`final$genotype`
pvals <- tukey_results[, "p adj"]
names(pvals) <- rownames(tukey_results)

cld_vec <- multcompLetters(pvals)$Letters
letter_df <- data.frame(
  genotype = names(cld_vec),
  cld = as.character(cld_vec),
  stringsAsFactors = FALSE
)

letter_df$genotype <- factor(letter_df$genotype, levels = totresval_genotype_order)
letter_df <- letter_df[order(letter_df$genotype), ]
print(letter_df)
#genotype cld
#        50   a
#         3   a
#         1 abc
#        31 abc
#        44   b
#       36 abc
#         7  bc
#        41   d
#        62 bcd
#        13  ac

