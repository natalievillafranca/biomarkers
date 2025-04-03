
# Install packages if needed - only do this once for each
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ANCOMBC")
#install.packages("pairwiseAdonis")
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")


########### Set Phylosq ######################
setwd('/Users/natalievillafranca/Desktop/biomarkers')

#reading in t0 nursery data 
seqtab.nochim1=readRDS("t0_seqtab2.rds")

taxa1=readRDS("t0_taxonomy.rds")

samdata1=read.csv("t0_sample2.csv")
rownames(samdata1)=samdata1$Sample_ID

#create phyloseq object
#BiocManager::install("phyloseq")
library(phyloseq)
library(ape)

ps1=phyloseq(otu_table(seqtab.nochim1, taxa_are_rows=FALSE), sample_data(samdata1), tax_table(taxa1))
ps1 #21925 taxa and 128 samples

#clean the phyloseq object as done previously (check out separate scripts for T0 and T1, respectively)
#prune this sample because I cannot account for its genotype as
#it was entered incorrectly into the system

ps1=subset_samples(ps1, sample_names(ps1) !="Gr7")

############################ Detour: Check mitochondria ###################################
#subset mitochondria 
mito=subset_taxa(ps1, Family =="Mitochondria") 
mito
mito=subset_samples(mito, sample_data(mito)$Sample_ID != "Gr7") #remove sample with unclear origin
mito=subset_samples(mito, sample_data(mito)$Genotype != "99") #remove tech replicate
mito=subset_taxa(mito, Phylum !="Ochrophyta") #protists
mito=subset_taxa(mito, Phylum !="Protalveolata") #protists
mito=subset_taxa(mito, Phylum !="Apicomplexa") #protists
mito=subset_samples(mito, sample_data(mito)$Genotype != "t12_99") #tech replicate
mito=subset_samples(mito, sample_data(mito)$Genotype != "t12_EDR  ") # FSW
mito=subset_samples(mito, sample_data(mito)$Genotype != "t12_ES  " ) #FSW
mito=subset_samples(mito, sample_data(mito)$Genotype != "t12_MS  " ) #FSW
mito

#Merge by type#
ps4<-merge_samples(mito, "Genotype")
ps4

library("ggplot2")
#plot to visualize distribution
#pmito=plot_bar(ps4, fill="Family")+  geom_bar(aes( fill=Family), stat="identity", position="stack")+ theme_bw()+
#  theme(text = element_text(size = 16))+
#  theme(axis.text.x = element_text(face="bold", 
#                                   size=10),
#        axis.text.y = element_text(face="bold",  
#                                   size=14))+ theme_bw() + theme(text = element_text(size = 16))+theme(legend.title = element_text(size = 14),
#                                                                                                       legend.text = element_text(size= 14)) + facet_wrap(~Timepoints, scales = "free_x")
#pmito

#if interested in creating a phylo tree or obtaining ASV 
#use Biostrings 

sequences <- Biostrings::DNAStringSet(taxa_names(ps4))
names(sequences) <- taxa_names(ps4)
ps.all <- merge_phyloseq(ps4, sequences)
ps.all
#use the ref_seq slots
rs <- refseq(ps.all)
rs # inspect
# Get strings with the full taxonomy for each OTU/ASV
tax <- tax_table(ps.all)
head(tax) # inspect
tax_strings <- apply(tax, 1, paste, collapse=";")
head(tax_strings) # inspect and BEAUTIFUL

# Create new names for the sequences in the refseq object; these will become
# the sequence headers in the fasta file. Here, I set these to be the OTU/ASV
# name and the tax string separated by one space (determined by the sep
# parameter)

taxa_names(ps4) <- paste0("ASV", seq(ntaxa(ps4)))
new_names <- paste(taxa_names(ps4), tax_strings, sep = " ")
head(new_names) # inspect
# Update the refeq names and save as a fasta file
names(rs) <- new_names
rs # inspect
#create fasta file with Biostrings!!!!
Biostrings::writeXStringSet(rs, "mitochondria.fasta")

######################### Back to our analysis  #######################################

#36613 taxa but what happens when
#filter out chloroplasts, mitochondria, protists and taxa counts 0 and under?

ps1=subset_taxa(ps1, (Class!="Chloroplast") | is.na(Class))
ps1 #21925 ASVs
ps1=subset_taxa(ps1, Family !="Mitochondria") #12495 ASVs taxa left
ps1
ps1=prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1 #now we have 12487 taxa and 127 samples samples,#if interested in chloroplasts or mitochondrial sequences, we can subset and explore later 

ps1=subset_taxa(ps1, Phylum !="Ochrophyta") #protists
ps1=subset_taxa(ps1, Phylum !="Protalveolata") #protists
ps1=subset_taxa(ps1, Phylum !="Apicomplexa") #protists

ps1=prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
#12487 taxa left


#gonna take out <5000 reads

#take out low reads
#let's toss out some samples < 5000 
#should have a total of  124 samples 
summary(sample_sums(ps1)> 5000) #106 samples left, 21 rejected
#write.csv(sample_sums(j2ps), "five000andabove_samplesums.csv")
ps1=subset_samples(ps1, sample_sums(ps1)>5000)
ps1 #12487 taxa and 106 samples
sample_data(ps1)

#write.csv(sample_data(joint_ps), "joint_samdata.csv")
#gonna take out water samples and anything not related to timepoint T0 and T12, which I labeled for clarity T665 and T666

#alternate to include water, 21554 ASVs and 191 samples
#j2ps=subset_samples(ps1, sample_data(ps1)$Sample_ID != "Gr7")
#j2ps=subset_samples(j2ps, sample_data(j2ps)$Genotype != "t12_99")

#to not include filtered seawater
j2ps=subset_samples(ps1, sample_data(ps1)$Timepoint !="T665")
j2ps  #103 samples
j2ps=subset_samples(j2ps, sample_data(j2ps)$Timepoint !="T666")
j2ps  #12487 taxa and 103 samples


#take out ASVs with 0 counts out in our samples and rarefy
ps=prune_taxa(taxa_sums(j2ps) > 0, j2ps)
length(get_taxa_unique(ps1, taxonomic.rank = "Genus"))
#[550]
#add a tree
library("ape")
random_tree=rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
#merge the tree
ps=merge_phyloseq(ps, random_tree)
ps  #11763 ASVs taxa without the water 

#FROM CARLY CODE
##########################   Rarefy   ###########################################
#find out how many reads in each sample before rarefying
set.seed(123)
ps_rare <- phyloseq::rarefy_even_depth(ps, rngseed = 123, replace = FALSE)
#2634OTUs were removed because they are no longer present in any sample after random subsampling
ps_rare #9129 ASVs 

plot(sort(sample_sums(ps_rare),TRUE),type="h", ylab="reads") 

##################### Nursery analysis  ####################
### What is initial community composition of nursery genets?
library(dplyr)
library(tidyr)
write.table(ps_rare %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
              arrange(OTU) %>% rename(ASV = OTU) %>% 
              select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Sample, Abundance) %>%
              spread(Sample, Abundance), 
            file = "ps_rare.relative_abundance.all.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

#quick check of dominant taxa
T <- ps_rare %>% 
  tax_glom(., "Genus") %>% #change taxonomic selection
  transform_sample_counts(function(x)100* x / sum(x)) %>% psmelt() %>% 
  arrange(OTU) %>% rename(OTUsID = OTU) %>% 
  select(OTUsID, Genus, Sample, Abundance) %>% #change taxonomic selection
  spread(Sample, Abundance)

T$Mean <- rowMeans(T[, c(3:ncol(T))])

FAM <- T[, c("Genus", "Mean" ) ] #change taxonomic selection

#order from major to minor 
FAM <- FAM[order(dplyr::desc(FAM$Mean)),]
rownames(FAM) <- NULL
head(FAM,n=10) #this shows most abundant classes from top to bottom
#  Synechococcus CC9902 37.876255
#                              MD3-55 26.998847
#                    NS5 marine group  7.858502
#             Candidatus Actinomarina  4.320972
#                              HIMB11  3.239160
summary(FAM)
sum(FAM$Mean)

tax_glom(ps_rare, taxrank = 'Class')

# rarefied
#filter out taxa that doesn't appear more than 5x in more than 1 sample. 
ps_rare=filter_taxa(ps_rare, function(x) sum(x > 5) > (0.005*length(x)), TRUE)
ps_rare ### 3700 ASVs 

#write out file of relative abundances by taxonomic length 
#for total counts use (function(x) {x/sum(x)})
#for relative abundance use (function(x)100* x / sum(x))
T <- ps_rare %>% 
  tax_glom(., "Genus") %>% #change taxonomic selection
  transform_sample_counts(function(x)100* x / sum(x)) %>% psmelt() %>% 
  arrange(OTU) %>% rename(OTUsID = OTU) %>% 
  select(OTUsID, Genus, Sample, Abundance) %>% #change taxonomic selection
  spread(Sample, Abundance)
T 
write.csv(T, "TaxonRelAbunByGenusSample_ps_rare.csv") #for a record or other analysis, like corrplots

#add a tree
random_tree=rtree(ntaxa(ps_rare), rooted=TRUE, tip.label=taxa_names(ps_rare))
#merge the tree
ps_rare=merge_phyloseq(ps_rare, random_tree)
ps_rare

## transform to relative abundance
trans_filt_ps_rare<-transform_sample_counts(ps_rare, function(x) x /sum(x))
sample_sums(trans_filt_ps_rare) #all sum to 1 as expected

##pretty plot color options
library(randomcoloR)
n <- 103
palette <- distinctColorPalette(n)
pie(rep(1, n), col=palette)

colorbfriend= c( "darkgreen", "deepskyblue4", "darkgoldenrod1", "aquamarine3", "coral2", "cadetblue", "plum4", "springgreen4", "aliceblue", "darkcyan", "darksalmon", "darkslateblue", "gold1", "honeydew2", "darkolivegreen4", "indianred3", "bisque1", "palevioletred4", "palegreen4", "maroon", "deepskyblue3", "darkred", "paleturquoise4", "burlywood2", "darkseagreen", "hotpink4", "lightcyan3", "darkslategray", "antiquewhite3", "goldenrod2", "royalblue4", "thistle", "forestgreen", "lightblue", "lightcoral", "lightseagreen", "peachpuff", "midnightblue", "rosybrown2", "mediumseagreen", "tan", "red4")
colorbfriend
####

ps2_rare<-merge_samples(trans_filt_ps_rare, "Genotype")
ps2_rare

p1=plot_bar(ps2_rare, fill="Phylum")+ scale_fill_manual(values=colorbfriend)
p1
phylumplot= p1+ geom_bar(aes(fill=Phylum), stat="identity", position="stack") 
phylumplot + theme_bw()+ 
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.text.y = element_text(size=14)) + 
  theme(text = element_text(size = 16))+theme(legend.title = element_text(size = 14))+facet_grid(~Genotype, scale="free_x", drop=TRUE)

#END CARLY CODE 



#################### Alpha Diversity  ############################################
alph <- estimate_richness(ps_rare, measures = c("Observed", "Shannon", "Simpson", "Chao1", "Fisher"))
alph
write.csv(alph, "alpha_diversity_ps_rare.csv") #for a record or other analysis, like corrplots

#plot richness by measure
plot_richness(ps_rare, x="Genotype")

p1= plot_richness(ps_rare, x="Genotype", measures= c("Shannon"))+ geom_violin(fill="grey90")+ geom_jitter(width=0.15, alpha=0.8,size=2)+ theme(axis.text.x = element_text( 
  size=20),
  axis.text.y = element_text(  
    size=20))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=24))+
  ylab("Alpha Diversity Measure")+ xlab ("Genet") +theme_bw()

p1

#is the data normal?
hist(alph$Shannon, main="Shannon diversity", breaks=10) #yes it is. 

#Now we want to know if genotypes are different from each other (a-diversity, Shannon)

#use Genotypes only subsetted data (see lines 172-180 above for more info)
#So just do ANOVA, followed by Tukey HSD. If data wasn't normal, you'd do a Kruskal Wallis test
#like below the aov, using vegan

my_factors = data.frame(sample_data(ps_rare))
summary(ps_rare)
aov_alpha_genets=aov(alph$Shannon ~ my_factors$Genotype)
summary(aov_alpha_genets) 
#yes 
#who is different? not very many. only 13-1, 62-1, 31-13, 62-31. 
Tukey <- TukeyHSD(aov_alpha_genets)
plot(Tukey , las=1 , col="brown")


#What about for different sites? 
aov_alpha_genets=aov(alph$Shannon ~ my_factors$Sites)
summary(aov_alpha_genets) #nope

######### what about beta diversity among nursery genets?

genotype_colors<-c("paleturquoise4", "mediumaquamarine", "seagreen", "lightgoldenrod", "goldenrod", "lightsalmon", "salmon", "rosybrown1","thistle", "rosybrown") 
# plot
set.seed(123)
#Ordination techniques, such as principal coordinates analysis (PCoA), 
#reduce the dimensionality of microbiome data sets so that a summary of the 
#beta diversity relationships can be visualized in two- or three-dimensional scatterplots.
ordinate<- ordinate(ps_rare, method = "PCoA", distance ="wunifrac") #weighted-UniFrac distance. 

#Weighted-UniFrac takes into account the relative abundance of species/taxa shared between samples, 
#whereas unweighted-UniFrac only considers presence/absence.

evals <- ordinate$values$Eigenvalues
ps_rare_betadiv=betadiv=plot_ordination(ps_rare, ordinate, color = "Genotype") 
# Reorder the Genotype factor levels in the desired order for legend
  ps_rare@sam_data$Genotype <- factor(ps_rare@sam_data$Genotype, 
                                      levels = c("36", "1", "50", "3", "44", "7", "31", "13", "62", "41"))
# Create the ordination plot with the reordered legend
ps_rare_betadiv <- plot_ordination(ps_rare, ordinate, color = "Genotype") + 
  scale_color_manual(values = genotype_colors) +
  labs(col = "Genet") + 
  geom_point(size = 4) + 
  stat_ellipse(aes(color = Genotype)) +
  theme_bw() + 
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14))


# Print the plot
print(ps_rare_betadiv)

#write.csv(ordinate$vectors, "beta_diversity_vectors_PCoA_wunifrac_filt_ps_rare_T0.csv") #for a record or other analysis, like corrplots


library(vegan)
set.seed(123)
dfbray=phyloseq::distance(ps_rare, method="bray") #Bray-Curtis distances are weighted by abundance
sampledf=data.frame(sample_data(ps_rare))
#betadisper in vegan implements PERMADISP for the analysis of multivariate homogeneity of group dispersions (variances)
#PERMDISP involves calculating the distance from each data point to its group centroid and then testing whether those distances differ among the groups.
beta=betadisper(dfbray, sampledf$Genotype)
beta
permutest(beta) #p=0.007, even dispersion (?)

#use PERMANOVA to test for differences in beta diversity, first confirm data satisfy test assumptions
#PERMANOVA compares the variation between groups to the variation within groups.
#adonis is an analysis of variance using distance matrices — 
#for partitioning distance matrices among sources of variation and fitting linear models to distance matrices
#a type of permanova 
adonis2(formula = dfbray ~ Genotype, data = sampledf)#p=0.001, yes sig differences between genets 
#but, have 10 levels, therefore many pairwise contrasts to perform
#solution: use pairwiseAdonis
vegan::adonis2(dfbray ~ phyloseq::sample_data(ps_rare)$Genotype) #verifying same result as above
library(pairwiseAdonis)
PERMANOVA <- pairwise.adonis(dfbray, phyloseq::sample_data(ps_rare)$Genotype) #now pairwise test
#weak significance when you compare genotypes to each other. 
write.csv(PERMANOVA, file ="permanova.csv")

################## Relative abundance ############### there is code from carlys above, 
#i think the difference is that carly's is from transformed data


#Merge samples based on a sample variable or factor -- merge by genotype 
ps2_rare<-merge_samples2(ps_rare, "Genotype")
ps2_rare

colorbfriend= c( "#E5F5F9", "#1D91C0", "#67001F", "#F7FCFD", "#CB181D", "#78C679", "#F46D43", "#A6CEE3", "#FD8D3C", "#A6D854", "#D4B9DA", "#6A51A3", "#7F0000", "#D9D9D9", "#FFF7BC", "#000000", "#F0F0F0", "#C7EAE5", "#003C30", "#F16913", "#FFF7FB", "#8C6BB1", "#C7E9B4", "#762A83", "#FC9272", "#AE017E", "#F7F7F7", "#DF65B0", "#EF3B2C", "#74C476", "#FF62BC")
colorbfriend


#plot relative abundance, phylum
p2=plot_bar(ps2_rare, fill="Phylum")
p2
phylumplot= p2+geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")+ scale_fill_manual(values = colorbfriend)
phylumplot + theme_bw()+ theme(text = element_text(size = 16))+ theme_bw()+
  scale_x_discrete(name ="Genotype", 
                   limits=c("50", "3", "1", "31", "44", "36", "7", "41", "62", "13")
) +
  theme(axis.text.x = element_text(face="bold", 
                                   size=10),
        axis.text.y = element_text(face="bold",  
                                   size=14))+ theme_bw() + theme(text = element_text(size = 16))+theme(legend.title = element_text(size = 14),
                                                                                                       legend.text = element_text(size= 14))


#plot top 100 ASVs
#preferably using the same color of the phyla for correspondent genera
#extract the colors 
#Get the scale for the class palette to match up with the colors in phyla 

library(scales)
n2 <- 18                                               # Higher amount of hex colors
hex_codes2 <- hue_pal()(n2)                             # Identify hex codes
show_col(hex_codes2) 
topcolors=c("#F564E3", "#00BF74")

top100 <- names(sort(taxa_sums(ps2_rare), decreasing=TRUE))[1:100]
ps.top100 <- transform_sample_counts(ps2_rare, function(OTU) OTU/sum(OTU))
ps.top100 <- prune_taxa(top100, ps.top100)
#plot genera
genus_top100=plot_bar(ps.top100, fill="Genus")+ scale_color_manual(values = topcolors)
genus_top100 = genus_top100 + 
  geom_bar(aes(x = factor(Sample, levels = c("50", "3", "1", "31", "44", "36", "7", "41", "62", "13")), 
               fill = Genus), 
           stat = "identity", position = "stack") + 
  theme_bw() +
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 14)) + 
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

genus_top100

