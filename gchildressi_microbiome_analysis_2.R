#Microbial community dynamics during key life-history transitions in the chemosymbiotic mussel, Gigantidas childressi

setwd("~/Downloads/Microbiome Analysis")

#installing and loading packages#
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

install.packages('devtools', repos='http://cran.rstudio.com/')
devtools::install_github("jbisanz/qiime2R")

if (!requireNamespace("BiocManager", quietly = TRUE)) 
install.packages("BiocManager") 
BiocManager::install("phyloseq")

packageVersion("DESeq2")
BiocManager::install("DESeq2", force=TRUE)

library(ggplot2)
library(devtools)
library(qiime2R)
library(phyloseq)
library(DEseq2)
library(vegan)
library(BiocManager)
library(DESeq2)
library(knitr)
library(pairwiseAdonis)
library(dendextend)


#Loading QIIME2 data into R w/ Phyloseq
physeq<-qza_to_phyloseq(
  features="filtered-table.qza",
  tree="rooted-tree.qza",
  "classification.qza",
  metadata = "magical_metadata_final.tsv")

#flagging this
#Remove blank samples
#physeq = subset_samples(physeq, individual != "0")

#Variance Stabilizing Transformation w/ DESeq2
diagdds = phyloseq_to_deseq2(physeq, ~ group)
diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)
dim(diagvst)
physeq.vst<-physeq
otu_table(physeq.vst) <- otu_table(diagvst, taxa_are_rows = TRUE)

#Re-ordering samples
sample_data(physeq)$group<-factor(sample_data(physeq)$group, levels = c("Gill","Juvenile","Pediveliger","Veliger", "Water"))
sample_data(physeq.vst)$group<-factor(sample_data(physeq.vst)$group, levels = c("Gill","Juvenile","Pediveliger","Veliger", "Water"))

# ALPHA DIVERSITY 
# Histogram of Read Counts 
hist(sample_sums(physeq), main="Histogram: Read Counts", xlab="Total Reads", 
     border="black", col="lightgreen", las=1, breaks=14)

# Strip plot of total ASVs by group 
theme_set(theme_classic())
plot_richness(physeq, x="group", measures = "Observed") 

#Rarefaction curve 
sample_data(physeq)$group
colors<-c("#E69F00","#F0E442","#999999","#E69F00","#009E73","#F0E442","#56B4E9","#009E73","#E69F00","#56B4E9","#009E73","#F0E442","#F0E442","#999999","#999999","#E69F00","#E69F00","#999999","#F0E442","#009E73","#999999")
rarecurve(t(otu_table(physeq)), step=50, cex=1, label = FALSE, col = colors, ylab = "ASV's",xlab="Library Size",xlim=c(0, 15000)) 

#Kruskal wallis test and pairwise wilcox
alphaObserved = estimate_richness(physeq, measures="Observed")
alpha.stats <- cbind(alphaObserved, sample_data(physeq))
kruskal.test(alpha.stats$Observed, alpha.stats$group, data = alpha.stats)
pairwise.wilcox.test(alpha.stats$Observed, alpha.stats$group, p.adjust.method = "BH")


# BETA DIVERSITY #
#PCoA
uunifrac <- phyloseq::distance(physeq.vst, method = "uunifrac")
sampledf <- data.frame(sample_data(physeq.vst))

# Permanova test on Beta Diversity 
adonis2(uunifrac ~ group, data = sampledf)
pairwise.adonis2(uunifrac ~ group, data = sampledf)

theme_set(theme_classic())
ordu1 = ordinate(physeq.vst, "PCoA", "unifrac", weighted=FALSE)
plot_ordination(physeq.vst, ordu1, axes = c(1,2), color = "group") + geom_point(size = 3)

#ANCOM-BC2 Analysis
library("ANCOMBC")
#Remove adult gill and water samples
physeq_juv_vel_pedi<- subset_samples(physeq,
                         sample_names(physeq) != "WA6100" &
                           sample_names(physeq) != "WA6098" &
                           sample_names(physeq) != "WA6099" &
                           sample_names(physeq) != "WA6101" &
                           sample_names(physeq) != "WA6102" &
                           sample_names(physeq) != "56229" &
                           sample_names(physeq) != "56230")

#re-order sample groups depending on desired contrasts
sample_data(physeq_juv_vel_pedi)$group<-factor(sample_data(physeq_juv_vel_pedi)$group, levels = c("Veliger","Pediveliger","Juvenile"))
group<-sample_data(physeq_gill_juv_vel_pedi)$group

#Run ANCOMBC test
out = ancombc(data = NULL, assay_name = NULL,
               tax_level = "Genus", phyloseq = physeq_juv_vel_pedi,
               formula = "group",
               p_adj_method = "holm", prv_cut = 0.05,
               group = "group", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
               n_cl = 1, verbose = TRUE)

res = out$res
res_global = out$res_global

#obtain a table of adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 6)




#Hierarchical cluster analysis

# This is the actual hierarchical clustering call, specifying average-link clustering
MM.hclust     <- hclust(uunifrac, method="average")
dend <- as.dendrogram(MM.hclust)
labels(dend) <- c("Niskin","Niskin", "Juvenile","Gill","Gill","Gill","Gill","Gill","Juvenile","Juvenile","Juvenile","Juvenile","Pediveliger","Pediveliger","Veliger","Pediveliger","Pediveliger","Pediveliger","Veliger","Veliger","Veliger")
dend %>% set("leaves_col", 1) %>% # adjust the leaves
  hang.dendrogram %>% # hang the leaves
  plot()

# Relative abundance plots
#Barplot

#remove niskins
physeq1<- subset_samples(physeq,
                  sample_names(physeq) != "56229" &
                    sample_names(physeq) != "56230")

#remove taxa with unknown families
physeq2 <- subset_taxa(physeq1, !is.na(Family) & !Family %in% c("", "uncharacterized"))

#transform values into relative abundances
physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )
otu_table(physeq3)
                                  
#sort top 50 ASVs#
topN = 50
most_abundant_taxa = sort(taxa_sums(physeq3), TRUE)[1:topN]
print(most_abundant_taxa)
physeq4 = prune_taxa(names(most_abundant_taxa), physeq3)
otu_table(physeq4)

## Family Level ##
glom <- tax_glom(physeq4, taxrank = 'Family')

glom # should list # taxa as # Family 
data_glom<- psmelt(glom) # create dataframe from phyloseq object

data_glom$Family <- as.character(data_glom$Family) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Family[data_glom$Abundance < 0.01] <- "< 1% abund."

# Clear out rows
data_glom[data_glom==0] <- NA
data_glom <- data_glom[complete.cases(data_glom),]


#Count # phyla to set color palette
Count = length(unique(data_glom$Family))
unique(data_glom$Family)

data_glom$Family <- factor(data_glom$Family, levels = c("Methylomonadaceae",
                                                        "Moraxellaceae",
                                                        "Pseudoalteromonadaceae",
                                                        "Sphingomonadaceae",
                                                        "Pseudomonadaceae", 
                                                        "Helicobacteraceae",
                                                        "Methylophagaceae",
                                                        "Oxalobacteraceae", 
                                                        "Cycloclasticaceae",
                                                        "Thioglobaceae", 
                                                        "Endozoicomonadaceae",
                                                        "Chloroplast",
                                                        "Bacillaceae",
                                                        "Alteromonadaceae",
                                                        "Spirochaetaceae",
                                                        "< 1% abund."))    

col_vector<- c("#7FC97F", "#BEAED4", "#FDC086","gold1" , "skyblue2",
               "#F0027F","#386CB0" , "#FF7F00", "#1B9E77", "brown",
               "#7570B3", "#FFFF99", "#66A61E", "black", "#A6761D",
               "gray70", "#A6CEE3", "#1F78B4", "#E31A1C")

spatial_plot <- ggplot(data=data_glom, aes(x=individual, y=Abundance, fill=Family)) + facet_grid(~group, scales = "free")
spatial_plot + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = col_vector) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=4)) +
  theme(legend.title=element_text(size=9)) +
  theme(legend.text=element_text(size=8)) 

# Heatmap 

# Define the ranks you want to include in the heatmap
myranks = c("Family","Genus")
mylabels = apply(tax_table(physeq4)[, myranks], 1, paste, sep="", collapse="_")
# Add concatenated labels as a new rank after strain
tax_table(physeq4) <- cbind(tax_table(physeq4), Taxonomy=mylabels)

sample_names(physeq4)
order<-c("WA6102","WA6101",
         "WA6100","WA6099","WA6098","56228","56227","56226","56225","56224","56223","56222","56221", "56220","56219","56218","56217","56216", "56215","56230","56229")

plot_heatmap(physeq4, taxa.label="Taxonomy", 
                   taxa.order="Phylum", low="royalblue2", high="magenta", na.value = "black", sample.order=order, sample.label="group")
