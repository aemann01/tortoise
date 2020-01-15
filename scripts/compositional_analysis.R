# Compositional data analysis
# Tortoise dataset

###########
#LIBRARIES
###########
library(tidyverse)
library(ggplot2)
library(ape)
library(factoextra)
library(RColorBrewer)
library(phyloseq)
library(UpSetR)
library(reshape2)
library(vegan)
library(ggtree)
library(DESeq2)
library(igraph)
library(microbiome)
library(philr)
library(SpiecEasi)

##########
#ENV SETUP
##########
PATH="/home/lymelab/lab_members/mann/tortoise"
setwd(PATH)
cols <- brewer.pal(6, "Set2")

###########
#LOAD DATA
###########
#load sample data
rawmetadata <- read_delim(file = file.path(PATH, "map_final_filt.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries
rownames(rawmetadata) <- rawmetadata$sample_id 
seqtab.filtered <- read.table("sequence_table.16s.filtered.txt", header=T, row.names=1)
system("sed -i 's/;/\t/g' taxonomy_L7.txt")
taxa <- read.table("taxonomy_L7.txt", header=F, sep="\t", row.names=1)
notinmeta <- setdiff(row.names(seqtab.filtered), rawmetadata$sample_id) #record samples absent in either metadata or OTU table
notinmeta
notinraw <- setdiff(rawmetadata$sample_id, row.names(seqtab.filtered))
notinraw
tree <- read.tree("rep_set.root.tre") #load representative tree

#filter out samples not in metadata
# seqtab.filtered <- subset(seqtab.filtered, !(row.names(seqtab.filtered) %in% notinmeta))
# dim(seqtab.filtered)
# dim(rawmetadata)

#create phyloseq object from "seqtab.filtered", "rawmetadata", "taxa", and "tree"
ps.dat <- phyloseq(otu_table(seqtab.filtered, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(as.matrix(taxa)), tree)

################
#CORE MICROBIOME
################
#across all necropsy and solar tortoises
ps.dat.trans <- microbiome::transform(ps.dat, "clr")
head(prevalence(ps.dat.trans, detection=1/100, sort=T))
#      ASV1     ASV12     ASV36      ASV5     ASV17      ASV9
# 1.0000000 0.8589744 0.7948718 0.7435897 0.7051282 0.6666667
head(prevalence(ps.dat.trans, detection=1/100, sort=T, count=T))
 # ASV1 ASV12 ASV36  ASV5 ASV17  ASV9
 #   78    67    62    58    55    52
core.taxa.standard <- core_members(ps.dat.trans, detection = 0, prevalence = 50/100)
ps.dat.core <- core(ps.dat.trans, detection=0, prevalence=.5)
prevalences <- seq(0.05, 1, 0.05)
detections <- 10^seq(log10(1e-3), log10(.2), length=10)
pdf("figs/core_all.pdf")
plot_core(ps.dat.trans, plot.type="heatmap", prevalences=prevalences, detections=detections, colours=rev(brewer.pal(5, "Spectral")), min.prevalence=.2, horizontal=T)
dev.off()

#necropsy only
nec.only <- subset_samples(ps.dat.trans, tortoise_group=="necropsy")
pdf("figs/core_nec.pdf")
plot_core(nec.only, plot.type="heatmap", prevalences=prevalences, detections=detections, colours=rev(brewer.pal(5, "Spectral")), min.prevalence=.2, horizontal=T)
dev.off()
head(prevalence(nec.only, detection=1/100, sort=T))
#      ASV1     ASV12     ASV36      ASV3      ASV5     ASV17
# 1.0000000 0.9148936 0.8936170 0.8936170 0.8085106 0.7872340
head(prevalence(nec.only, detection=1/100, sort=T, count=T))
 # ASV1 ASV12 ASV36  ASV3  ASV5 ASV17
 #   47    43    42    42    38    37
#solar only
sol.only <- subset_samples(ps.dat.trans, tortoise_group=="solar")
pdf("figs/core_sol.pdf")
plot_core(sol.only, plot.type="heatmap", prevalences=prevalences, detections=detections, colours=rev(brewer.pal(5, "Spectral")), min.prevalence=.2, horizontal=T)
dev.off()
head(prevalence(sol.only, detection=1/100, sort=T))
#      ASV1    ASV111     ASV71    ASV135     ASV69     ASV91
# 1.0000000 0.9032258 0.8709677 0.8387097 0.8064516 0.8064516
head(prevalence(sol.only, detection=1/100, sort=T, count=T))
# ASV1 ASV111  ASV71 ASV135  ASV69  ASV91
#   31     28     27     26     25     25

#########################
#ALPHA DIVERSITY BY GROUP
#########################
pdf("figs/adiv_group.pdf")
plot_richness(ps.dat, x="tortoise_group", measures=c("Observed", "Shannon")) + theme_minimal() + geom_boxplot(alpha=0.0)
dev.off()
adiv <- estimate_richness(ps.dat)
# get p value for solar vs necropsy
wilcox.test(adiv[grepl("N", rownames(adiv)),]$Observed, adiv[grepl("S", rownames(adiv)),]$Observed)

# 	Wilcoxon rank sum test with continuity correction

# data:  adiv[grepl("N", rownames(adiv)), ]$Observed and adiv[grepl("S", rownames(adiv)), ]$Observed
# W = 195.5, p-value = 5.401e-08
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(adiv[grepl("N", rownames(adiv)),]$Shannon, adiv[grepl("S", rownames(adiv)),]$Shannon)

# 	Wilcoxon rank sum test

# data:  adiv[grepl("N", rownames(adiv)), ]$Shannon and adiv[grepl("S", rownames(adiv)), ]$Shannon
# W = 290, p-value = 3.195e-06
# alternative hypothesis: true location shift is not equal to 0

################
#BETA DIVERSITY
################
philr.dat <- transform_sample_counts(ps.dat, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
otu.table <- otu_table(philr.dat)
phy_tree(philr.dat) <- makeNodeLabel(phy_tree(philr.dat), method="number", prefix="n")
philr.t <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")
philr.dist <- dist(philr.t, method="euclidean")
pca <- prcomp(as.matrix(philr.dist))
pdf("figs/philr_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("figs/philr_pca.pdf")
fviz_pca_ind(pca, habillage=rawmetadata$tortoise_group, geom="point", pointshape=19, point.size=1, mean.point=F) + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(c(-100,110)) + ylim(c(-100,110))
dev.off()

ps.rare <- rarefy_even_depth(ps.dat, sample.size=15000, rngseed=456)
wuni <- UniFrac(ps.rare, weighted=T)
pca <- prcomp(wuni)
pdf("figs/wunifrac_pca.pdf")
fviz_pca_ind(pca, habillage=rawmetadata$tortoise_group, geom="point", pointshape=19, point.size=1, mean.point=F) + scale_color_brewer(palette="Dark2") + theme_minimal()
dev.off()

#############
#RAREFACTION
#############
pdf("figs/rarefaction_curve.pdf")
rarecurve(seqtab.filtered, label=FALSE)
dev.off()

#######
#DESEQ
#######
ps.deseq <- phyloseq_to_deseq2(ps.dat, ~ Mycoplasmatacea)
desq <- DESeq(ps.deseq, test="Wald", fitType="parametric")
res <- results(desq, cooksCutoff=F)
sigtab <- res[which(res$padj < 0.01),]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps.dat)[rownames(sigtab),], "matrix"))

x <- tapply(sigtab$log2FoldChange, sigtab$V3, function(x) max(x))
x <- sort(x, TRUE)
sigtab$V3 <- factor(as.character(sigtab$V3), levels=names(x))


x <- tapply(sigtab$log2FoldChange, sigtab$V6, function(x) max(x))
x <- sort(x, TRUE)
sigtab$V6 <- factor(as.character(sigtab$V6), levels=names(x))
ggplot(sigtab, aes(x=V6, y=log2FoldChange, color=V3)) + geom_point(size=6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + coord_flip()

#necropsy only
sub <- subset_samples(ps.dat, tortoise_group=="necropsy")
ps.deseq <- phyloseq_to_deseq2(sub, ~ Mycoplasmatacea)
desq <- DESeq(ps.deseq, test="Wald", fitType="parametric")
res <- results(desq, cooksCutoff=F)
sigtab <- res[which(res$padj < 0.01),]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps.dat)[rownames(sigtab),], "matrix"))

x <- tapply(sigtab$log2FoldChange, sigtab$V3, function(x) max(x))
x <- sort(x, TRUE)
sigtab$V3 <- factor(as.character(sigtab$V3), levels=names(x))


x <- tapply(sigtab$log2FoldChange, sigtab$V6, function(x) max(x))
x <- sort(x, TRUE)
sigtab$V6 <- factor(as.character(sigtab$V6), levels=names(x))
ggplot(sigtab, aes(x=V6, y=log2FoldChange, color=V3)) + geom_point(size=6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + coord_flip()

#solar only
sub <- subset_samples(ps.dat, tortoise_group=="solar")
ps.deseq <- phyloseq_to_deseq2(sub, ~ Mycoplasmatacea)
desq <- DESeq(ps.deseq, test="Wald", fitType="parametric")
res <- results(desq, cooksCutoff=F)
sigtab <- res[which(res$padj < 0.01),]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps.dat)[rownames(sigtab),], "matrix"))

x <- tapply(sigtab$log2FoldChange, sigtab$V3, function(x) max(x))
x <- sort(x, TRUE)
sigtab$V3 <- factor(as.character(sigtab$V3), levels=names(x))

x <- tapply(sigtab$log2FoldChange, sigtab$V6, function(x) max(x))
x <- sort(x, TRUE)
sigtab$V6 <- factor(as.character(sigtab$V6), levels=names(x))
ggplot(sigtab, aes(x=V6, y=log2FoldChange, color=V3)) + geom_point(size=6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + coord_flip()

####################
#SOLAR ONLY ANALYSES
####################
ps.dat <- subset_samples(ps.dat, tortoise_group=="solar")
rawmetadata <- sample_data(ps.dat)

#################
#ALPHA DIVERSITY 
#################
#no significant difference between groups
plot_richness(ps.dat, x="urtd_clinical", measures=c("Observed", "Shannon")) + theme_minimal() + geom_boxplot(alpha=0.0)

################
#BETA DIVERSITY
################
philr.dat <- transform_sample_counts(ps.dat, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
otu.table <- otu_table(philr.dat)
phy_tree(philr.dat) <- makeNodeLabel(phy_tree(philr.dat), method="number", prefix="n")
philr.t <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")
philr.dist <- dist(philr.t, method="euclidean")
pca <- prcomp(as.matrix(philr.dist))
#no metadata category clusters
fviz_pca_ind(pca, habillage=rawmetadata$, geom="point", pointshape=19, point.size=1, mean.point=F) + scale_color_brewer(palette="Dark2") + theme_minimal()

ps.rare <- rarefy_even_depth(ps.dat, sample.size=15000, rngseed=456)
wuni <- UniFrac(ps.rare, weighted=T)
pca <- prcomp(wuni)
#no metadata category clusters
fviz_pca_ind(pca, habillage=rawmetadata$urtd_clinical, geom="point", pointshape=19, point.size=1, mean.point=F) + scale_color_brewer(palette="Dark2") + theme_minimal()

#############
#RAREFACTION
#############
rarecurve(otu_table(ps.dat), label=F)

#####################
#INTERACTION NETWORK
#####################
#filter taxa
ps.dat.filt <- filter_taxa(ps.dat, function(x) sum(x > 5) > (0.25*length(x)), T)
spiec.out <- spiec.easi(ps.dat.filt, method='glasso', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(ps.dat.filt)))

pdf("figs/microbe_network.pdf")
plot_network(spiec.graph, ps.dat.filt, type='taxa', color="V3")
dev.off()

#collapse by genus 
ps.g <- tax_glom(ps.dat.filt, taxrank=rank_names(ps.dat.filt)[5], NArm=T, bad_empty=c(NA, "", " ", "\t"))
spiec.out <- spiec.easi(ps.g, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(ps.g)))

#how many positive and negative edges inferred from network?
betaMat <- as.matrix(symBeta(getOptBeta(spiec.out)))
positive <- length(betaMat[betaMat>0])/2
#64
negative <- length(betaMat[betaMat<0])/2
#23
total <- length(betaMat[betaMat!=0])/2
#87

#visualize network with pos and neg edges with different colors
# Col vector up to 74 color samples
col_vector74 = c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2","#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#CCCCCC","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")
asv.ids <- colnames(spiec.out[[1]]$data)
edges <- E(spiec.graph)
edge_cols <- ifelse(betaMat>0, 'black', 'red')[upper.tri(betaMat) & betaMat!=0]
E(spiec.graph)$color <- edge_cols

#How many nodes connected at specific rank
nb_nodes <- vcount(spiec.graph)
tax_table(ps.g) <- tax_table(ps.g)[,getrank]
asv_ids <- V(spiec.graph)$name
idx <- which(row.names(tax_table(ps.g)) %in% asv_ids)
taxa <- as.character(tax_table(ps.g)[,getrank])[idx]

library(intergraph)
library(ggnet)
library(network)

ig2 <- asNetwork(spiec.graph)
network.vertex.names(ig2) <- taxa
net <- ig2
net %v% getrank <- as.character(taxa)
y <- col_vector74[1:nb_nodes]
names(y) <- levels(as.factor(taxa))

#Plot the network
pdf("figs/microbiome_network_genus.pdf")
ggnet2(net, color = getrank, palette = y, size = 6, edge.size=1, edge.color="color", edge.alpha = 0.5, label = FALSE)
dev.off()

#write spiec-easi graph to file to load into cytoscape
write.graph(spiec.graph, file="spieceasi.ncol.txt", format="ncol")
system("sed -i 's/ /,/g' spieceasi.ncol.txt")






