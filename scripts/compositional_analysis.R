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
library(intergraph)
library(ggnet)
library(network)
library(Matrix)

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
ps.dat.filt <- filter_taxa(ps.dat, function(x) sum(x > 1) > (0.75*length(x)), T)
spiec.out <- spiec.easi(ps.dat.filt, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(ps.dat.filt)))

pdf("figs/microbe_network.pdf")
plot_network(spiec.graph, ps.dat.filt, type='taxa', color="V3")
dev.off()

#collapse by genus 
# ps.g <- tax_glom(ps.dat.filt, taxrank=rank_names(ps.dat.filt)[5], NArm=T, bad_empty=c(NA, "", " ", "\t"))
# spiec.out <- spiec.easi(ps.g, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
# spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(ps.g)))

#how many positive and negative edges inferred from network?
betaMat <- as.matrix(symBeta(getOptBeta(spiec.out)))
positive <- length(betaMat[betaMat>0])/2
positive
negative <- length(betaMat[betaMat<0])/2
negative
total <- length(betaMat[betaMat!=0])/2
total

#visualize network with pos and neg edges with different colors
# Col vector up to 74 color samples
getrank <- "V3"
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

ig2 <- asNetwork(spiec.graph)
network.vertex.names(ig2) <- taxa
net <- ig2
net %v% getrank <- as.character(taxa)
y <- col_vector74[1:nb_nodes]
names(y) <- levels(as.factor(taxa))

#Plot the network
pdf("figs/microbiome_network_posneg.pdf")
ggnet2(net, color = getrank, palette = y, size = 6, edge.size=1, edge.color="color", edge.alpha = 0.5, label = FALSE)
dev.off()

#write spiec-easi graph to file 
write.graph(spiec.graph, file="spieceasi.ncol.txt", format="ncol")

#get edge weight distributions
elist.mb <- summary(symBeta(getOptBeta(spiec.out), mode='maxabs'))
pdf("edge_weights.pdf")
hist(elist.mb[,3], xlab="edge weights", main="")
dev.off()

#degree statistics
dd.mb <- degree.distribution(spiec.graph)
pdf("degree_distribution.pdf")
plot(dd.mb, xlab="Degree", main="Degree Distributions", ylab="Frequency", type='b')
dev.off()

#which nodes have the highest number of connections?
sort(igraph::degree(spiec.graph))
# ASV19 ASV319  ASV11 ASV185 ASV172
#      2      2      2      2      2      2      3      3      3      3      3
# ASV326   ASV1 ASV234 ASV347 ASV114  ASV86 ASV202 ASV193 ASV343
#      3      3      3      3      3      3      4      5      5

########################
#NECROPSY ONLY ANALYSES
########################
ps.dat <- subset_samples(ps.dat, tortoise_group=="necropsy")
rawmetadata <- sample_data(ps.dat)

#################
#ALPHA DIVERSITY 
#################
pdf("figs/adiv_necropsy_urtd_clinical.pdf")
plot_richness(ps.dat, x="urtd_clinical", measures=c("Observed", "Shannon")) + theme_minimal() + geom_boxplot(alpha=0.0)
dev.off()

pdf("figs/adiv_necropsy_nasal_discharge.pdf")
plot_richness(ps.dat, x="nasal_discharge", measures=c("Observed", "Shannon")) + theme_minimal() + geom_boxplot(alpha=0.0)
dev.off()

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
pdf("philr_necropsy_urtd_clinical.pdf")
fviz_pca_ind(pca, habillage=rawmetadata$urtd_clinical, geom="point", pointshape=19, point.size=1, mean.point=T) + scale_color_brewer(palette="Dark2") + theme_minimal()
dev.off()

ps.rare <- rarefy_even_depth(ps.dat, sample.size=15000, rngseed=456)
wuni <- UniFrac(ps.rare, weighted=T)
pca <- prcomp(wuni)
pdf("wuni_necropsy_urtd_clinical.pdf")
fviz_pca_ind(pca, habillage=rawmetadata$urtd_clinical, geom="point", pointshape=19, point.size=1, mean.point=T) + scale_color_brewer(palette="Dark2") + theme_minimal()
dev.off()

#####################
#INTERACTION NETWORK
#####################
#filter taxa
ps.dat.filt <- filter_taxa(ps.dat, function(x) sum(x > 1) > (0.25*length(x)), T)
spiec.out <- spiec.easi(ps.dat.filt, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(ps.dat.filt)))

pdf("figs/microbe_network_necropsy.pdf")
plot_network(spiec.graph, ps.dat.filt, type='taxa', color="V3")
dev.off()

#collapse by genus 
# ps.g <- tax_glom(ps.dat.filt, taxrank=rank_names(ps.dat.filt)[5], NArm=T, bad_empty=c(NA, "", " ", "\t"))
# spiec.out <- spiec.easi(ps.g, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
# spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(ps.g)))

#how many positive and negative edges inferred from network?
betaMat <- as.matrix(symBeta(getOptBeta(spiec.out)))
positive <- length(betaMat[betaMat>0])/2
positive
negative <- length(betaMat[betaMat<0])/2
negative
total <- length(betaMat[betaMat!=0])/2
total

#visualize network with pos and neg edges with different colors
# Col vector up to 74 color samples
getrank <- "V3"
col_vector74 = c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2","#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#CCCCCC","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")
asv.ids <- colnames(spiec.out[[1]]$data)
edges <- E(spiec.graph)
edge_cols <- ifelse(betaMat>0, 'black', 'red')[upper.tri(betaMat) & betaMat!=0]
E(spiec.graph)$color <- edge_cols

#How many nodes connected at specific rank
nb_nodes <- vcount(spiec.graph)
tax_table(ps.dat.filt) <- tax_table(ps.dat.filt)[,getrank]
asv_ids <- V(spiec.graph)$name
idx <- which(row.names(tax_table(ps.dat.filt)) %in% asv_ids)
taxa <- as.character(tax_table(ps.dat.filt)[,getrank])[idx]

ig2 <- asNetwork(spiec.graph)
network.vertex.names(ig2) <- taxa
net <- ig2
net %v% getrank <- as.character(taxa)
y <- col_vector74[1:nb_nodes]
names(y) <- levels(as.factor(taxa))

#Plot the network
pdf("figs/microbiome_network_posneg.pdf")
ggnet2(net, color = getrank, palette = y, size = 6, edge.size=1, edge.color="color", edge.alpha = 0.5, label = FALSE)
dev.off()

#write spiec-easi graph to file 
write.graph(spiec.graph, file="spieceasi.ncol.txt", format="ncol")

#get edge weight distributions
elist.mb <- summary(symBeta(getOptBeta(spiec.out), mode='maxabs'))
pdf("edge_weights.pdf")
hist(elist.mb[,3], xlab="edge weights", main="")
dev.off()

#degree statistics
dd.mb <- degree.distribution(spiec.graph)
pdf("degree_distribution.pdf")
plot(dd.mb, xlab="Degree", main="Degree Distributions", ylab="Frequency", type='b')
dev.off()

#which nodes have the highest number of connections?
sort(igraph::degree(spiec.graph))
#  ASV14   ASV8  ASV20  ASV34  ASV88  ASV12  ASV19   ASV7 ASV296  ASV39 ASV128
#      0      0      0      0      0      0      0      0      0      0      0
#   ASV4  ASV70 ASV112  ASV31  ASV28  ASV58 ASV134 ASV132   ASV6 ASV367  ASV41
#      0      0      0      0      0      0      0      0      0      0      1
# ASV117 ASV447  ASV92  ASV79  ASV37  ASV11  ASV78 ASV303 ASV278  ASV36   ASV5
#      1      1      1      1      1      1      1      1      1      1      1
# ASV162  ASV53 ASV100 ASV105   ASV9  ASV48 ASV151  ASV95 ASV173 ASV211 ASV195
#      1      1      1      1      1      1      1      2      2      2      2
#  ASV40  ASV35 ASV167  ASV21  ASV60 ASV111 ASV244 ASV148   ASV3 ASV135   ASV2
#      2      2      2      2      2      2      2      3      3      3      3
#  ASV32 ASV183 ASV141 ASV170 ASV200 ASV116  ASV18 ASV137 ASV468 ASV292 ASV226
#      3      3      3      3      3      3      4      4      4      4      4
#   ASV1 ASV219  ASV69  ASV96 ASV114 ASV185 ASV146 ASV191  ASV71  ASV17 ASV171
#      4      4      4      4      4      5      5      5      5      6      7
# ASV207  ASV86  ASV91
#      8     10     12

#Strength -- weighted measure of degree, total number of interactions of each ASV with anyother ASV
sort(igraph::strength(spiec.graph))
#  ASV14   ASV8  ASV20  ASV34  ASV88  ASV12  ASV19   ASV7 ASV296  ASV39 ASV128
#      0      0      0      0      0      0      0      0      0      0      0
#   ASV4  ASV70 ASV112  ASV31  ASV28  ASV58 ASV134 ASV132   ASV6 ASV367  ASV41
#      0      0      0      0      0      0      0      0      0      0      1
# ASV117 ASV447  ASV92  ASV79  ASV37  ASV11  ASV78 ASV303 ASV278  ASV36   ASV5
#      1      1      1      1      1      1      1      1      1      1      1
# ASV162  ASV53 ASV100 ASV105   ASV9  ASV48 ASV151  ASV95 ASV173 ASV211 ASV195
#      1      1      1      1      1      1      1      2      2      2      2
#  ASV40  ASV35 ASV167  ASV21  ASV60 ASV111 ASV244 ASV148   ASV3 ASV135   ASV2
#      2      2      2      2      2      2      2      3      3      3      3
#  ASV32 ASV183 ASV141 ASV170 ASV200 ASV116  ASV18 ASV137 ASV468 ASV292 ASV226
#      3      3      3      3      3      3      4      4      4      4      4
#   ASV1 ASV219  ASV69  ASV96 ASV114 ASV185 ASV146 ASV191  ASV71  ASV17 ASV171
#      4      4      4      4      4      5      5      5      5      6      7
# ASV207  ASV86  ASV91
#      8     10     12

#Closeness -- how many steps are required to access every other node from a given node -- higher values mean less centrality
sort(igraph::closeness(spiec.graph, normalized=T))
#      ASV14       ASV8      ASV20      ASV34      ASV88      ASV12      ASV19
# 0.01250000 0.01250000 0.01250000 0.01250000 0.01250000 0.01250000 0.01250000
#       ASV7     ASV296      ASV39     ASV128       ASV4      ASV70     ASV112
# 0.01250000 0.01250000 0.01250000 0.01250000 0.01250000 0.01250000 0.01250000
#      ASV31      ASV28      ASV58     ASV134     ASV132       ASV6     ASV367
# 0.01250000 0.01250000 0.01250000 0.01250000 0.01250000 0.01250000 0.01250000
#      ASV53      ASV48     ASV162      ASV79       ASV5       ASV9     ASV111
# 0.01265823 0.01265823 0.01331760 0.01332209 0.01332209 0.01332209 0.01332659
#       ASV1      ASV78     ASV117      ASV41      ASV40     ASV105     ASV173
# 0.01333108 0.03045490 0.03077522 0.03093187 0.03104126 0.03105346 0.03137411
#     ASV151     ASV100     ASV278     ASV303      ASV92      ASV35     ASV244
# 0.03137411 0.03141153 0.03144904 0.03148665 0.03149920 0.03149920 0.03153693
#      ASV36      ASV95      ASV60     ASV195      ASV21     ASV135      ASV37
# 0.03157474 0.03158737 0.03160000 0.03162530 0.03162530 0.03166333 0.03167602
#       ASV2     ASV141       ASV3     ASV226     ASV167     ASV211     ASV447
# 0.03167602 0.03168873 0.03177796 0.03179074 0.03180354 0.03182917 0.03184200
#      ASV11     ASV170     ASV200     ASV116     ASV183     ASV148      ASV32
# 0.03185484 0.03185484 0.03188055 0.03190630 0.03194501 0.03197086 0.03198381
#     ASV185     ASV292      ASV96     ASV114     ASV468     ASV219      ASV71
# 0.03199676 0.03203569 0.03203569 0.03206169 0.03207471 0.03210077 0.03211382
#     ASV137      ASV69     ASV191      ASV18      ASV17     ASV146     ASV171
# 0.03212688 0.03212688 0.03213995 0.03220546 0.03220546 0.03223174 0.03231084
#      ASV86     ASV207      ASV91
# 0.03235053 0.03248355 0.03249691

#eigenvector centrality -- measure of being well connected to the well connected
sort(igraph::eigen_centrality(spiec.graph)$vector)
#       ASV162        ASV48        ASV53         ASV9        ASV79         ASV5
# 0.0005916426 0.0006460567 0.0006460680 0.0007076142 0.0007076288 0.0007076491
#       ASV111         ASV1        ASV14         ASV8        ASV20        ASV34
# 0.0013084042 0.0022318107 0.0027634920 0.0027634920 0.0027634920 0.0027634920
#        ASV88        ASV12        ASV19         ASV7       ASV296        ASV39
# 0.0027634920 0.0027634920 0.0027634920 0.0027634920 0.0027634920 0.0027634920
#       ASV128         ASV4        ASV70       ASV112        ASV31        ASV28
# 0.0027634920 0.0027634920 0.0027634920 0.0027634920 0.0027634920 0.0027634920
#        ASV58       ASV134       ASV132         ASV6       ASV367        ASV78
# 0.0027634920 0.0027634920 0.0027634920 0.0027634920 0.0027634920 0.0059168780
#       ASV117        ASV41       ASV105        ASV40       ASV100       ASV173
# 0.0101639565 0.0116553150 0.0150083645 0.0282728390 0.0403865510 0.0501885341
#       ASV303       ASV244        ASV95       ASV151        ASV60        ASV35
# 0.0509160121 0.0582050578 0.0645103730 0.0676910959 0.0682981306 0.0727882784
#       ASV135       ASV278         ASV2       ASV195        ASV92       ASV141
# 0.0759090808 0.0857029265 0.0890407616 0.0901894903 0.1024611261 0.1088348243
#         ASV3        ASV37        ASV36        ASV21       ASV170       ASV447
# 0.1205043504 0.1272651927 0.1296302591 0.1355472293 0.1551474670 0.1636619932
#        ASV11       ASV200        ASV96       ASV211       ASV226       ASV148
# 0.1818743290 0.1958439795 0.2078328105 0.2113234553 0.2425081573 0.2450031394
#       ASV137       ASV167        ASV71        ASV32       ASV191       ASV183
# 0.2520007103 0.2611255498 0.2629189657 0.2843834281 0.2863340291 0.2889438333
#       ASV116       ASV185       ASV292       ASV219       ASV114       ASV468
# 0.3078803126 0.3485881439 0.3601209380 0.3671420944 0.4037313756 0.4419018522
#        ASV18        ASV69       ASV171        ASV17       ASV146       ASV207
# 0.4542236676 0.5284693722 0.6580636072 0.6693846002 0.6857502011 0.8467137755
#        ASV91        ASV86
# 0.9411880254 1.0000000000

#which asvs are connected to the Mycoplasma asvs?
igraph::neighbors(spiec.graph, v=which(V(spiec.graph)$name=="ASV3"))
# + 3/80 vertices, named, from b304816:
# [1] ASV18  ASV135 ASV2
# ASV2 is corynebacterium, both are possible pneumonia causing

#define communities in the graph
igraph::components(spiec.graph)
# $membership
#  ASV14   ASV8  ASV20  ASV34  ASV41  ASV88  ASV18  ASV12  ASV19   ASV7  ASV95
#      1      2      3      4      5      6      5      7      8      9      5
# ASV173 ASV148 ASV117 ASV207   ASV3 ASV447 ASV296 ASV211 ASV195  ASV92  ASV79
#      5      5      5      5      5      5     10      5      5      5     11
#  ASV37  ASV39  ASV40 ASV128   ASV4  ASV11  ASV78  ASV70 ASV112  ASV35  ASV17
#      5     12      5     13     14      5      5     15     16      5      5
# ASV303 ASV135 ASV167 ASV137 ASV468 ASV278  ASV36  ASV31   ASV2  ASV32  ASV21
#      5      5      5      5      5      5      5     17      5      5      5
# ASV185 ASV292   ASV5 ASV183 ASV146 ASV226  ASV60   ASV1  ASV28  ASV58 ASV134
#      5      5     11      5      5      5      5     11     18     19     20
# ASV132 ASV162  ASV53 ASV100   ASV6 ASV111 ASV141 ASV219  ASV91 ASV191  ASV71
#     21     11     22      5     23     11      5      5      5      5      5
#  ASV69 ASV170 ASV171  ASV96 ASV105   ASV9 ASV367 ASV244 ASV200  ASV48 ASV114
#      5      5      5      5      5     11     24      5      5     22      5
# ASV116  ASV86 ASV151
#      5      5      5

# $csize
#  [1]  1  1  1  1 51  1  1  1  1  1  6  1  1  1  1  1  1  1  1  1  1  2  1  1

# $no
# [1] 24

#within largest cluster are there subpopulations?
giant <- igraph::decompose(spiec.graph)[[5]]
comm <- cluster_infomap(giant)
modularity(comm)
# [1] 0.4849108
#plot communities
pdf("figs/network_subcommunity_node_color.pdf")
par(mar=c(0,0,0,0));plot(comm, giant)
dev.off()
#color nodes by sub community
V(giant)$color <- membership(comm)
pdf("figs/network_subcommunity_elipse_color.pdf")
par(mar=c(0,0,0,0));plot(giant)
dev.off()

#k core decomposition -- define what asvs are the core the the network vs the periphery
#core network
which(coreness(spiec.graph)==3)
# ASV207  ASV17 ASV468 ASV185 ASV292 ASV183 ASV146 ASV226 ASV219  ASV91  ASV69
#     15     33     38     45     46     48     49     50     63     64     67
# ASV171 ASV114 ASV116  ASV86
#     69     77     78     79

which(coreness(spiec.graph)==2)
 # ASV18  ASV95 ASV148   ASV3 ASV211 ASV195  ASV35 ASV135 ASV167 ASV137   ASV2
 #     7     11     13     16     19     20     32     35     36     37     42
 # ASV32  ASV60 ASV141 ASV191  ASV71 ASV170  ASV96 ASV200
 #    43     51     62     65     66     68     70     75

which(coreness(spiec.graph)==1)
#  ASV41 ASV173 ASV117 ASV447  ASV92  ASV79  ASV37  ASV40  ASV11  ASV78 ASV303
#      5     12     14     17     21     22     23     25     28     29     34
# ASV278  ASV36  ASV21   ASV5   ASV1 ASV162  ASV53 ASV100 ASV111 ASV105   ASV9
#     39     40     44     47     52     57     58     59     61     71     72
# ASV244  ASV48 ASV151
#     74     76     80

#periphery of the network
which(coreness(spiec.graph)==0)
 # ASV14   ASV8  ASV20  ASV34  ASV88  ASV12  ASV19   ASV7 ASV296  ASV39 ASV128
 #     1      2      3      4      6      8      9     10     18     24     26
 #  ASV4  ASV70 ASV112  ASV31  ASV28  ASV58 ASV134 ASV132   ASV6 ASV367
 #    27     30     31     41     53     54     55     56     60     73

# Visualizing network structure
V(spiec.graph)$coreness <- coreness(spiec.graph)
par(mfrow=c(1, 3), mar=c(0.1,0.1,1,0.1))
set.seed(777); fr <- layout_with_fr(spiec.graph)
pdf("figs/core_network_structure.pdf")
for (k in 1:3){
  V(spiec.graph)$color <- ifelse(V(spiec.graph)$coreness>=k, "orange", "grey")
  plot(spiec.graph, main=paste0(k, '-core shell'), layout=fr)
}
dev.off()

#############
#RAREFACTION
#############
pdf("figs/rarefaction_curve_necropsy.pdf")
rarecurve(otu_table(ps.dat), label=F)
dev.off()

#picrust
#most connected core asvs functional analysis
system("awk 'NF{NF-=2}1' FS='\t' OFS='\t' sequence_taxonomy_table.16s.merged.txt > sequence_table.forpicrust.txt")
system("picrust2_pipeline.py -s core3.fa -i sequence_table.forpicrust.txt -o picrust -p 6")

