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
library(ranacapa)

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

#necropsy only
nec.only <- subset_samples(ps.dat.trans, tortoise_group=="necropsy")
head(prevalence(nec.only, detection=1/100, sort=T))
#      ASV1     ASV12     ASV36      ASV3      ASV5     ASV17
# 1.0000000 0.9148936 0.8936170 0.8936170 0.8085106 0.7872340
head(prevalence(nec.only, detection=1/100, sort=T, count=T))
 # ASV1 ASV12 ASV36  ASV3  ASV5 ASV17
 #   47    43    42    42    38    37
#solar only
sol.only <- subset_samples(ps.dat.trans, tortoise_group=="solar")
head(prevalence(sol.only, detection=1/100, sort=T))
#      ASV1    ASV111     ASV71    ASV135     ASV69     ASV91
# 1.0000000 0.9032258 0.8709677 0.8387097 0.8064516 0.8064516
head(prevalence(sol.only, detection=1/100, sort=T, count=T))
# ASV1 ASV111  ASV71 ASV135  ASV69  ASV91
#   31     28     27     26     25     25

#################
#ALPHA DIVERSITY 
#################
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

ps.nec <- subset_samples(ps.dat, tortoise_group=="necropsy")

pdf("figs/adiv_necropsy_urtd_clinical.pdf")
plot_richness(ps.nec, x="urtd_clinical", measures=c("Observed", "Shannon")) + theme_minimal() + geom_boxplot(alpha=0.0)
dev.off()

# get p value for urtd clinical within dtcc vs non
ps.nec.p <- subset_samples(ps.nec, urtd_clinical=="1")
ps.nec.n <- subset_samples(ps.nec, urtd_clinical=="0")

adiv.nec.p <- estimate_richness(ps.nec.p)
adiv.nec.n <- estimate_richness(ps.nec.n)

wilcox.test(adiv.nec.p$Observed, adiv.nec.n$Observed)

# 	Wilcoxon rank sum test with continuity correction

# data:  adiv.nec.p$Observed and adiv.nec.n$Observed
# W = 80.5, p-value = 0.001641
# alternative hypothesis: true location shift is not equal to 0

# Warning message:
# In wilcox.test.default(adiv.nec.p$Observed, adiv.nec.n$Observed) :
#   cannot compute exact p-value with ties

wilcox.test(adiv.nec.p$Shannon, adiv.nec.n$Shannon)

# 	Wilcoxon rank sum test

# data:  adiv.nec.p$Shannon and adiv.nec.n$Shannon
# W = 90, p-value = 0.002683
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
pdf("figs/philr_screeplot_group.pdf")
screeplot(pca)
dev.off()
pdf("figs/philr_pca_group.pdf")
fviz_pca_ind(pca, habillage=rawmetadata$tortoise_group, geom="point", pointshape=19, point.size=1, mean.point=F) + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(c(-100,110)) + ylim(c(-100,110))
dev.off()

ps.rare <- rarefy_even_depth(ps.dat, sample.size=15000, rngseed=456)
wuni <- UniFrac(ps.rare, weighted=T)
pca <- prcomp(wuni)
pdf("figs/wunifrac_pca_group.pdf")
fviz_pca_ind(pca, habillage=rawmetadata$tortoise_group, geom="point", pointshape=19, point.size=1, mean.point=F) + scale_color_brewer(palette="Dark2") + theme_minimal()
dev.off()

#permanova
metadata <- as(sample_data(ps.dat), "data.frame")
adonis(philr.dist ~ tortoise_group, data=metadata)

# Call:
# adonis(formula = philr.dist ~ tortoise_group, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# tortoise_group  1    4518.2  4518.2  15.226 0.16691  0.001 ***
# Residuals      76   22551.8   296.7         0.83309
# Total          77   27070.0                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#DTCC tortoises only
philr.dat <- transform_sample_counts(ps.nec, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
otu.table <- otu_table(philr.dat)
phy_tree(philr.dat) <- makeNodeLabel(phy_tree(philr.dat), method="number", prefix="n")
philr.t <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")
philr.dist <- dist(philr.t, method="euclidean")
pca <- prcomp(as.matrix(philr.dist))

rawmetadata.nec <- sample_data(ps.nec)

pdf("figs/philr_necropsy_urtd_clinical.pdf")
fviz_pca_ind(pca, habillage=rawmetadata.nec$urtd_clinical, geom="point", pointshape=19, point.size=1, mean.point=T) + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(c(-75,110)) + ylim(c(-75,110))
dev.off()

ps.rare <- rarefy_even_depth(ps.nec, sample.size=15000, rngseed=456)
wuni <- UniFrac(ps.rare, weighted=T)
pca <- prcomp(wuni)
pdf("figs/wuni_necropsy_urtd_clinical.pdf")
fviz_pca_ind(pca, habillage=rawmetadata.nec$urtd_clinical, geom="point", pointshape=19, point.size=1, mean.point=T) + scale_color_brewer(palette="Dark2") + theme_minimal()
dev.off()

metadata <- as(sample_data(ps.nec), "data.frame")
adonis(philr.dist ~ urtd_clinical, data=metadata)

# Call:
# adonis(formula = philr.dist ~ urtd_clinical, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#               Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
# urtd_clinical  1    3108.4 3108.41  6.5936 0.1278  0.002 **
# Residuals     45   21214.1  471.42         0.8722
# Total         46   24322.5                 1.0000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ sex, data=metadata)

# Call:
# adonis(formula = philr.dist ~ sex, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# sex        2    1253.2  626.61  1.1951 0.05152  0.276
# Residuals 44   23069.3  524.30         0.94848
# Total     46   24322.5                 1.00000

adonis(philr.dist ~ agegroup, data=metadata)

# Call:
# adonis(formula = philr.dist ~ agegroup, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# agegroup   2    3179.6 1589.78  3.3084 0.13073  0.001 ***
# Residuals 44   21143.0  480.52         0.86927
# Total     46   24322.5                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#############
#RAREFACTION
#############
pdf("figs/rarefaction_curve.pdf")
ggrare(ps.dat, color="tortoise_group") + theme_minimal()
dev.off()

#######
#DESEQ
#######
nec.only <- subset_samples(ps.dat, tortoise_group=="necropsy")
ps.deseq <- phyloseq_to_deseq2(nec.only, ~ urtd_clinical)
desq <- DESeq(ps.deseq, test="Wald", fitType="parametric")
res <- results(desq, cooksCutoff=F)
sigtab <- res[which(res$padj < 0.01),]
dim(sigtab)
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps.dat)[rownames(sigtab),], "matrix"))

x <- tapply(sigtab$log2FoldChange, sigtab$V3, function(x) max(x))
x <- sort(x, TRUE)
sigtab$V3 <- factor(as.character(sigtab$V3), levels=names(x))

x <- tapply(sigtab$log2FoldChange, sigtab$V7, function(x) max(x))
x <- sort(x, TRUE)
sigtab$V7 <- factor(as.character(sigtab$V7), levels=names(x))
pdf("figs/deseq_necropsy_urtd_clinical.pdf")
ggplot(sigtab, aes(x=V7, y=log2FoldChange, color=V3)) + geom_point(size=6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + coord_flip() + theme_minimal()
dev.off()

#####################
#INTERACTION NETWORK
#####################
#filter taxa
ps.dat.filt <- filter_taxa(ps.dat, function(x) sum(x > 10) > (0.05*length(x)), T)
# spiec.out <- spiec.easi(ps.dat, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
# spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(ps.dat)))

#collapse by genus 
ps.g <- tax_glom(ps.dat.filt, taxrank=rank_names(ps.dat.filt)[5], NArm=T, bad_empty=c(NA, "", " ", "\t"))
spiec.out <- spiec.easi(ps.g, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
spiec.graph <- adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(ps.g)))

pdf("figs/microbe_network.pdf")
plot_network(spiec.graph, ps.dat.filt, type='taxa', color="V3")
dev.off()

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
tax_table(ps.dat) <- tax_table(ps.dat)[,getrank]
asv_ids <- V(spiec.graph)$name
idx <- which(row.names(tax_table(ps.dat)) %in% asv_ids)
taxa <- as.character(tax_table(ps.dat)[,getrank])[idx]

ig2 <- asNetwork(spiec.graph)
network.vertex.names(ig2) <- taxa
net <- ig2
net %v% getrank <- as.character(taxa)
y <- col_vector74[1:nb_nodes]
names(y) <- levels(as.factor(taxa))

#Plot the network
pdf("figs/microbiome_network_posneg.pdf")
ggnet2(net, color = getrank, alpha=0.75, palette = y, size = 6, edge.size=1, edge.color="color", edge.alpha = 0.5, label = FALSE)
dev.off()

#write spiec-easi graph to file 
write.graph(spiec.graph, file="spieceasi.ncol.txt", format="ncol")

#get edge weight distributions
elist.mb <- summary(symBeta(getOptBeta(spiec.out), mode='maxabs'))
pdf("figs/edge_weights.pdf")
hist(elist.mb[,3], xlab="edge weights", main="")
dev.off()

#degree statistics
dd.mb <- degree.distribution(spiec.graph)
pdf("figs/degree_distribution.pdf")
plot(dd.mb, xlab="Degree", main="Degree Distributions", ylab="Frequency", type='b')
dev.off()

sort(igraph::degree(spiec.graph))

#Closeness -- how many steps are required to access every other node from a given node -- higher values mean less centrality
sort(igraph::closeness(spiec.graph, normalized=T))

#eigenvector centrality -- measure of being well connected to the well connected
sort(igraph::eigen_centrality(spiec.graph)$vector)

#which asvs are connected to the Mycoplasma asvs?
igraph::neighbors(spiec.graph, v=which(V(spiec.graph)$name=="ASV3"))
# + 3/212 vertices, named, from 2ec4094:
# [1] ASV135 ASV13  ASV69
igraph::neighbors(spiec.graph, v=which(V(spiec.graph)$name=="ASV240"))
# + 6/212 vertices, named, from 2ec4094:
# [1] ASV573 ASV186 ASV848 ASV223 ASV446 ASV781
igraph::neighbors(spiec.graph, v=which(V(spiec.graph)$name=="ASV1682"))
# + 6/212 vertices, named, from 2ec4094:
# [1] ASV573 ASV186 ASV848 ASV223 ASV446 ASV781

#define communities in the graph
igraph::components(spiec.graph)

#within largest cluster are there subpopulations?
giant <- igraph::decompose(spiec.graph)[[1]]
comm <- cluster_infomap(giant)
modularity(comm)
#plot communities
pdf("figs/network_subcommunity_elipse_color.pdf")
par(mar=c(0,0,0,0));plot(comm, giant)
dev.off()
#color nodes by sub community
V(giant)$color <- membership(comm)
pdf("figs/network_subcommunity_node_color.pdf")
par(mar=c(0,0,0,0));plot(giant)
dev.off()

#k core decomposition -- define what asvs are the core the the network vs the periphery
#core network
which(coreness(spiec.graph)==4)
 # ASV455  ASV471  ASV207  ASV254  ASV258  ASV295   ASV17   ASV93 ASV1995  ASV226
 #      6      39      43      50      52      60      83     120     131     134
 # ASV276  ASV269 ASV2290  ASV285  ASV243  ASV248   ASV60  ASV174  ASV286  ASV487
 #    138     139     141     142     143     144     147     189     190     191
 # ASV318   ASV69  ASV261  ASV151
 #    194     195     208     210

#periphery of the network
which(coreness(spiec.graph)==0)
   # ASV7 ASV1358  ASV296  ASV559   ASV27 ASV1021  ASV973  ASV356
   #   27      71      76      92     161     162     174     177

# Visualizing network structure
V(spiec.graph)$coreness <- coreness(spiec.graph)
par(mfrow=c(1, 3), mar=c(0.1,0.1,1,0.1))
set.seed(777); fr <- layout_with_fr(spiec.graph)

for (k in 1:3){
  V(spiec.graph)$color <- ifelse(V(spiec.graph)$coreness>=k, "orange", "grey")
  plot(spiec.graph, main=paste0(k, '-core shell'), layout=fr)
}


