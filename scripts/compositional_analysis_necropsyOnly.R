# Compositional data analysis
# Tortoise dataset

###########
#LIBRARIES
###########
library(tidyverse)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(dendextend)
library(compositions)
library(ape)
library(factoextra)
library(RColorBrewer)
library(philr)
library(phyloseq)
library(reshape2)
library(vegan)
library(phylofactor)
library(ggtree)
library(DESeq2)
library(igraph)
library(SpiecEasi)
library(beepr)
library(corrplot)

##########
#ENV SETUP
##########
PATH="/Users/mann/github/tortoise"
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
seqtab.filtered <- read.table("sequence_table.16s.filtered.txt", header=T, row.names=1)
taxa <- read.table("assigntax/rep_set_fix_tax_assignments_6levels.txt", header=F, sep="\t", row.names=1)
notinmeta <- setdiff(row.names(seqtab.filtered), rawmetadata$sample_id) #record samples absent in either metadata or OTU table
notinmeta
notinraw <- setdiff(rawmetadata$sample_id, row.names(seqtab.filtered))
notinraw
tree <- read.tree("rep_set.root.tre")

#filter out samples not in metadata
seqtab.filtered <- subset(seqtab.filtered, !(row.names(seqtab.filtered) %in% notinmeta))
dim(seqtab.filtered)
dim(rawmetadata)
rownames(rawmetadata) <- rawmetadata$sample_id 

ps.dat <- phyloseq(otu_table(seqtab.filtered, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(as.matrix(taxa)), tree)

ps.dat <- subset_samples(ps.dat, tortoise_group=="necropsy")
rawmetadata <- rawmetadata[rawmetadata$tortoise_group == "necropsy",]
rownames(rawmetadata) <- rawmetadata$sample_id 

##############################
#PHILR DISTANCE (PHYLOGENETIC)
##############################
philr.dat <- transform_sample_counts(ps.dat, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
is.rooted(phy_tree(philr.dat)) #check that tree is rooted
is.binary.tree(phy_tree(philr.dat)) #check that multichotomies are resolved in tree
phy_tree(philr.dat) <- makeNodeLabel(phy_tree(philr.dat), method="number", prefix="n")
otu.table <- otu_table(philr.dat)
tree <- phy_tree(philr.dat)
metadata <- sample_data(philr.dat)
tax <- tax_table(philr.dat)
philr.t <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")

# PCA
philr.dist <- dist(philr.t, method="euclidean")
pca <- prcomp(as.matrix(philr.dist))
pdf("figs/philr_screeplot_necOnly.pdf")
screeplot(pca)
dev.off()
pdf("figs/philr_pca_necOnly.pdf")
fviz_pca_ind(pca, habillage=rawmetadata$mycoplasma_100reads, geom="point", addEllipses=F, pointshape=19, point.size=1, mean.point=T, repel=T) + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-110,160) + ylim(-110,160)
dev.off()

# PCoA of Weighted UniFrac for good measure
ps.rare <- rarefy_even_depth(ps.dat, rngseed=778)
wuni <- UniFrac(ps.rare, weighted=T)
hc <- hclust(wuni, method="complete")
df2 <- data.frame(cluster=cutree(hc,10), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
p1 <- ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
merge <- merge(df2, rawmetadata, by.x=c("states"), by.y=c("sample_id"))
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$cluster))) + geom_tile() + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
pdf("figs/wunifrac_dendrogram_necOnly.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()

#pca plot of rarefied data
pca <- prcomp(wuni)
pdf("figs/wunifrac_screeplot_necOnly.pdf")
screeplot(pca)
dev.off()
pdf("figs/wunifrac_pca_necOnly.pdf")
fviz_pca_ind(pca, habillage=merge$tortoise_group, geom="point", pointsize=3, mean.point=FALSE) + labs(title="Weighted Unifrac Distance") + theme_classic()
dev.off()

#############################
#DESEQ DIFFERENTIAL ABUNDACE
#############################
ps.deseq <- phyloseq_to_deseq2(ps.dat, ~mycoplasma_100reads)
dseq.results <- DESeq(ps.deseq, test="Wald", fitType="parametric")
# order by p value, filter out na entries
res <- results(dseq.results, cooksCutoff=F)
sigtab <- as.data.frame(res[which(res$padj < 0.001),])
dim(sigtab)
# merge with your taxonomy to plot log fold change 
sigtab <- merge(sigtab, taxa, by=0, all.x=T)
dim(sigtab)
write.table(sigtab, as.character("deseq_results.txt"), row.names=F, sep="\t", quote=F)
# plot ordered by phylum
x <- tapply(sigtab$log2FoldChange, sigtab$V3, function(x) max(x))
x <- sort(x, TRUE)
sigtab$V3 = factor(as.character(sigtab$V3), levels=names(x))
# then order by genus
x <- tapply(sigtab$log2FoldChange, sigtab$V7, function(x) max(x))
x <- sort(x, TRUE)
sigtab$V7 <- factor(as.character(sigtab$V7), levels=names(x))
pdf("figs/logfold_change_necOnly.pdf")
ggplot(sigtab, aes(x=V7, y=log2FoldChange, color=V3)) + geom_point(size=3) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + coord_flip() + theme_minimal()
dev.off()

################
#ALPHA DIVERSITY 
################
pdf("figs/adiv_necOnly.pdf")
plot_richness(ps.dat, x="mycoplasma_100reads", measures=c("Observed", "Shannon")) + theme_minimal() + geom_boxplot(alpha=0.0)
dev.off()
adiv <- estimate_richness(ps.dat)
adiv <- merge(adiv, rawmetadata, by=0)
wilcox.test(adiv[adiv$mycoplasma_100reads == "yes",]$Observed, adiv[adiv$mycoplasma_100reads == "no",]$Observed)

# 	Wilcoxon rank sum test with continuity correction

# data:  adiv[adiv$mycoplasma_100reads == "yes", ]$Observed and adiv[adiv$mycoplasma_100reads == "no", ]$Observed
# W = 168, p-value = 0.03451
# alternative hypothesis: true location shift is not equal to 0
wilcox.test(adiv[adiv$mycoplasma_100reads == "yes",]$Shannon, adiv[adiv$mycoplasma_100reads == "no",]$Shannon)

# 	Wilcoxon rank sum test

# data:  adiv[adiv$mycoplasma_100reads == "yes", ]$Shannon and adiv[adiv$mycoplasma_100reads == "no", ]$Shannon
# W = 208, p-value = 0.2149
# alternative hypothesis: true location shift is not equal to 0

##########
#PERMANOVA
##########
# is there a difference in microbial diversity across sample type?
metadata <- as(sample_data(ps.dat), "data.frame")
adonis(philr.dist ~ mycoplasma_100reads, data=metadata)

# Call:
# adonis(formula = philr.dist ~ mycoplasma_100reads, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# mycoplasma_100reads  1    5068.3  5068.3  11.655 0.20572  0.001 ***
# Residuals           45   19568.5   434.9         0.79428
# Total               46   24636.8                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##############
#PHYLOFACTOR
##############
OTUTable <- as.matrix(t(otu_table(ps.dat)))
filt.list <- colnames(OTUTable)
filtmap <- rawmetadata[rawmetadata$sample_id %in% filt.list,]
filtmap <- filtmap[match(filt.list, filtmap$sample_id),]
x <- as.factor(filtmap$mycoplasma_100reads) # CHANGE ME to your variable of interest
tree <- phy_tree(philr.dat)
tax <- read.table("assigntax/rep_set_fix_tax_assignments_phylofactor.txt", sep="\t", header=T)
common.otus <- which(rowSums(OTUTable>0)>10)
OTUTable <- OTUTable[common.otus,]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(OTUTable)))
PF <- PhyloFactor(OTUTable, tree, x, nfactors=3)
PF$Data <- PF$Data[PF$tree$tip.label,]
gtree <- pf.tree(PF,layout="rectangular")
pdf("figs/phylofactor_tree_rectangular_necOnly.pdf")
gtree$ggplot + geom_tiplab(size=2, align=T)
dev.off()


y <- t(PF$basis[,1]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
pdf("figs/factor1_boxp_necOnly.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[1]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 1')
dev.off()
wilcox.test(dat[dat$V1 == "yes",]$V2, dat[dat$V1 == "no",]$V2)

# 	Wilcoxon rank sum test

# data:  dat[dat$V1 == "yes", ]$V2 and dat[dat$V1 == "no", ]$V2
# W = 528, p-value = 3.442e-12
# alternative hypothesis: true location shift is not equal to 0

y <- t(PF$basis[,2]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
pdf("figs/factor2_boxp_necOnly.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[2]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 2')
dev.off()

wilcox.test(dat[dat$V1 == "yes",]$V2, dat[dat$V1 == "no",]$V2)

# 	Wilcoxon rank sum test

# data:  dat[dat$V1 == "yes", ]$V2 and dat[dat$V1 == "no", ]$V2
# W = 417, p-value = 0.0007773
# alternative hypothesis: true location shift is not equal to 0


y <- t(PF$basis[,3]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
pdf("figs/factor3_boxp_necOnly.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[3]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 3')
dev.off()

wilcox.test(dat[dat$V1 == "yes",]$V2, dat[dat$V1 == "no",]$V2)

# 	Wilcoxon rank sum test

# data:  dat[dat$V1 == "yes", ]$V2 and dat[dat$V1 == "no", ]$V2
# W = 436, p-value = 0.0001261
# alternative hypothesis: true location shift is not equal to 0

####################
#TAXONOMY BARCHART
####################
# order by similar composition 
hc <- hclust(dist(philr.t), method="complete")
df2 <- data.frame(cluster=cutree(hc,5), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
orderlist <- as.vector(dend_data$labels$label)

# any organism that makes up less than 1.5% of the total dataset collapsed into "other"
dat <- read.table("collapsed_otu_table_for_barchart.rdat.txt", header=T, sep="\t")
datmelt <- melt(dat)
datmelt$variable <- factor(datmelt$variable, levels=orderlist)
cols <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#000000")

pdf("figs/taxonomy_barchart_cluster_necOnly.pdf")
ggplot(datmelt, aes(fill=datmelt$L6, x=datmelt$variable, y=datmelt$value)) + geom_bar(stat="identity", position="fill") + theme_minimal() + theme(axis.text.x=element_text(angle=90)) + coord_flip() + scale_fill_manual(values=cols)
dev.off()

#unclustered
datmelt <- melt(dat)
pdf("figs/taxonomy_barchart_necOnly.pdf")
ggplot(datmelt, aes(fill=datmelt$L6, x=datmelt$variable, y=datmelt$value)) + geom_bar(stat="identity", position="fill") + theme_minimal() + theme(axis.text.x=element_text(angle=90)) + coord_flip() + scale_fill_manual(values=cols)
dev.off()

#####################
#CO OCCURANCE NETWORK
#####################
ps.filt <- filter_taxa(ps.dat, function(x) sum(x>1) > (0.25*length(x)), TRUE) # remove taxa not found at least 1 times in 25% of the samples to limit rare ASVs
res.cor <- cor(otu_table(ps.filt))
corrplot(res.cor)
mb.net <- spiec.easi(ps.filt, method="mb", lambda.min.ratio=1e-2, nlambda=20, icov.select.params=list(rep.num=100)) ; beep(4)
ig2.mb <- adj2igraph(mb.net$refit$stars, vertex.attr=list(name=taxa_names(ps.dat)))
plot_network(ig2.mb, ps.dat, type='taxa', color='V3')




















