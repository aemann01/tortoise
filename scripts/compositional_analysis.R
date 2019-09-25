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
library(UpSetR)
library(reshape2)
library(vegan)
library(phylofactor)
library(ggtree)
library(DESeq2)
library(igraph)

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
rawmetadata <- read_delim(file = file.path(PATH, "map.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries
rownames(rawmetadata) <- rawmetadata$sample_id 
seqtab.filtered <- read.table("sequence_table.16s.filtered.txt", header=T, row.names=1)
taxa <- read.table("assigntax/rep_set_fix_tax_assignments_6levels.txt", header=F, sep="\t", row.names=1)
notinmeta <- setdiff(row.names(seqtab.filtered), rawmetadata$sample_id) #record samples absent in either metadata or OTU table
notinmeta
notinraw <- setdiff(rawmetadata$sample_id, row.names(seqtab.filtered))
notinraw
tree <- read.tree("rep_set.filt.tre") #load representative tree

################
#NORMALIZE DATA
################
#Centered log-ratio transformation -- for Aitchison only
seqtab.clr <- clr(seqtab.filtered)
write.table(as.data.frame(seqtab.clr), "sequence_table.16s.clr.txt", quote=F, sep="\t", col.names=NA)

#############################
#AITCHISON DISTANCE (LINEAR)
#############################
# Biplot
#have to read back in clr sequence table
comp <- read.table("sequence_table.16s.clr.txt", header=T, row.names=1)
pca <- prcomp(as.matrix(comp))
pdf("figs/aitchison_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("figs/aitchison_pca.pdf")
fviz_pca_biplot(pca, habillage=rawmetadata$tortoise_group, geom="point", addEllipses=T, ellipse.type="convex", pointshape=19, point.size=1, mean.point=F, select.var=list(cos2=0.90)) + labs(title="Aitchison Distance", ylab="PC1 (9.3%)", xlab="PC2 (11.3%)") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-20,25) + ylim(-20,25)
dev.off()

##############################
#PHILR DISTANCE (PHYLOGENETIC)
##############################
#create phyloseq object from "seqtab.filtered", "rawmetadata", "taxa", and "tree"
#tree was rooted in figtree at midpoint
tree <- read.tree("rep_set.root.tre")
ps.dat <- phyloseq(otu_table(seqtab.filtered, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(as.matrix(taxa)), tree)

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
pdf("figs/philr_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("figs/philr_pca.pdf")
fviz_pca_ind(pca, habillage=rawmetadata$tortoise_group, geom="point", addEllipses=F, pointshape=19, point.size=1, mean.point=T, repel=T) + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-110,160) + ylim(-110,160)
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
pdf("figs/wunifrac_dendrogram.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()

#pca plot of rarefied data
pca <- prcomp(wuni)
pdf("figs/wunifrac_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("figs/wunifrac_pca_group.pdf")
fviz_pca_ind(pca, habillage=merge$tortoise_group, geom="point", pointsize=3, mean.point=FALSE) + labs(title="Weighted Unifrac Distance") + theme_classic()
dev.off()

#############################
#DESEQ DIFFERENTIAL ABUNDACE
#############################
# positive log fold == enriched in solar, negative log cold == enriched in necropsy
ps.deseq <- phyloseq_to_deseq2(ps.dat, ~tortoise_group)
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
pdf("figs/logfold_change.pdf")
ggplot(sigtab, aes(x=V7, y=log2FoldChange, color=V3)) + geom_point(size=3) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + coord_flip() + theme_minimal()
dev.off()

############
#UPSET PLOT
############
map <- as.matrix(read.table("map.txt", header=T, sep="\t", row.names=1))
merged <- merge(seqtab.filtered, map, by="row.names")
n <- ncol(seqtab.filtered) + 1
agg <- aggregate(merged[,2:n], by=list(merged$tortoise_group), FUN=sum)
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
#convert to presence absence table -- ignore warnining message, still works
agg[agg>1] <- 1
#transpose again
agg <- setNames(data.frame(t(agg[,-1])), agg[,1])
#upsetR 
pdf("figs/upset.pdf", onefile=F)
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Total ASVs per sample type", mb.ratio = c(0.55, 0.45))
dev.off()

#########################
#ALPHA DIVERSITY BY GROUP
#########################
pdf("figs/adiv_group.pdf")
plot_richness(ps.dat, x="tortoise_group", measures=c("Observed", "Shannon")) + theme_minimal() + geom_boxplot(alpha=0.0)
dev.off()
pdf("figs/adiv_age.pdf")
plot_richness(ps.dat, x="agegroup", color="tortoise_group", measures=c("Observed", "Shannon")) + theme_minimal() + geom_boxplot(alpha=0.0)
dev.off()
adiv <- estimate_richness(ps.dat)
# get p value for solar vs necropsy
wilcox.test(adiv[grepl("N", rownames(adiv)),]$Observed, adiv[grepl("S", rownames(adiv)),]$Observed)

# 	Wilcoxon rank sum test with continuity correction

# data:  adiv[grepl("N", rownames(adiv)), ]$Observed and adiv[grepl("S", rownames(adiv)), ]$Observed
# W = 261.5, p-value = 2.215e-08
# alternative hypothesis: true location shift is not equal to 0
wilcox.test(adiv[grepl("N", rownames(adiv)),]$Shannon, adiv[grepl("S", rownames(adiv)),]$Shannon)

# 	Wilcoxon rank sum test with continuity correction

# data:  adiv[grepl("N", rownames(adiv)), ]$Shannon and adiv[grepl("S", rownames(adiv)), ]$Shannon
# W = 373, p-value = 3.159e-06
# alternative hypothesis: true location shift is not equal to 0

##########
#PERMANOVA
##########
# is there a difference in microbial diversity across sample type?
metadata <- as(sample_data(ps.dat), "data.frame")
adonis(philr.dist ~ tortoise_group, data=metadata)

# Call:
# adonis(formula = philr.dist ~ tortoise_group, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#                Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
# tortoise_group  1      7423  7423.0  21.931 0.1977  0.001 ***
# Residuals      89     30124   338.5         0.8023
# Total          90     37547                 1.0000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ agegroup, data=metadata)

# Call:
# adonis(formula = philr.dist ~ agegroup, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# agegroup   2      2987 1493.27  3.8022 0.07954  0.002 **
# Residuals 88     34561  392.74         0.92046
# Total     90     37547                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ nasal_discharge, data=metadata)

# Call:
# adonis(formula = philr.dist ~ nasal_discharge, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#                 Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# nasal_discharge  2      3534 1767.22  4.5723 0.09413  0.002 **
# Residuals       88     34013  386.51         0.90587
# Total           90     37547                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ agegroup_group, data=metadata)

# Call:
# adonis(formula = philr.dist ~ agegroup_group, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# agegroup_group  4      9178 2294.52  6.9558 0.24444  0.001 ***
# Residuals      86     28369  329.87         0.75556
# Total          90     37547                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ bcs, data=metadata)

# Call:
# adonis(formula = philr.dist ~ bcs, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# bcs        1      5141  5141.2   14.12 0.13693  0.001 ***
# Residuals 89     32406   364.1         0.86307
# Total     90     37547                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


adonis(philr.dist ~ disposition, data=metadata)

# Call:
# adonis(formula = philr.dist ~ disposition, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# disposition  3      7931 2643.51  7.7654 0.21121  0.001 ***
# Residuals   87     29617  340.42         0.78879
# Total       90     37547                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ histo_rhinitis, data=metadata)

# Call:
# adonis(formula = philr.dist ~ histo_rhinitis, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# histo_rhinitis  2      8451  4225.3  12.779 0.22507  0.001 ***
# Residuals      88     29097   330.6         0.77493
# Total          90     37547                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ histo_rhinitis_chronic, data=metadata)

# Call:
# adonis(formula = philr.dist ~ histo_rhinitis_chronic, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#                        Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# histo_rhinitis_chronic  2      7310  3654.9  10.637 0.19468  0.001 ***
# Residuals              88     30237   343.6         0.80532
# Total                  90     37547                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ nasal_discharge_character, data=metadata)

# Call:
# adonis(formula = philr.dist ~ nasal_discharge_character, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#                           Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
# nasal_discharge_character  4      5497 1374.20  3.6873 0.1464  0.001 ***
# Residuals                 86     32050  372.68         0.8536
# Total                     90     37547                 1.0000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ nasal_erosion, data=metadata)

# Call:
# adonis(formula = philr.dist ~ nasal_erosion, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# nasal_erosion  1      2300 2299.87  5.8072 0.06125  0.008 **
# Residuals     89     35247  396.04         0.93875
# Total         90     37547                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##############
#PHYLOFACTOR
##############
OTUTable <- as.matrix(t(seqtab.filtered))
filt.list <- colnames(OTUTable)
filtmap <- rawmetadata[rawmetadata$sample_id %in% filt.list,]
filtmap <- filtmap[match(filt.list, filtmap$sample_id),]
x <- as.factor(filtmap$tortoise_group) # CHANGE ME to your variable of interest
tree <- phy_tree(philr.dat)
tax <- read.table("assigntax/rep_set_fix_tax_assignments_phylofactor.txt", sep="\t", header=T)
common.otus <- which(rowSums(OTUTable>0)>10)
OTUTable <- OTUTable[common.otus,]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(OTUTable)))
PF <- PhyloFactor(OTUTable, tree, x, nfactors=3)
PF$Data <- PF$Data[PF$tree$tip.label,]
gtree <- pf.tree(PF,layout="rectangular")
pdf("figs/phylofactor_tree_rectangular.pdf")
gtree$ggplot + geom_tiplab(size=2, align=T)
dev.off()
gtree <- pf.tree(PF)
pdf("figs/phylofactor_tree_circular.pdf")
gtree$ggplot
dev.off()

y <- t(PF$basis[,1]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
pdf("figs/factor1_boxp.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[1]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 1')
dev.off()
wilcox.test(dat[dat$V1 == "necropsy",]$V2, dat[dat$V1 == "solar",]$V2)

# 	Wilcoxon rank sum test with continuity correction

# data:  dat[dat$V1 == "necropsy", ]$V2 and dat[dat$V1 == "solar", ]$V2
# W = 1705, p-value = 8.828e-11
# alternative hypothesis: true location shift is not equal to 0

y <- t(PF$basis[,2]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
pdf("figs/factor2_boxp.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[2]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 2')
dev.off()

wilcox.test(dat[dat$V1 == "necropsy",]$V2, dat[dat$V1 == "solar",]$V2)

# 	Wilcoxon rank sum test with continuity correction

# data:  dat[dat$V1 == "necropsy", ]$V2 and dat[dat$V1 == "solar", ]$V2
# W = 186, p-value = 4.78e-10
# alternative hypothesis: true location shift is not equal to 0


y <- t(PF$basis[,2]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
pdf("figs/factor2_boxp.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[2]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 2')
dev.off()

wilcox.test(dat[dat$V1 == "Necropsy",]$V2, dat[dat$V1 == "Solar",]$V2)


# 	Wilcoxon rank sum test with continuity correction

# data:  dat[dat$V1 == "Necropsy", ]$V2 and dat[dat$V1 == "Solar", ]$V2
# W = 193, p-value = 3.033e-10
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
pdf("figs/taxonomy_barchart_cluster.pdf")
ggplot(datmelt, aes(fill=datmelt$L6, x=datmelt$variable, y=datmelt$value)) + geom_bar(stat="identity", position="fill") + theme_classic() + theme(axis.text.x=element_text(angle=90)) + coord_flip() + scale_fill_manual(values=cols)
dev.off()

#unclustered
datmelt <- melt(dat)
pdf("figs/taxonomy_barchart.pdf")
ggplot(datmelt, aes(fill=datmelt$L6, x=datmelt$variable, y=datmelt$value)) + geom_bar(stat="identity", position="fill") + theme_classic() + theme(axis.text.x=element_text(angle=90)) + coord_flip() + scale_fill_manual(values=cols)
dev.off()

