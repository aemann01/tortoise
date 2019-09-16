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
rownames(rawmetadata) <- rawmetadata$SampleID 
seqtab.filtered <- read.table("sequence_table.16s.filtered.txt", header=T, row.names=1)
taxa <- read.table("assigntax/rep_set_fix_tax_assignments.txt", header=F, sep="\t", row.names=1)
notinmeta <- setdiff(row.names(seqtab.filtered), rawmetadata$SampleID) #record samples absent in either metadata or OTU table
notinmeta
notinraw <- setdiff(rawmetadata$SampleID, row.names(seqtab.filtered))
notinraw
tree <- read.tree("rep_set.filt.tre") #load representative tree

################
#NORMALIZE DATA
################
#Centered log-ratio transformation
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
fviz_pca_biplot(pca, habillage=merge$SampleType, geom="point", addEllipses=T, ellipse.type="convex", pointshape=19, point.size=1, mean.point=F, select.var=list(cos2=0.90)) + labs(title="Aitchison Distance", ylab="PC1 (9.3%)", xlab="PC2 (11.3%)") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-20,25) + ylim(-20,25)
dev.off()

##############################
#PHILR DISTANCE (PHYLOGENETIC)
##############################
#create phyloseq object from "seqtab.filtered", "rawmetadata", "taxa", and "tree"
#tree was rooted in figtree at midpoint
tree <- read.tree("rep_set.root.tre")
ps.dat <- phyloseq(otu_table(seqtab.filtered, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(as.matrix(taxa[1])), tree)

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
fviz_pca_ind(pca, habillage=merge$SampleType, geom="point", addEllipses=F, pointshape=19, point.size=1, mean.point=T, repel=T) + labs(title="PhILR Distance") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-110,160) + ylim(-110,160)
dev.off()

############
#Upset plot
############
map <- as.matrix(read.table("map.txt", header=T, sep="\t", row.names=1))
merged <- merge(seqtab.filtered, map, by="row.names")
n <- ncol(seqtab.filtered) + 1
agg <- aggregate(merged[,2:n], by=list(merged$SampleType), FUN=sum)
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

#########
#Adonis 
#########
# is there a difference in microbial diversity across sample type?
metadata <- as(sample_data(ps.dat), "data.frame")
adonis(philr.dist ~ SampleType, data=metadata)
# Call:
# adonis(formula = philr.dist ~ SampleType, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# SampleType  1      7512  7511.9  22.837 0.19886  0.001 ***
# Residuals  92     30262   328.9         0.80114
# Total      93     37774                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##############
#Phylofactor
##############
OTUTable <- as.matrix(t(seqtab.filtered))
filt.list <- colnames(OTUTable)
filtmap <- rawmetadata[rawmetadata$SampleID %in% filt.list,]
filtmap <- filtmap[match(filt.list, filtmap$SampleID),]
x <- as.factor(filtmap$SampleType) # CHANGE ME to your variable of interest
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

wilcox.test(dat[dat$V1 == "Necropsy",]$V2, dat[dat$V1 == "Solar",]$V2)


# 	Wilcoxon rank sum test with continuity correction

# data:  dat[dat$V1 == "Necropsy", ]$V2 and dat[dat$V1 == "Solar", ]$V2
# W = 1772, p-value = 1.62e-10
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




###STILL TO DO###
####################
#Taxonomy barchart
####################
orderlist <- hc$labels
dat <- read.table("collapsed_simplified_taxonomy_for_plot.txt", header=T, sep="\t")
datmelt <- melt(dat)
datmelt$variable <- factor(datmelt$variable, levels=orderlist)
pdf("figs/taxonomy_barchart.pdf")
ggplot(datmelt, aes(fill=datmelt$taxonomy, x=datmelt$variable, y=datmelt$value)) + geom_bar(stat="identity", position="fill") + theme_classic() + theme(axis.text.x=element_text(angle=90))
dev.off()
