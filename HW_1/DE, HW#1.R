library(GEOquery)
library(Biobase)
library(ggplot2)
library(reshape2)
library(limma)
library(MASS)
library(dplyr)
library(fgsea)
library(sva)
library(stringr) #optional
library(pathview) #optional
library(gage) #optional
library(ggrepel)

gse122121 <- getGEO("GSE122121", AnnotGPL = TRUE)[[1]]

#look at the data - already in logarithmic form
head(exprs(gse122121))

#quntile_normalization
exprs(gse122121) <- normalizeBetweenArrays(exprs(gse122121), method="quantile")

#get specific columns according to given task
gse122121 <- gse122121[, c(4:9, 13:18)]
pData(gse122121)$title

#processing the data
pData(gse122121)$rep <- gsub(".*(rep\\d)$", "\\1", pData(gse122121)$title)
pData(gse122121) <- pData(gse122121)[, c("characteristics_ch1.1", "characteristics_ch1.2", "rep")]
colnames(pData(gse122121)) <- c("genotype", "infection", "replicate")

#look now at the phenotype data
head(pData(gse122121))
colnames(fData(gse122121))

#look at the features
fData(gse122121) <- fData(gse122121)[, c("ID", "Gene symbol", "Gene ID")]
head(fData(gse122121))

#weed out non-specific (with ///), unknown(with empty fields) and duplicated genes 
# + count genes mean expression and order them in accordance to it 
gse122121 <- gse122121[!grepl("///", fData(gse122121)$`Gene symbol`), ]
gse122121 <- gse122121[fData(gse122121)$`Gene symbol` != "", ]
fData(gse122121)$mean_expression <- apply(exprs(gse122121), 1, mean)
gse122121 <- gse122121[order(fData(gse122121)$mean_expression, decreasing = TRUE), ]
gse122121 <- gse122121[!duplicated(fData(gse122121)$`Gene ID`), ]

#how many genes do we have at this time?
dim(gse122121)

#take significant 12000 genes (just a traditional number)
gse122121 <- gse122121[seq_len(12000), ]
#calculate limma (choose significant factors and so on)
gse122121.design <- model.matrix(~0+genotype, data=pData(gse122121))
colnames(gse122121.design) <- c("Ripk3KO", "WT")
fit <- lmFit(gse122121, gse122121.design)
fit2 <- contrasts.fit(fit, makeContrasts(WT - Ripk3KO, levels=gse122121.design))
fit2 <- eBayes(fit2, trend = T)
de <- topTable(fit2, adjust.method="BH", number=Inf, sort.by = "P")

#____PCA

#use gene names as row names in an expression matrix 
symbols = fData(gse122121)[,2]
expr_mat = exprs(gse122121)
rownames(expr_mat) = symbols
#make PCA plots
pcas <- prcomp(t(expr_mat), scale. = T)
plotData <- cbind(pcas$x[, 1:2], pData(gse122121))

ggplot(plotData, aes(x=PC1, y=PC2, color=genotype, shape=infection)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)

#PC1_variance (genes that are most influenced by PC1) [optional]
rotation <- pcas$rotation
PC1GenesDown <- head(rownames(rotation[order(rotation[, 1]), ]), 10)
PC1GenesUp <- tail(rownames(rotation[order(rotation[, 1]), ]), 10)
print(PC1GenesDown)
print(PC1GenesUp)

#PC2_variance (genes that are most influenced by PC2) [optional]
rotation <- pcas$rotation
PC2GenesDown <- head(rownames(rotation[order(rotation[, 2]), ]), 10)
PC2GenesUp <- tail(rownames(rotation[order(rotation[, 2]), ]), 10)
print(PC2GenesDown)
print(PC2GenesUp)

#____heatmap (It's not necessary)

#blueWhiteRed <- colorRampPalette(c("#3859A8", "#EEEEEE", "#EE2930"))(10)
#pheatmap(exprs(gse122121), scale="row", color=blueWhiteRed, border_color = NA, kmeans_k = 8,
 #        annotation_col = pData(gse122121), cluster_cols = F)

#pheatmap(expr_mat[c(PC1GenesDown, PC1GenesUp), ], 
#         scale="row", color=blueWhiteRed, border_color = NA,
#         annotation_col = pData(temp), cluster_cols = F)
#
#pheatmap(expr_mat[c(PC2GenesDown, PC2GenesUp), ], 
#         scale="row", color=blueWhiteRed, border_color = NA,
#         annotation_col = pData(temp), cluster_cols = F)

#____Batch (take batch effects into consideration and check how PCA plot changed
# after that)

batch <- pData(gse122121)$replicate
modcombat <- model.matrix(~1+genotype, data=pData(gse122121))
combat_gse122121 = ComBat(dat=expr_mat, batch=batch, mod=modcombat)

#look at PCA
pcas <- prcomp(t(combat_gse122121), scale. = T)
plotData <- cbind(pcas$x[, 1:2], pData(gse122121))
ggplot(plotData, aes(x=PC1, y=PC2, color=genotype, shape=infection)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)

#____Volcano plot

ggplot(de, aes(x=logFC, y=-log10(adj.P.Val), color=adj.P.Val < 0.05)) +
  geom_point() + theme_bw()

# adjust the number of adj.P.Val in order to get a picture with good resolution
ggplot(de, aes(x=logFC, y=-log10(adj.P.Val), color=adj.P.Val < 0.0001)) +
    geom_point() + theme_bw() + scale_color_manual(values=c("black", "red")) +
  geom_text_repel(data=de %>% dplyr::filter(adj.P.Val < 0.00001), aes(label=Gene.symbol, color=NULL))

#____GSEA, fgsea

#genes which expression significantly differs in samples
upRegulatedGenes <- de %>% dplyr::filter(adj.P.Val < 0.05 & logFC > 0) %>% pull("Gene.symbol")
length(upRegulatedGenes)
downRegulatedGenes <- de %>% dplyr::filter(adj.P.Val < 0.05 & logFC < 0) %>% pull("Gene.symbol")
length(downRegulatedGenes)

#randomGeneSet <- kegg.gs$`mmu04142 Lysosome`
#randomGeneSet <- randomGeneSet[randomGeneSet %in% rownames(de)]

#upload metabolic pathways for mouse (as HW I got data from experiment with mice) 
kg.mouse<- kegg.gsets("mouse")
kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]

# I calculate it but don't know why
stats <- de$t
names(stats) <- de$Gene.ID
fgseaResults <- fgseaMultilevel(kegg.gs, stats, minSize = 15, maxSize = 500)
head(fgseaResults, 3)

# 5 significant pathways
topPathwaysUp <- fgseaResults[ES > 0, ][head(order(pval), n=5), pathway]
topPathwaysDown <- fgseaResults[ES < 0, ][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#plot a FGSEA plot
plotGseaTable(kegg.gs[topPathways], stats, fgseaResults, gseaParam = 0.5)

#get a list of top 5 pathways identidicators [optional]
topFivePathways <- str_extract(topPathways[1:5], 'mmu[0-9]+')

#stats <- de$t
#names(stats) <- de$Gene.symbol
#nonRandomGeneSet <- kegg.gs$`mmu04060 Cytokine-cytokine receptor interaction`
#plotEnrichment(nonRandomGeneSet, stats)

#____GSEA, gage, pathview [optional, I did it just for getting illustrative 
# metabolic maps as plots ]

de.fc <- de$logFC
names(de.fc) <- de$Gene.ID
exp.fc <- de.fc
out.suffix = 'limma'

fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 &
  + !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 &
  + !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)

#pv.out.list <- sapply(path.ids2[1:5], function(pid) pathview(gene.data = exp.fc, pathway.id = pid, species = "mmu", out.suffix=out.suffix))
mydat <- sapply(topFivePathways, function(pid) pathview(gene.data = exp.fc, pathway.id = pid, species = "mmu", out.suffix=out.suffix))

