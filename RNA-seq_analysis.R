library(ggplot2)
library(dplyr)
library(extrafont)
library(gridExtra)
library(RGraphics)
library(data.table)
library(ggpubr)
library(viridis)
library(scales)
#library(ggbeeswarm)
library(gmodels)
library(DESeq2)

#BiocManager::install("DESeq2")
#BiocManager::install("apeglm")
library(apeglm)

#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")
library(rlang)
library(clusterProfiler)
library(enrichplot)
font_import()

organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

library(MetBrewer)
library(msigdbr)
library(ggh4x)
library(ggrepel)
library(ggdendro)
library(cowplot)
library(ggbeeswarm)
library(stringr)


cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbviridis <- c("#440154FF", "#31688EFF", "#35B779FF", "#E69F00")

palette <- met.brewer(name="Veronese", n=7, type="discrete")

met.brewer(name="Veronese", n=7, type="discrete")

palette2 <- c(palette[7], palette[3])
pallete4 <- c(palette[1], palette[3], palette[5],palette[7])
pallete5 <- c(palette[1], palette[2], palette[3], palette[7], palette[5])
pallete6 <- c(palette[1], palette[7], palette[2],palette[3], palette[6], palette[4])


library(showtext)

font_add(family = "Arial", regular = "C:\\WINDOWS\\Fonts\\arial.ttf")

showtext.auto()
showtext_opts(dpi = 300)

#Theme Information
#windowsFonts("Arial Black" = windowsFont("Arial Black"))

theme_science <- function (base_size = 7, base_family = "Arial") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=0.3), axis.text = element_text(colour = "black", size = base_size),
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(), plot.title = element_text(size = base_size),
          #axis.line.x = element_line(colour= "black", size=0.7),  axis.line.y = element_line(colour= "black", size=0.7),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), strip.text = element_text(size = base_size), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), legend.background=element_blank(),  
          strip.background = element_rect(colour = "white", size = 0.5), legend.key = element_blank())
}


makeStars <- function(x){
  stars <- c("***", "**", "*", "ns")
  vec <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

scaleFUN <- function(x) sprintf("%.2f", x)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbviridis <- c("#440154FF", "#31688EFF", "#35B779FF", "#E69F00")


#read ATF4 targets from https://pubmed.ncbi.nlm.nih.gov/30624206/
atf4_targets <- read.delim("elife-42940-supp1-v1.txt")

atf4_targets <- atf4_targets[1:69,2]

#read core_isgs 
Core_ISGs <- read.delim("Core_ISGs.txt")

ISGS <- as.character(Core_ISGs$ISGs)


nfkb <- read.delim("HALLMARK_TNFA_SIGNALING_VIA_NFKB.v2022.1.Hs.grp", header = FALSE)

nfkb <- as.character(nfkb$V1)

#read data
df_r <- fread("hit-counts/raw_counts.csv")

colnames(df_r) 

#rename columns to add replicate numbers to name

colnames(df_r) <- gsub("X1\\.|X2\\.|X7\\.|X8\\.|X13\\.|X14\\.|X15\\.|X16\\.", "R1_", colnames(df_r))
colnames(df_r) <- gsub("X3\\.|X4\\.|X9\\.|X10\\.|X17\\.|X18\\.|X19\\.|X20\\.", "R2_", colnames(df_r))
colnames(df_r) <- gsub("X5\\.|X6\\.|X11\\.|X12\\.|X21\\.|X22\\.|X23\\.|X24\\.", "R3_", colnames(df_r))

colnames(df_r) 

#select only the HCC1806 data
data_hcc1806 <- as.matrix(df_r[,2:7])

#make a row.names column based on ID
row.names(data_hcc1806) <- df_r$ID

#review data structure
head(data_hcc1806)

#make coldata df
coldata_hcc1806 <- DataFrame(condition = factor(c("NT","PACT", "NT", "PACT","NT","PACT")), 
                          row.names=as.character(colnames(data_hcc1806)))
#review coldata df
coldata_hcc1806

#make dds object
dds_hcc1806 <- DESeqDataSetFromMatrix(countData= data_hcc1806, 
                                        colData = coldata_hcc1806, 
                                        design = ~ condition)
#review dds object
dds_hcc1806

#remove genes with low counts
smallestGroupSize <- 3
keep <- rowSums(counts(dds_hcc1806) >= 10) >= smallestGroupSize
dds_hcc1806 <- dds_hcc1806[keep,]

#set reference level for comparisons
dds_hcc1806$condition <- relevel(dds_hcc1806$condition, ref = "NT")

#perform deseq analysis
dds_hcc1806 <- DESeq(dds_hcc1806)

#make results object
hcc1806 <- results(dds_hcc1806)

#view results comparisons
resultsNames(dds_hcc1806)

#perform logfoldchange shrinkage for low abundance RNAs
hcc1806_shrink <- lfcShrink(dds_hcc1806, coef = "condition_PACT_vs_NT", type = "apeglm")

#make MA plots
plotMA(hcc1806)
plotMA(hcc1806_shrink)

#how many genes are signficantly altered
sum(hcc1806$padj <0.01, na.rm=TRUE)
sum(hcc1806_shrink$padj <0.01, na.rm=TRUE)

#make data frames of results
hcc1806_df <- as.data.frame(hcc1806)
hcc1806_shrink_df <- as.data.frame(hcc1806_shrink)

#make ID column
hcc1806_df$ID <- row.names(hcc1806_df)
hcc1806_shrink_df$ID <- row.names(hcc1806_shrink_df)

write.csv(hcc1806_shrink_df, "hcc1806_deseq2_lfcshrink.csv")

#merge with counts file to get gene symbols and count information
hcc1806_df <- merge(hcc1806_df, df_r, by = "ID")
hcc1806_shrink_df <- merge(hcc1806_shrink_df, df_r, by = "ID")

#make column for signficance for labeling volcano plot
hcc1806_shrink_df$Sig <- ifelse(hcc1806_shrink_df$padj<0.05 & hcc1806_shrink_df$log2FoldChange < -1 | hcc1806_shrink_df$log2FoldChange > 1, "<0.05", ">0.05")

hcc1806_shrink_df$an_rank <- -log10(hcc1806_shrink_df$padj)*hcc1806_shrink_df$log2FoldChange

ggplot(hcc1806_shrink_df, aes(log2FoldChange, -log10(padj))) + geom_point(size = 0.3) +
  theme_science()  + labs(x = "Log2 Fold Change", y = "-log10 (FDR)", title = "HCC1806\nsgC1 vs sgPACT-2")  +
  annotate('rect', xmax = 1, xmin = -1, ymax = -log10(min(na.omit(hcc1806_shrink_df$padj))), 
           ymin = -log10(0.05), fill = "white", alpha = 0.8) +
  annotate('rect', xmax = 5, xmin = -5, ymax = -log10(0.05), 
           ymin = -log10(max(na.omit(hcc1806_shrink_df$padj)))-1, fill = "white", alpha = 0.8) +
  scale_x_continuous(limits = c(-2.5,5), breaks = c(-2.5,0,2.5,5)) + 
  geom_point(data=head(hcc1806_shrink_df[order(-hcc1806_shrink_df$an_rank),],7), aes(log2FoldChange, -log10(padj), color = `Gene.name`), size = 0.6) +
  geom_text_repel(data=head(hcc1806_shrink_df[order(-hcc1806_shrink_df$an_rank),],7), 
                  aes(log2FoldChange, -log10(padj), label=`Gene.name`, colour = `Gene.name`), family = 'Arial', size = 2, alpha = 1) +
  scale_colour_manual(values = palette, guide = 'none') + 
  geom_vline(xintercept = -1, linetype="dashed", colour="grey")+ geom_vline(xintercept = 1, linetype="dashed", colour="grey")  +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey")
ggsave("volcano_hcc1806.tiff", height=2.7, width=1.66, unit="in", dpi = 300)


#make ranking for GSE
hcc1806_shrink_df$rank <- -log10(hcc1806_shrink_df$padj)*sign(hcc1806_shrink_df$log2FoldChange)

#make genes with rank INF = NA
hcc1806_shrink_df <- do.call(data.frame, lapply(hcc1806_shrink_df, function(x) replace(x, is.infinite(x), NA)))

#make genelist based on rank values
gene_list_hcc1806 <- hcc1806_shrink_df$rank

# name the vector
names(gene_list_hcc1806) <- hcc1806_shrink_df$ID

# omit any NA values 
gene_list_hcc1806<-na.omit(gene_list_hcc1806)

# sort the list in decreasing order (#required for clusterProfiler)
gene_list_hcc1806 = sort(gene_list_hcc1806, decreasing = TRUE)

#perform gse analysis
gse_hcc1806 <- gseGO(geneList=gene_list_hcc1806, 
                        ont ="ALL", 
                        keyType = "ENSEMBL", 
                        minGSSize = 3, 
                        maxGSSize = 800, 
                        pvalueCutoff = 0.05, 
                        verbose = TRUE, 
                        OrgDb = get('org.Hs.eg.db'), 
                        pAdjustMethod = "fdr", 
                        eps = 0)
#load DOSE package
require(DOSE)

#make dotplot of GO terms
dotplot(gse_hcc1806, showCategory=20, split=".sign") + facet_grid(.~.sign)

#make data frame of GSE results
gse_hcc1806 <- as.data.frame(gse_hcc1806@result)

write.csv(gse_hcc1806, "gse_hcc1806_go.csv")

#rank by enrichment and adjusted p value
gse_hcc1806$rank2 <- gse_hcc1806$NES*-log10(gse_hcc1806$p.adjust)

#sort by decending rank
gse_hcc1806 <- dplyr::arrange(gse_hcc1806, desc(-rank2))

#subset for biological pathway or molecular function and make dotplots
gse_hcc1806_bp <- subset(gse_hcc1806, ONTOLOGY == "BP")

ggplot(gse_hcc1806_bp, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 8) + 
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count", title = "HCC1806, GO Biological Process") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank())
ggsave("go_bp_hcc1806.tiff", height = 9, width = 7, units = "in")

gse_hcc1806_mf <- subset(gse_hcc1806, ONTOLOGY == "MF")

ggplot(gse_hcc1806_mf, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 8) + 
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count", title = "HCC1806, GO Molecular Function") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank())
ggsave("go_mf_hcc1806.tiff", height = 5, width = 7, units = "in")

#get hallmark gene sets

H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, ensembl_gene)
head(H_t2g)

#make an ATF4 targets gene set
gene_ids <- data.frame(gs_name = "ATF4_target", ensembl_gene = df_r$ID, symbol = df_r$Gene.name)

ATF4_t2g <- subset(gene_ids, symbol %in% toupper(atf4_targets))

ATF4_t2g$symbol <- NULL

#combine hallmark and ATF4 targets
H_t2g <- rbind(H_t2g, ATF4_t2g)

#perform GSEA for hallmark and ATF4 targets
set.seed(1)
H_hcc1806 <- GSEA(gene_list_hcc1806, TERM2GENE = H_t2g, nPerm =10000)
head(H_hcc1806)

#make dotplot of GO terms
dotplot(H_hcc1806, showCategory=20, split=".sign") + facet_grid(.~.sign)


#make nice dotplot for HM 
#for HCC1806 see mb468 below

gseaplot2(H_hcc1806, geneSetID = 1, title = paste("HCC1806 ", H_hcc1806$Description[1], "\n Adjusted p-value ", round(H_hcc1806$p.adjust[1], 3), sep = ""),
          subplots = 1:2, base_size = 7) 
ggsave("gsea_NFKB_HCC1806.tiff", height = 2.5, width = 4, units = "in", dpi = 300)

gseaplot2(H_hcc1806, geneSetID = 4, title = paste("HCC1806 ", H_hcc1806$Description[4], "\n Adjusted p-value ", round(H_hcc1806$p.adjust[4], 3), sep = ""),
          subplots = 1:2, base_size = 7)  
ggsave("gsea_ATF4_HCC1806.tiff", height = 2.5, width = 4, units = "in", dpi = 300)

#make data frame of GSE results
H_hcc1806_df <- as.data.frame(H_hcc1806@result)

write.csv(H_hcc1806_df, "gse_hcc1806_hm.csv")

#rank by enrichment and adjusted p value
H_hcc1806_df$rank2 <- H_hcc1806_df$NES*-log10(H_hcc1806_df$p.adjust)

#sort by decending rank
H_hcc1806_df <- dplyr::arrange(H_hcc1806_df, desc(-rank2))

H_hcc1806_df$Description <- gsub("HALLMARK_", "", H_hcc1806_df$Description)
H_hcc1806_df$Description <- gsub("_", " ", H_hcc1806_df$Description)
H_hcc1806_df$Description <- gsub("target", "TARGET", H_hcc1806_df$Description)


ggplot(H_hcc1806_df, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 8) + 
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank())
ggsave("hallmark_hcc1806.tiff", height = 3, width = 5, units = "in")




######HERE BEGINS THE mb468 ANALYSIS
#select only the mb468 data
data_mb468 <- as.matrix(df_r[,8:13])

#make a row.names column based on ID
row.names(data_mb468) <- df_r$ID

#review data structure
head(data_mb468)

#make coldata df
coldata_mb468 <- DataFrame(condition = factor(c("NT","PACT", "NT", "PACT","NT","PACT")), 
                             row.names=as.character(colnames(data_mb468)))
#review coldata df
coldata_mb468

#make dds object
dds_mb468 <- DESeqDataSetFromMatrix(countData= data_mb468, 
                                      colData = coldata_mb468, 
                                      design = ~ condition)
#review dds object
dds_mb468

#remove genes with low counts
smallestGroupSize <- 3
keep <- rowSums(counts(dds_mb468) >= 10) >= smallestGroupSize
dds_mb468 <- dds_mb468[keep,]

#set reference level for comparisons
dds_mb468$condition <- relevel(dds_mb468$condition, ref = "NT")

#perform deseq analysis
dds_mb468 <- DESeq(dds_mb468)

#make results object
mb468 <- results(dds_mb468)

#view results comparisons
resultsNames(dds_mb468)

#perform logfoldchange shrinkage for low abundance RNAs
mb468_shrink <- lfcShrink(dds_mb468, coef = "condition_PACT_vs_NT", type = "apeglm")

#make MA plots
plotMA(mb468)
plotMA(mb468_shrink)

#how many genes are signficantly altered
sum(mb468$padj <0.01, na.rm=TRUE)
sum(mb468_shrink$padj <0.01, na.rm=TRUE)

#make data frames of results
mb468_df <- as.data.frame(mb468)
mb468_shrink_df <- as.data.frame(mb468_shrink)

#make ID column
mb468_df$ID <- row.names(mb468_df)
mb468_shrink_df$ID <- row.names(mb468_shrink_df)

write.csv(mb468_shrink_df, "mb468_deseq2_lfcshrink.csv")

#merge with counts file to get gene symbols and count information
mb468_df <- merge(mb468_df, df_r, by = "ID")
mb468_shrink_df <- merge(mb468_shrink_df, df_r, by = "ID")

#make column for signficance for labeling volcano plot
mb468_shrink_df$Sig <- ifelse(mb468_shrink_df$padj<0.05, "<0.05", ">0.05")

mb468_shrink_df$an_rank <- -log10(mb468_shrink_df$padj)*mb468_shrink_df$log2FoldChange

ggplot(mb468_shrink_df, aes(log2FoldChange, -log10(padj))) + geom_point(size = 0.3) +
  theme_science()  + labs(x = "Log2 Fold Change", y = "-log10 (FDR)", title = "MDA-MB-468\nsgNTA vs sgPACT-2")   +
  annotate('rect', xmax = 1, xmin = -1, ymax = -log10(min(na.omit(mb468_shrink_df$padj))), 
           ymin = -log10(0.05), fill = "white", alpha = 0.8) +
  annotate('rect', xmax = 5, xmin = -5, ymax = -log10(0.05), 
           ymin = -log10(max(na.omit(mb468_shrink_df$padj)))-1, fill = "white", alpha = 0.8) +
  scale_x_continuous(limits = c(-2.5,5), breaks = c(-2.5,0,2.5,5)) + 
  geom_point(data=head(mb468_shrink_df[order(-mb468_shrink_df$an_rank),],7), aes(log2FoldChange, -log10(padj), color = `Gene.name`), size = 0.6) +
  geom_text_repel(data=head(mb468_shrink_df[order(-mb468_shrink_df$an_rank),],7), 
                  aes(log2FoldChange, -log10(padj), label=`Gene.name`, colour = `Gene.name`), family = 'Arial', size = 2, alpha = 1) +
  scale_colour_manual(values = palette, guide = 'none') + 
  geom_vline(xintercept = -1, linetype="dashed", colour="grey")+ geom_vline(xintercept = 1, linetype="dashed", colour="grey")  +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey")
ggsave("volcano_mb468.tiff", height=2.7, width=1.66, unit="in", dpi = 300)


# 

#make ranking for GSE
mb468_shrink_df$rank <- -log10(mb468_shrink_df$padj)*sign(mb468_shrink_df$log2FoldChange)

#make genes with rank INF = NA
mb468_shrink_df <- do.call(data.frame, lapply(mb468_shrink_df, function(x) replace(x, is.infinite(x), NA)))


#make genelist based on rank values
gene_list_mb468 <- mb468_shrink_df$rank

# name the vector
names(gene_list_mb468) <- mb468_shrink_df$ID

# omit any NA values 
gene_list_mb468<-na.omit(gene_list_mb468)

# sort the list in decreasing order (#required for clusterProfiler)
gene_list_mb468 = sort(gene_list_mb468, decreasing = TRUE)

#perform gse analysis
gse_mb468 <- gseGO(geneList=gene_list_mb468, 
                     ont ="ALL", 
                     keyType = "ENSEMBL", 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = get('org.Hs.eg.db'), 
                     pAdjustMethod = "fdr", 
                     eps = 0)
#load DOSE package
require(DOSE)

#make dotplot of GO terms
dotplot(gse_mb468, showCategory=20, split=".sign") + facet_grid(.~.sign)

#make data frame of GSE results
gse_mb468 <- as.data.frame(gse_mb468@result)

write.csv(gse_mb468, "gse_mb468_go.csv")

#rank by enrichment and adjusted p value
gse_mb468$rank2 <- gse_mb468$NES*-log10(gse_mb468$p.adjust)

#sort by decending rank
gse_mb468 <- dplyr::arrange(gse_mb468, desc(-rank2))

#subset for biological pathway or molecular function and make dotplots
gse_mb468_bp <- subset(gse_mb468, ONTOLOGY == "BP")

ggplot(gse_mb468_bp, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 8) + 
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count", title = "MDA-MB-468, GO Biological Process") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank())
ggsave("go_bp_mb468.tiff", height = 5, width = 7, units = "in")

gse_mb468_mf <- subset(gse_mb468, ONTOLOGY == "MF")

ggplot(gse_mb468_mf, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 8) + 
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count", title = "MDA-MB-468, GO Molecular Function") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank())
ggsave("go_mf_mb468.tiff", height = 5, width = 7, units = "in")


set.seed(1)
H_mb468 <- GSEA(gene_list_mb468, TERM2GENE = H_t2g, nPerm =10000)
head(H_mb468)

#make dotplot of GO terms
dotplot(H_mb468, showCategory=20, split=".sign") + facet_grid(.~.sign)

gseaplot2(H_mb468, geneSetID = 1, title = paste("MDA-MB-468 ", H_mb468$Description[1], "\n Adjusted p-value ", round(H_mb468$p.adjust[1], 3), sep = ""),
          subplots = 1:2, base_size = 7)  
ggsave("gsea_NFKB_mb468.tiff", height = 2.5, width = 4, units = "in", dpi = 300)


gseaplot2(H_mb468, geneSetID = 2, title = paste("MDA-MB-468 ", H_mb468$Description[2], "\n Adjusted p-value ", round(H_mb468$p.adjust[2], 3), sep = ""), 
          subplots = 1:2, base_size = 7) 
ggsave("gsea_ATF4_mb468.tiff", height = 2.5, width = 4, units = "in", dpi = 300)

#make data frame of GSE results
H_mb468_df <- as.data.frame(H_mb468@result)

write.csv(H_mb468_df, "gse_mb468_hm.csv")

#rank by enrichment and adjusted p value
H_mb468_df$rank2 <- H_mb468_df$NES*-log10(H_mb468_df$p.adjust)

#sort by decending rank
H_mb468_df <- dplyr::arrange(H_mb468_df, desc(-rank2))

H_mb468_df$Description <- gsub("HALLMARK_", "", H_mb468_df$Description)
H_mb468_df$Description <- gsub("_", " ", H_mb468_df$Description)
H_mb468_df$Description <- gsub("target", "TARGET", H_mb468_df$Description)

ggplot(H_mb468_df, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 8) + 
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank())
ggsave("hallmark_mb468.tiff", height = 3, width = 5, units = "in")

#combine GSEA results for mb468 and HCC1806 and make plots

H_mb468_df$cell_line = "MDA-MB-468"
H_hcc1806_df$cell_line = "HCC1806"

Hallmark_df <- rbind(H_mb468_df, H_hcc1806_df)

ggplot(Hallmark_df, aes(NES, reorder(Description, NES), size = setSize, colour = `p.adjust`)) + 
  geom_point() + theme_science(base_size = 7) + geom_vline(xintercept = 0, colour = "black", linewidth = 0.3) +
  labs(x = "Normalized Enrichment Score", colour = "FDR", size = "Gene Count") +
  scale_colour_gradientn(colors = met.brewer(name="Veronese", type = "continuous")) + 
  theme(axis.title.y = element_blank(), 
        panel.grid.major.y = element_line(color = "grey", size = 0.3), axis.ticks.y = element_line(colour = "grey"), 
        axis.line.y = element_blank()) + 
  facet_grid2(cols = vars(cell_line)) + scale_size_continuous(range = c(0.3, 2)) 
ggsave("hallmark_mb468_hcc1806.tiff", height = 2.8, width = 4.2, units = "in", dpi = 300)

#repeat GSEA but keep all gene sets

H_mb468_all <- GSEA(gene_list_mb468, TERM2GENE = H_t2g, pvalueCutoff = 1, nPerm =10000)
H_mb468_all$Description

H_hcc1806_all <- GSEA(gene_list_hcc1806, TERM2GENE = H_t2g, pvalueCutoff = 1, nPerm =10000)
H_hcc1806_all$Description

H_hcc1806_all$p.adjust[28]

#make gsea plots for IFNA

gseaplot2(H_hcc1806_all, geneSetID = 28, title = paste("HCC1806 ", H_hcc1806_all$Description[28], "\n Adjusted p-value ", round(H_hcc1806_all$p.adjust[28], 3), sep = ""), 
          subplots = 1:2, base_size = 7) 
ggsave("gsea_IFNA_HCC1806.tiff", height = 2.5, width = 4, units = "in", dpi = 300)

gseaplot2(H_mb468_all, geneSetID = 23, title = paste("MDA-MB-468 ", H_mb468_all$Description[23], "\nAdjusted p-value ", round(H_mb468_all$p.adjust[23], 3), sep = ""), 
          subplots = 1:2, base_size = 7) 
ggsave("gsea_IFNA_mb468.tiff", height = 2.5, width = 4, units = "in", dpi = 300)



# Convert the DESeq transformed object to a data frame
vst_hcc1806 <- vst(dds_hcc1806)
vst_hcc1806 <- assay(vst_hcc1806)
vst_hcc1806 <- as.data.frame(vst_hcc1806)

vst_mb468 <- vst(dds_mb468)
vst_mb468 <- assay(vst_mb468)
vst_mb468 <- as.data.frame(vst_mb468)

#merge mb468 and 1806 data

vst_all <- merge(vst_mb468, vst_hcc1806, by ='row.names', all = TRUE)

row.names(vst_all) <- vst_all$Row.names

vst_all$Row.names <- NULL

#transform to get z-scores

vst_all <- as.data.frame(t(scale(t(vst_all))))

#get gene names for subsetting
vst_all$ID <- row.names(vst_all)

df_r_idents <- df_r[,c("ID", "Gene.name")]

vst_all <- merge(vst_all, df_r_idents, by = "ID")

row.names(vst_all) <- vst_all$id

vst_all$id <- NULL

vst_all <- na.omit(vst_all)

#get Hallmark IFNA genes and subset data

ifna <- H_t2g$ensembl_gene[H_t2g$gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE"]

vst_all$ISG <- grepl(paste(ifna, collapse = "|"), vst_all$ID)

vst_all_ISG <- subset(vst_all, ISG == TRUE)

row.names(vst_all_ISG) <- vst_all_ISG$Gene.name

vst_all_ISG$ID <- NULL
vst_all_ISG$ISG <- NULL

vst_all_ISG <- na.omit(vst_all_ISG)

#make long version of data.frame

vst_all_ISG_long <- tidyr::gather(vst_all_ISG, sample, value, 1:12)

vst_all_ISG$Gene.name <- NULL

#create dendogram for clustering of genes

vst_all_ISG_dendro <- as.dendrogram(hclust(d = dist(x = vst_all_ISG)))

# Create dendro
dendro_plot <- ggdendrogram(data = vst_all_ISG_dendro, rotate = TRUE)

# Preview the plot
print(dendro_plot)

#order of genes by clustering result

vst_all_ISG_order <- order.dendrogram(vst_all_ISG_dendro)

vst_all_ISG$Gene.name <- row.names(vst_all_ISG)

#set factor level of gene names by order of clustering

vst_all_ISG_long$Gene.name <- factor(x = vst_all_ISG_long$Gene.name,
                                              levels = vst_all_ISG$Gene.name[vst_all_ISG_order], 
                                              ordered = TRUE)


#adjust sample names for plotting
vst_all_ISG_long$sample <- gsub("\\.2", "-2", vst_all_ISG_long$sample)

vst_all_ISG_long$sample <- gsub("\\.", "_", vst_all_ISG_long$sample)

out <- strsplit(vst_all_ISG_long$sample, "_", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)


vst_all_ISG_long$replicate <- gsub("R", "", out$V1)

vst_all_ISG_long$sgRNA <- ifelse(grepl("sgPACT-2", out$V3), "sgPACT-2", "sgC1 or\nsgNTA")

vst_all_ISG_long$cl <- out$V2

vst_all_ISG_long$Cell_line <- ifelse(grepl("MB468", vst_all_ISG_long$cl), "MDA-MB-468", "HCC1806")


# Make a heatmap
p1_isg <- ggplot(vst_all_ISG_long, aes(x=sgRNA, y=Gene.name, fill=value, group = replicate)) + 
  geom_raster(position = "dodge") + 
  scale_fill_gradient2(low = palette[1], mid = 'white', high = palette[5]) + theme_science(base_size = 7) +
  geom_vline(xintercept = c(1.5,3.5), colour = "black", linetype = "dashed", linewidth = 0.4) +
  labs(y = "IFN stimulated genes", fill = "RNA\nz-score") +
  guides(fill=guide_colorbar(title.hjust=0.5, title.vjust = 0.5)) +
  theme(plot.margin = unit(c(0, 0.2, 0, 0), "cm"), axis.title.x = element_blank(), strip.text = element_blank(), 
        legend.position = "right", axis.text.y = element_text(colour = "white"), 
        axis.line.y = element_line(colour = "white"), 
        axis.ticks.y = element_line(colour = "white"), legend.key.width = unit(0.2, "cm")) +
  scale_x_discrete() + facet_wrap(~Cell_line, ncol = 2) 

p2_isg <- ggplot(vst_all_ISG_long, aes(sgRNA, value)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA, linewidth = 0.3) + geom_quasirandom(alpha = 0.2, size = 0.3) +  theme_science(base_size = 7) + 
  labs(y = "RNA z-score\n", x = "") + geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey", linewidth = 0.4) +
  theme(plot.margin = unit(c(0, 0.2, 0, 0), "cm"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank(), 
        axis.ticks.x = element_blank()) + facet_wrap(~Cell_line, ncol = 2) 


tiff("isg_heatmap_boxplot.tiff", height = 2.8, width = 2.6, units = "in", res = 300)
plot_grid(p2_isg, p1_isg, nrow = 2, rel_heights = c(2,5), align = "v", axis = "lr") 
dev.off()


####repeat as above for nfkb targets

nfkb <- H_t2g$ensembl_gene[H_t2g$gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]

vst_all$NFKB <- grepl(paste(nfkb, collapse = "|"), vst_all$ID)

vst_all_NFKB <- subset(vst_all, NFKB == TRUE)

row.names(vst_all_NFKB) <- vst_all_NFKB$Gene.name

vst_all_NFKB$ID <- NULL
vst_all_NFKB$NFKB <- NULL

vst_all_NFKB <- na.omit(vst_all_NFKB)

vst_all_NFKB_long <- tidyr::gather(vst_all_NFKB, sample, value, 1:12)

vst_all_NFKB$Gene.name <- NULL

vst_all_NFKB_dendro <- as.dendrogram(hclust(d = dist(x = vst_all_NFKB)))

# Create dendro
dendro_plot <- ggdendrogram(data = vst_all_NFKB_dendro, rotate = TRUE)

# Preview the plot
print(dendro_plot)

vst_all_NFKB_order <- order.dendrogram(vst_all_NFKB_dendro)

vst_all_NFKB$Gene.name <- row.names(vst_all_NFKB)

vst_all_NFKB_long$Gene.name <- factor(x = vst_all_NFKB_long$Gene.name,
                                     levels = vst_all_NFKB$Gene.name[vst_all_NFKB_order], 
                                     ordered = TRUE)



vst_all_NFKB_long$sample <- gsub("\\.2", "-2", vst_all_NFKB_long$sample)

vst_all_NFKB_long$sample <- gsub("\\.", "_", vst_all_NFKB_long$sample)

out <- strsplit(vst_all_NFKB_long$sample, "_", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)


vst_all_NFKB_long$replicate <- gsub("R", "", out$V1)

vst_all_NFKB_long$sgRNA <- ifelse(grepl("sgPACT-2", out$V3), "sgPACT-2", "sgC1 or\nsgNTA")

vst_all_NFKB_long$cl <- out$V2

vst_all_NFKB_long$Cell_line <- ifelse(grepl("MB468", vst_all_NFKB_long$cl), "MDA-MB-468", "HCC1806")

# Make a heatmap
p1_NFKB <- ggplot(vst_all_NFKB_long, aes(x=sgRNA, y=Gene.name, fill=value, group = replicate)) + 
  geom_raster(position = "dodge") + 
  scale_fill_gradient2(low = palette[1], mid = 'white', high = palette[5]) + theme_science(base_size = 7) +
  geom_vline(xintercept = c(1.5,3.5), colour = "black", linetype = "dashed", linewidth = 0.4) +
  labs(y = "NF-kB target genes", fill = "RNA\nz-score") +
  guides(fill=guide_colorbar(title.hjust=0.5, title.vjust = 0.5)) +
  theme(plot.margin = unit(c(0, 0.2, 0, 0), "cm"), axis.title.x = element_blank(), strip.text = element_blank(), 
        legend.position = "right", axis.text.y = element_text(colour = "white"), 
        axis.line.y = element_line(colour = "white"), 
        axis.ticks.y = element_line(colour = "white"), legend.key.width = unit(0.2, "cm")) +
  scale_x_discrete() + facet_wrap(~Cell_line, ncol = 2) 

p2_NFKB <- ggplot(vst_all_NFKB_long, aes(sgRNA, value)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA, linewidth = 0.3) + geom_quasirandom(alpha = 0.2, size = 0.3) +  theme_science(base_size = 7) + 
  labs(y = "RNA z-score\n", x = "") + geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey", linewidth = 0.4) +
  theme(plot.margin = unit(c(0, 0.2, 0, 0), "cm"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank(), 
        axis.ticks.x = element_blank()) + facet_wrap(~Cell_line, ncol = 2) 


tiff("NFKB_heatmap_boxplot.tiff", height = 2.8, width = 2.6, units = "in", res = 300)
plot_grid(p2_NFKB, p1_NFKB, nrow = 2, rel_heights = c(2,5), align = "v", axis = "lr") 
dev.off()


####repeat as above for atf4 targets

vst_all$ATF4 <- grepl(paste(ATF4_t2g$ensembl_gene, collapse = "|"), vst_all$ID)

vst_all_ATF4 <- subset(vst_all, ATF4 == TRUE)

row.names(vst_all_ATF4) <- vst_all_ATF4$Gene.name

vst_all_ATF4$ID <- NULL
vst_all_ATF4$ATF4 <- NULL

vst_all_ATF4 <- na.omit(vst_all_ATF4)

vst_all_ATF4_long <- tidyr::gather(vst_all_ATF4, sample, value, 1:12)

vst_all_ATF4$Gene.name <- NULL

vst_all_ATF4_dendro <- as.dendrogram(hclust(d = dist(x = vst_all_ATF4)))

# Create dendro
dendro_plot <- ggdendrogram(data = vst_all_ATF4_dendro, rotate = TRUE)

# Preview the plot
print(dendro_plot)

vst_all_ATF4_order <- order.dendrogram(vst_all_ATF4_dendro)

vst_all_ATF4$Gene.name <- row.names(vst_all_ATF4)

vst_all_ATF4_long$Gene.name <- factor(x = vst_all_ATF4_long$Gene.name,
                                     levels = vst_all_ATF4$Gene.name[vst_all_ATF4_order], 
                                     ordered = TRUE)



vst_all_ATF4_long$sample <- gsub("\\.2", "-2", vst_all_ATF4_long$sample)

vst_all_ATF4_long$sample <- gsub("\\.", "_", vst_all_ATF4_long$sample)

out <- strsplit(vst_all_ATF4_long$sample, "_", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)


vst_all_ATF4_long$replicate <- gsub("R", "", out$V1)

vst_all_ATF4_long$sgRNA <- ifelse(grepl("sgPACT-2", out$V3), "sgPACT-2", "sgC1 or\nsgNTA")

vst_all_ATF4_long$cl <- out$V2

vst_all_ATF4_long$Cell_line <- ifelse(grepl("MB468", vst_all_ATF4_long$cl), "MDA-MB-468", "HCC1806")

# Make a heatmap
p1_ATF4 <- ggplot(vst_all_ATF4_long, aes(x=sgRNA, y=Gene.name, fill=value, group = replicate)) + 
  geom_raster(position = "dodge") + 
  scale_fill_gradient2(low = palette[1], mid = 'white', high = palette[5]) + theme_science(base_size = 7) +
  geom_vline(xintercept = c(1.5,3.5), colour = "black", linetype = "dashed", linewidth = 0.4) +
  labs(y = "ATF4 target genes", fill = "RNA\nz-score") +
  guides(fill=guide_colorbar(title.hjust=0.5, title.vjust = 0.5)) +
  theme(plot.margin = unit(c(0, 0.2, 0, 0), "cm"), axis.title.x = element_blank(), strip.text = element_blank(), 
        legend.position = "right", axis.text.y = element_text(colour = "white"), 
        axis.line.y = element_line(colour = "white"), 
        axis.ticks.y = element_line(colour = "white"), legend.key.width = unit(0.2, "cm")) +
  scale_x_discrete() + facet_wrap(~Cell_line, ncol = 2) 

p2_ATF4 <- ggplot(vst_all_ATF4_long, aes(sgRNA, value)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA, linewidth = 0.3) + geom_quasirandom(alpha = 0.2, size = 0.3) +  theme_science(base_size = 7) + 
  labs(y = "RNA z-score\n", x = "") + geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey", linewidth = 0.4) +
  theme(plot.margin = unit(c(0, 0.2, 0, 0), "cm"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank(), 
        axis.ticks.x = element_blank()) + facet_wrap(~Cell_line, ncol = 2) 


tiff("ATF4_heatmap_boxplot.tiff", height = 2.8, width = 2.6, units = "in", res = 300)
plot_grid(p2_ATF4, p1_ATF4, nrow = 2, rel_heights = c(2,5), align = "v", axis = "lr") 
dev.off()