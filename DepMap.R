# Load necessary libraries

library(ggplot2)
library(dplyr)
library(extrafont)
library(gridExtra)
library(data.table)
library(ggpubr)
library(ggrepel)
library(ggbeeswarm)
library(forcats)
library(gridExtra)
library(cowplot)
library(patchwork)
library(biomaRt)
library(matrixStats)
library(purrr)
library(boot)


library(showtext)

font_add(family = "Arial", regular = "C:\\WINDOWS\\Fonts\\arial.ttf")

showtext.auto()
showtext_opts(dpi = 300)



# Define a custom theme for ggplot2
theme_science <- function (base_size = 8, base_family = "Arial") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=0.3), axis.text = element_text(colour = "black", size = base_size),
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(), plot.title = element_text(size = base_size),
          #axis.line.x = element_line(colour= "black", size=0.7),  axis.line.y = element_line(colour= "black", size=0.7),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank())
}



scaleFUN <- function(x) sprintf("%.2f", x)

makeStars <- function(x){
  stars <- c("***", "**", "*", "ns")
  vec <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

# Define a function for formatting numeric values
scaleFUN <- function(x) sprintf("%.2f", x)

# Set color palette
cbPalette <- c("#440154FF", "#31688EFF", "#35B779FF", "#E69F00", "#CC79A7")
palette2 <- c(cbPalette[1], cbPalette[2])

# Set working directory
setwd("/Users/kacottre/Desktop/DepMap")

# Read cell line annotations
#cell line annotations from Marcotte et al., 2016, https://github.com/neellab/bfg/tree/gh-pages/data/annotations
cell_line_subtypes <- read.delim("cell_line_subtypes_corrected.txt")

#read depmap cell line information
cell_line_info <- fread("DepMap_24Q2/Model.csv")

cell_line_ids <- cell_line_info[,c("ModelID", "CCLEName")]

# Read DEMETER2 scores
#read in DEMETER2 scores from DepMap portal
demeter2 <- fread("DepMap_24Q4/RNAi_(Achilles+DRIVE+Marcotte,_DEMETER2)_subsetted.csv")

str(demeter2)

demeter2 <- data.frame(demeter2)

#read in CHRONOS scores from DepMap portal, get just ADAR scores
CHRONOS <- fread("DepMap_24Q4/CRISPR_(DepMap_Public_24Q4+Score,_Chronos)_subsetted.csv")

str(CHRONOS)

CHRONOS <- data.frame(CHRONOS)

# Merge DEMETER2 and CHRONOS scores
Dependency <- merge(CHRONOS, demeter2, by = "V1", all = TRUE)



#ADAR = (103)
# Create an empty matrix to store results
CORS_ADAR_CHRONOS <- matrix(nrow=length(colnames(CHRONOS)), ncol=3)
# Iterate over columns of CHRONOS (excluding the first column)
for (i in 2:length(CHRONOS)) {
  a <- cor.test(CHRONOS$ADAR, CHRONOS[,i], method = "pearson", na.action = "na.exclude")
  # Store the results in the CORS_ADAR_CHRONOS matrix
  CORS_ADAR_CHRONOS[i,] <- c(colnames(CHRONOS)[i], a$estimate, a$p.value)
}


# Clean the data by removing rows with missing values
CORS_ADAR_CHRONOS <- as.data.frame(na.omit(CORS_ADAR_CHRONOS))

# Replace zeros with a small value (2E-20) to avoid issues in subsequent calculations
#CORS_ADAR_CHRONOS[CORS_ADAR_CHRONOS == 0]  <- 2E-20
CORS_ADAR_CHRONOS <- subset(CORS_ADAR_CHRONOS, !V1 == "ADAR")

# Convert character columns to numeric and populate with values from cor.test
CORS_ADAR_CHRONOS$Pearson <- as.numeric(as.character(CORS_ADAR_CHRONOS$V2))
CORS_ADAR_CHRONOS$p_value <- as.numeric(as.character(CORS_ADAR_CHRONOS$V3))
# Adjust p-values using the false discovery rate (FDR) method
CORS_ADAR_CHRONOS$FDR <- p.adjust(CORS_ADAR_CHRONOS$p_value, method = "fdr")
# Remove rows where the gene name is "ADAR"
CORS_ADAR_CHRONOS$V1 <- gsub("PRKRA", "PACT", CORS_ADAR_CHRONOS$V1)

ggplot(CORS_ADAR_CHRONOS, aes(Pearson, -log10(FDR))) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "ADAR1 CHRONOS Score vs All CHRONOS Scores") +
  geom_point(alpha = 0.07) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) + 
  geom_point(data=head(CORS_ADAR_CHRONOS[order(CORS_ADAR_CHRONOS$FDR),],5), aes(Pearson, -log10(FDR), color = V1)) + 
  geom_text_repel(data=head(CORS_ADAR_CHRONOS[order(CORS_ADAR_CHRONOS$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=V1, colour = V1), family = 'Arial', size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1.5), plot.title.position = "plot") + scale_color_manual(values = cbPalette, guide = "none")
ggsave("PACT/CORS_adar_CHRONOS.eps", height = 2, width = 2.7, units = "in", dpi = 600, device=cairo_ps)

#PRKRA = (8575)
#repeat above for correlation with PACT dependency
CORS_PACT_CHRONOS <- matrix(nrow=length(colnames(CHRONOS)), ncol=3) 
for (i in 2:length(CHRONOS)) {
  a <- cor.test(CHRONOS$PRKRA, CHRONOS[,i], method = "pearson", na.action = "na.exclude")
  CORS_PACT_CHRONOS[i,] <- c(colnames(CHRONOS)[i], a$estimate, a$p.value)
}

CORS_PACT_CHRONOS <- as.data.frame(na.omit(CORS_PACT_CHRONOS))

#CORS_PACT_CHRONOS[CORS_PACT_CHRONOS == 0]  <- 2E-20

CORS_PACT_CHRONOS$Pearson <- as.numeric(as.character(CORS_PACT_CHRONOS$V2))
CORS_PACT_CHRONOS$p_value <- as.numeric(as.character(CORS_PACT_CHRONOS$V3))
CORS_PACT_CHRONOS$FDR <- p.adjust(CORS_PACT_CHRONOS$p_value, method = "fdr")
CORS_PACT_CHRONOS <- subset(CORS_PACT_CHRONOS, !V1 == "PRKRA")


ggplot(CORS_PACT_CHRONOS, aes(Pearson, -log10(FDR))) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "PACT CHRONOS Score vs All CHRONOS Scores") +
  geom_point(alpha = 0.07) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_PACT_CHRONOS[order(CORS_PACT_CHRONOS$FDR),],5), aes(Pearson, -log10(FDR), color = V1)) + 
  geom_text_repel(data=head(CORS_PACT_CHRONOS[order(CORS_PACT_CHRONOS$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=V1, colour = V1), family = 'Arial', size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1.5), plot.title.position = "plot") + scale_color_manual(values = cbPalette, guide = "none")
ggsave("PACT/CORS_PACT_CHRONOS.eps", height = 2, width = 2.7, units = "in", dpi = 600, device=cairo_ps)



#ADAR = (103)
# Create an empty matrix to store results
CORS_ADAR_DEMETER2 <- matrix(nrow=length(colnames(demeter2)), ncol=3)
# Iterate over columns of demeter2 (excluding the first column)
for (i in 2:length(demeter2)) {
  a <- cor.test(demeter2$ADAR, demeter2[,i], method = "pearson", na.action = "na.exclude")
  # Store the results in the CORS_ADAR_DEMETER2 matrix
  CORS_ADAR_DEMETER2[i,] <- c(colnames(demeter2)[i], a$estimate, a$p.value)
}

# Clean the data by removing rows with missing values
CORS_ADAR_DEMETER2 <- as.data.frame(na.omit(CORS_ADAR_DEMETER2))

# Replace zeros with a small value (2E-20) to avoid issues in subsequent calculations
#CORS_ADAR_DEMETER2[CORS_ADAR_DEMETER2 == 0]  <- 2E-20
# Remove rows where the gene name is "ADAR"
CORS_ADAR_DEMETER2 <- subset(CORS_ADAR_DEMETER2, !V1 == "ADAR")

# Convert character columns to numeric and populate with values from cor.test
CORS_ADAR_DEMETER2$Pearson <- as.numeric(as.character(CORS_ADAR_DEMETER2$V2))
CORS_ADAR_DEMETER2$p_value <- as.numeric(as.character(CORS_ADAR_DEMETER2$V3))
# Adjust p-values using the false discovery rate (FDR) method
CORS_ADAR_DEMETER2$FDR <- p.adjust(CORS_ADAR_DEMETER2$p_value, method = "fdr")

CORS_ADAR_DEMETER2$V1 <- gsub("PRKRA", "PACT", CORS_ADAR_DEMETER2$V1)

#make volcano plot
ggplot(CORS_ADAR_DEMETER2, aes(Pearson, -log10(FDR))) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "ADAR1 DEMETER2 Score vs All DEMETER2 Scores") +
  geom_point(alpha = 0.07) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_ADAR_DEMETER2[order(CORS_ADAR_DEMETER2$FDR),],5), aes(Pearson, -log10(FDR), color = V1)) + 
  geom_text_repel(data=head(CORS_ADAR_DEMETER2[order(CORS_ADAR_DEMETER2$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=V1, colour = V1), family = 'Arial', size = 2) +
  theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + scale_color_manual(values = cbPalette, guide = "none")
ggsave("PACT/CORS_adar_DEMETER2.eps", height = 2, width = 2.5, units = "in", dpi = 600, device=cairo_ps)

#PRKRA = (8575)
#repeat above for correlation with PACT dependency
CORS_PACT_DEMETER2 <- matrix(nrow = length(colnames(demeter2)), ncol=3) 
for (i in 2:length(demeter2)) {
  a <- cor.test(demeter2$PRKRA, demeter2[,i], method = "pearson", na.action = "na.exclude")
  CORS_PACT_DEMETER2[i,] <- c(colnames(demeter2)[i], a$estimate, a$p.value)
}

CORS_PACT_DEMETER2 <- as.data.frame(na.omit(CORS_PACT_DEMETER2))

#CORS_PACT_DEMETER2[CORS_PACT_DEMETER2 == 0]  <- 2E-20
CORS_PACT_DEMETER2 <- subset(CORS_PACT_DEMETER2, !V1 == "PRKRA")

CORS_PACT_DEMETER2$Pearson <- as.numeric(as.character(CORS_PACT_DEMETER2$V2))
CORS_PACT_DEMETER2$p_value <- as.numeric(as.character(CORS_PACT_DEMETER2$V3))
CORS_PACT_DEMETER2$FDR <- p.adjust(CORS_PACT_DEMETER2$p_value, method = "fdr")



ggplot(CORS_PACT_DEMETER2, aes(Pearson, -log10(FDR))) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "PACT DEMETER2 Score vs All DEMETER2 Scores") +
  geom_point(alpha = 0.07) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_PACT_DEMETER2[order(CORS_PACT_DEMETER2$FDR),],5), aes(Pearson, -log10(FDR), color = V1)) + 
  geom_text_repel(data=head(CORS_PACT_DEMETER2[order(CORS_PACT_DEMETER2$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=V1, colour = V1), family = 'Arial', size = 2) +
  theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot", legend.title=element_blank()) + scale_color_manual(values = cbPalette, guide = "none")
ggsave("PACT/CORS_PACT_DEMETER2.eps", height = 2, width = 2.5, units = "in", dpi = 600, device=cairo_ps)

#make list of dsRBD containing genes
dsrbd <- c("STAU1", "STAU2", "DROSHA", "PACT", "ADARB1", "ILF3", "ADARB2", "CDKN2AIP", 
           "SON", "TARBP2", "DGCR8", "STRBP", "EIF2AK2", "MRPL44", "DICER1",
           "ADAD2", "ADAD1", "DHX9", "ADAR")
#subset dependency correlation dfs to just dsRBD genes
CORS_ADAR_CHRONOS_dsrbd <- subset(CORS_ADAR_CHRONOS, V1 %in% dsrbd)
CORS_ADAR_DEMETER2_dsrbd <- subset(CORS_ADAR_DEMETER2, V1 %in% dsrbd)

CORS_PACT_CHRONOS_dsrbd <- subset(CORS_PACT_CHRONOS, V1 %in% dsrbd)
CORS_PACT_DEMETER2_dsrbd <- subset(CORS_PACT_DEMETER2, V1 %in% dsrbd)

#make identifier column for type of dependency score
CORS_PACT_DEMETER2_dsrbd <- data.frame(CORS_PACT_DEMETER2_dsrbd[,c("Pearson", "FDR")], Gene_score = paste(CORS_PACT_DEMETER2_dsrbd$V1, "(DEMETER2)", sep = " "), Score = "DEMETER2")
CORS_PACT_CHRONOS_dsrbd <- data.frame(CORS_PACT_CHRONOS_dsrbd[,c("Pearson", "FDR")], Gene_score = paste(CORS_PACT_CHRONOS_dsrbd$V1, "(CHRONOS)", sep = " "), Score = "CHRONOS")

#combine correlation dfs for dsRBDs
CORS_PACT_dependency_dsrbd <- rbind(CORS_PACT_CHRONOS_dsrbd, CORS_PACT_DEMETER2_dsrbd)
CORS_PACT_dependency_dsrbd$Gene_score <- gsub("PRKRA", "PACT", CORS_PACT_dependency_dsrbd$Gene_score)

ggplot(CORS_PACT_dependency_dsrbd, aes(Pearson, -log10(FDR), shape = Score)) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "PACT dependency vs dsRBD genes dependency") +
  geom_point(alpha = 0.07) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_PACT_dependency_dsrbd[order(CORS_PACT_dependency_dsrbd$FDR),],5), aes(Pearson, -log10(FDR), colour = Gene_score)) + 
  geom_text_repel(data=head(CORS_PACT_dependency_dsrbd[order(CORS_PACT_dependency_dsrbd$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=Gene_score, colour = Gene_score), family = 'Arial', size = 2) +
  theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot", legend.position = "bottom", legend.title=element_blank()) + 
  scale_colour_manual(values = cbPalette) + guides(colour="none", shape = "none")
ggsave("PACT/CORS_PACT_DC_dsrbd.eps", height = 2, width = 2.5, units = "in", dpi = 600, device=cairo_ps)

#make identifier column for type of dependency score
CORS_PACT_DEMETER2 <- data.frame(CORS_PACT_DEMETER2[,c("Pearson", "FDR")], Gene_score = paste(CORS_PACT_DEMETER2$V1, "(DEMETER2)", sep = " "), Score = "DEMETER2")
CORS_PACT_CHRONOS <- data.frame(CORS_PACT_CHRONOS[,c("Pearson", "FDR")], Gene_score = paste(CORS_PACT_CHRONOS$V1, "(CHRONOS)", sep = " "), Score = "CHRONOS")

#combine correlation dfs
CORS_PACT_dependency <- rbind(CORS_PACT_CHRONOS, CORS_PACT_DEMETER2)
CORS_PACT_dependency$Gene_score <- gsub("ADAR", "ADAR1", CORS_PACT_dependency$Gene_score)

ggplot(CORS_PACT_dependency, aes(Pearson, -log10(FDR), shape = Score)) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "PACT dependency vs All genes depednency") +
  geom_point(alpha = 0.07) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_PACT_dependency[order(CORS_PACT_dependency$FDR),],5), aes(Pearson, -log10(FDR), colour = Gene_score)) + 
  geom_text_repel(data=head(CORS_PACT_dependency[order(CORS_PACT_dependency$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=Gene_score, colour = Gene_score), family = 'Arial', size = 2) +
  theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + 
  scale_colour_manual(values = cbPalette) + guides(colour="none", shape = "none")
ggsave("PACT/CORS_PACT_DC.eps", height = 2, width = 2.5, units = "in", dpi = 600, device=cairo_ps)

#make identifier column for type of dependency score
CORS_ADAR_DEMETER2 <- data.frame(CORS_ADAR_DEMETER2[,c("Pearson", "FDR")], Gene_score = paste(CORS_ADAR_DEMETER2$V1, "(DEMETER2)", sep = " "), Score = "DEMETER2")
CORS_ADAR_CHRONOS <- data.frame(CORS_ADAR_CHRONOS[,c("Pearson", "FDR")], Gene_score = paste(CORS_ADAR_CHRONOS$V1, "(CHRONOS)", sep = " "), Score = "CHRONOS")

#combine correlation dfs
CORS_ADAR_dependency <- rbind(CORS_ADAR_CHRONOS, CORS_ADAR_DEMETER2)
CORS_ADAR_dependency$Gene_score <- gsub("PRKRA", "PACT", CORS_ADAR_dependency$Gene_score)

ggplot(CORS_ADAR_dependency, aes(Pearson, -log10(FDR), shape = Score)) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "ADAR dependency vs All genes dependency") +
  geom_point(alpha = 0.07) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_ADAR_dependency[order(CORS_ADAR_dependency$FDR),],5), aes(Pearson, -log10(FDR), colour = Gene_score)) + 
  geom_text_repel(data=head(CORS_ADAR_dependency[order(CORS_ADAR_dependency$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=Gene_score, colour = Gene_score), family = 'Arial', size = 2) +
  theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot", legend.position = "bottom", legend.title=element_blank()) + 
  scale_colour_manual(values = cbPalette) + guides(colour="none", shape = "none")
ggsave("PACT/CORS_ADAR_DC.eps", height = 2, width = 2.5, units = "in", dpi = 600, device=cairo_ps)

#make identifier column for type of dependency score
CORS_ADAR_DEMETER2_dsrbd <- data.frame(CORS_ADAR_DEMETER2_dsrbd[,c("Pearson", "FDR")], Gene_score = paste(CORS_ADAR_DEMETER2_dsrbd$V1, "(DEMETER2)", sep = " "), Score = "DEMETER2")
CORS_ADAR_CHRONOS_dsrbd <- data.frame(CORS_ADAR_CHRONOS_dsrbd[,c("Pearson", "FDR")], Gene_score = paste(CORS_ADAR_CHRONOS_dsrbd$V1, "(CHRONOS)", sep = " "), Score = "CHRONOS")

#combine correlation dfs for dsRBDs
CORS_ADAR_dependency_dsrbd <- rbind(CORS_ADAR_CHRONOS_dsrbd, CORS_ADAR_DEMETER2_dsrbd)
CORS_ADAR_dependency_dsrbd$Gene_score <- gsub("PRKRA", "PACT", CORS_ADAR_dependency_dsrbd$Gene_score)

ggplot(CORS_ADAR_dependency_dsrbd, aes(Pearson, -log10(FDR), shape = Score)) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "ADAR dependency vs dsRBD genes dependency") +
  geom_point(alpha = 0.07) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_ADAR_dependency_dsrbd[order(CORS_ADAR_dependency_dsrbd$FDR),],5), aes(Pearson, -log10(FDR), colour = Gene_score)) + 
  geom_text_repel(data=head(CORS_ADAR_dependency_dsrbd[order(CORS_ADAR_dependency_dsrbd$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=Gene_score, colour = Gene_score), family = 'Arial', size = 2) +
  theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot", legend.position = "bottom", legend.title=element_blank()) + 
  scale_colour_manual(values = cbPalette) + guides(colour="none", shape = "none")
ggsave("PACT/CORS_ADAR_DC_dsrbd.eps", height = 2, width = 2.5, units = "in", dpi = 600, device=cairo_ps)


#subset for BRCA cell lines
cell_line_info_breast <- subset(cell_line_info, OncotreeLineage == "Breast")


#subset Dependency for BRCA cell lines
CHRONOS_breast <- subset(CHRONOS, V1 %in% cell_line_info_breast$ModelID)


#ADAR = (103)
# Create an empty matrix to store results
CORS_ADAR_CHRONOS_breast <- matrix(nrow = length(colnames(CHRONOS_breast)), ncol=3)
# Iterate over columns of CHRONOS_breast (excluding the first column)
for (i in 2:length(CHRONOS_breast)) {
  a <- cor.test(CHRONOS_breast$ADAR, CHRONOS_breast[,i], method = "pearson", na.action = "na.exclude")
  # Store the results in the CORS_ADAR_CHRONOS_breast matrix
  CORS_ADAR_CHRONOS_breast[i,] <- c(colnames(CHRONOS_breast)[i], a$estimate, a$p.value)
}

# Clean the data by removing rows with missing values
CORS_ADAR_CHRONOS_breast <- as.data.frame(na.omit(CORS_ADAR_CHRONOS_breast))

# Replace zeros with a small value (2E-20) to avoid issues in subsequent calculations
#CORS_ADAR_CHRONOS_breast[CORS_ADAR_CHRONOS_breast == 0]  <- 2E-20
# Remove rows where the gene name is "ADAR"
CORS_ADAR_CHRONOS_breast <- subset(CORS_ADAR_CHRONOS_breast, !V1 == "ADAR")

# Convert character columns to numeric and populate with values from cor.test
CORS_ADAR_CHRONOS_breast$Pearson <- as.numeric(as.character(CORS_ADAR_CHRONOS_breast$V2))
CORS_ADAR_CHRONOS_breast$p_value <- as.numeric(as.character(CORS_ADAR_CHRONOS_breast$V3))
# Adjust p-values using the false discovery rate (FDR) method
CORS_ADAR_CHRONOS_breast$FDR <- p.adjust(CORS_ADAR_CHRONOS_breast$p_value, method = "fdr")



#make volcano plot
ggplot(CORS_ADAR_CHRONOS_breast, aes(Pearson, -log10(FDR))) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "ADAR1 CHRONOS Score vs All CHRONOS Scores\n BRCA Cell Lines") +
  geom_point(alpha = 0.07) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_ADAR_CHRONOS_breast[order(CORS_ADAR_CHRONOS_breast$FDR),],5), aes(Pearson, -log10(FDR)), color = cbPalette[2]) + 
  geom_text_repel(data=head(CORS_ADAR_CHRONOS_breast[order(CORS_ADAR_CHRONOS_breast$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=V1), family = 'Arial', size = 2, colour = cbPalette[2]) +
  theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot")
ggsave("PACT/CORS_adar_CHRONOS_breast.eps", height = 2, width = 2.5, units = "in", dpi = 600, device=cairo_ps)


#repeat above for correlation with PACT dependency
CORS_PACT_CHRONOS_breast <- matrix(nrow = length(colnames(CHRONOS_breast)), ncol=3) 
for (i in 2:length(CHRONOS_breast)) {
  a <- cor.test(CHRONOS_breast$PRKRA, CHRONOS_breast[,i], method = "pearson", na.action = "na.exclude")
  CORS_PACT_CHRONOS_breast[i,] <- c(colnames(CHRONOS_breast)[i], a$estimate, a$p.value)
}

CORS_PACT_CHRONOS_breast <- as.data.frame(na.omit(CORS_PACT_CHRONOS_breast))

#CORS_PACT_CHRONOS_breast[CORS_PACT_CHRONOS_breast == 0]  <- 2E-20
#remove row where the gene is PRKRA
CORS_PACT_CHRONOS_breast <- subset(CORS_PACT_CHRONOS_breast, !V1 == "PRKRA")

CORS_PACT_CHRONOS_breast$Pearson <- as.numeric(as.character(CORS_PACT_CHRONOS_breast$V2))
CORS_PACT_CHRONOS_breast$p_value <- as.numeric(as.character(CORS_PACT_CHRONOS_breast$V3))
CORS_PACT_CHRONOS_breast$FDR <- p.adjust(CORS_PACT_CHRONOS_breast$p_value, method = "fdr")


ggplot(CORS_PACT_CHRONOS_breast, aes(Pearson, -log10(FDR))) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "PACT CHRONOS Score vs All CHRONOS Scores\n BRCA Cell Lines") +
  geom_point(alpha = 0.07) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_PACT_CHRONOS_breast[order(CORS_PACT_CHRONOS_breast$FDR),],5), aes(Pearson, -log10(FDR)), color = cbPalette[2]) + 
  geom_text_repel(data=head(CORS_PACT_CHRONOS_breast[order(CORS_PACT_CHRONOS_breast$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=V1), family = 'Arial', size = 2, colour = cbPalette[2]) +
  theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot")
ggsave("PACT/CORS_PACT_CHRONOS_breast.eps", height = 2, width = 2.5, units = "in", dpi = 600, device=cairo_ps)



#subset Dependency for BRCA cell lines
demeter2_breast <- subset(demeter2, V1 %in% cell_line_info_breast$ModelID)


#ADAR = (103)
# Create an empty matrix to store results
CORS_ADAR_DEMETER2_breast <- matrix(nrow=length(colnames(demeter2_breast)), ncol=3)
# Iterate over columns of demeter2_breast (excluding the first column)
for (i in 2:length(demeter2_breast)) {
  a <- cor.test(demeter2_breast$ADAR, demeter2_breast[,i], method = "pearson", na.action = "na.exclude")
  # Store the results in the CORS_ADAR_DEMETER2_breast matrix
  CORS_ADAR_DEMETER2_breast[i,] <- c(colnames(demeter2_breast)[i], a$estimate, a$p.value)
}

# Clean the data by removing rows with missing values
CORS_ADAR_DEMETER2_breast <- as.data.frame(na.omit(CORS_ADAR_DEMETER2_breast))

# Replace zeros with a small value (2E-20) to avoid issues in subsequent calculations
#CORS_ADAR_DEMETER2_breast[CORS_ADAR_DEMETER2_breast == 0]  <- 2E-20
# Remove rows where the gene name is "ADAR"
CORS_ADAR_DEMETER2_breast <- subset(CORS_ADAR_DEMETER2_breast, !V1 == "ADAR")

# Convert character columns to numeric and populate with values from cor.test
CORS_ADAR_DEMETER2_breast$Pearson <- as.numeric(as.character(CORS_ADAR_DEMETER2_breast$V2))
CORS_ADAR_DEMETER2_breast$p_value <- as.numeric(as.character(CORS_ADAR_DEMETER2_breast$V3))
# Adjust p-values using the false discovery rate (FDR) method
CORS_ADAR_DEMETER2_breast$FDR <- p.adjust(CORS_ADAR_DEMETER2_breast$p_value, method = "fdr")
# Remove rows where the gene name is "ADAR"
CORS_ADAR_DEMETER2_breast <- subset(CORS_ADAR_DEMETER2_breast, !V1 == "ADAR")


#make volcano plot
ggplot(CORS_ADAR_DEMETER2_breast, aes(Pearson, -log10(FDR))) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "ADAR1 DEMETER2 Score vs All DEMETER2 Scores\n BRCA Cell Lines") +
  geom_point(alpha = 0.07) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_ADAR_DEMETER2_breast[order(CORS_ADAR_DEMETER2_breast$FDR),],5), aes(Pearson, -log10(FDR)), color = cbPalette[2]) + 
  geom_text_repel(data=head(CORS_ADAR_DEMETER2_breast[order(CORS_ADAR_DEMETER2_breast$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=V1), family = 'Arial', size = 2, colour = cbPalette[2]) +
  theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot")
ggsave("PACT/CORS_adar_DEMETER2_breast.eps", height = 2, width = 2.5, units = "in", dpi = 600, device=cairo_ps)

#PRKRA = (8575)
#repeat above for correlation with PACT dependency
CORS_PACT_DEMETER2_breast <- matrix(nrow=length(colnames(demeter2_breast)), ncol=3) 
for (i in 2:length(demeter2_breast)) {
  a <- cor.test(demeter2_breast$PRKRA, demeter2_breast[,i], method = "pearson", na.action = "na.exclude")
  CORS_PACT_DEMETER2_breast[i,] <- c(colnames(demeter2_breast)[i], a$estimate, a$p.value)
}

CORS_PACT_DEMETER2_breast <- as.data.frame(na.omit(CORS_PACT_DEMETER2_breast))

#CORS_PACT_DEMETER2_breast[CORS_PACT_DEMETER2_breast == 0]  <- 2E-20
#remove rows where gene name is PRKRA
CORS_PACT_DEMETER2_breast <- subset(CORS_PACT_DEMETER2_breast, !V1 == "PRKRA")

CORS_PACT_DEMETER2_breast$Pearson <- as.numeric(as.character(CORS_PACT_DEMETER2_breast$V2))
CORS_PACT_DEMETER2_breast$p_value <- as.numeric(as.character(CORS_PACT_DEMETER2_breast$V3))
CORS_PACT_DEMETER2_breast$FDR <- p.adjust(CORS_PACT_DEMETER2_breast$p_value, method = "fdr")
CORS_PACT_DEMETER2_breast <- subset(CORS_PACT_DEMETER2_breast, !V1 == "PRKRA")


ggplot(CORS_PACT_DEMETER2_breast, aes(Pearson, -log10(FDR))) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "PACT DEMETER2 Score vs All DEMETER2 Scores\n BRCA Cell Lines") +
  geom_point(alpha = 0.07) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_PACT_DEMETER2_breast[order(CORS_PACT_DEMETER2_breast$FDR),],5), aes(Pearson, -log10(FDR)), color = cbPalette[2]) + 
  geom_text_repel(data=head(CORS_PACT_DEMETER2_breast[order(CORS_PACT_DEMETER2_breast$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=V1), family = 'Arial', size = 2, colour = cbPalette[2]) +
  theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot")
ggsave("PACT/CORS_PACT_DEMETER2_breast.eps", height = 2, width = 2.5, units = "in", dpi = 600, device=cairo_ps)





#merge dependency data and cell line information
CHRONOS_breast <- merge(CHRONOS_breast, cell_line_info, by.x = "V1", by.y = "ModelID")
CHRONOS_breast$LegacySubSubtype <- factor(CHRONOS_breast$LegacySubSubtype, c("ERneg_HER2neg", "ERpos_HER2neg", "ERneg_HER2pos", "ERpos_HER2pos", ""))

#make scatter plots
ggplot(subset(CHRONOS_breast, !LegacySubSubtype == ""), aes(ADAR, PRKRA, colour = LegacySubSubtype)) + geom_point(size = 0.5) + 
  geom_smooth(method = "lm", linewidth = 0.2) + theme_science() + 
  stat_cor(method = "pearson", label.x.npc = 0.17, label.y.npc = 0.4, size = 2) + 
  scale_colour_manual(breaks = c("ERneg_HER2neg", "ERpos_HER2neg", "ERneg_HER2pos", "ERpos_HER2pos"), 
                      labels = c("TNBC", "ER+", "HER2+", "HER2+/ER+"), values =cbPalette) + 
  labs(x = "ADAR1 CHRONOS Score", y = "PACT CHRONOS Score", colour = "Subtype")
ggsave("PACT/ADAR_PRKRA_breast_subtype_CHRONOS.eps", height = 2, width = 3, units = "in", dpi = 600, device=cairo_ps)

ggplot(CHRONOS_breast, aes(ADAR, PRKRA)) + geom_point(size = 0.5) + 
  geom_smooth(method = "lm", colour = cbPalette[2], linewidth = 0.2) + theme_science() + 
  stat_cor(method = "pearson", size = 3, label.x.npc = 0, label.y.npc = 0.07) + 
  labs(x = "ADAR1 CHRONOS Score", y = "PACT CHRONOS Score")
ggsave("PACT/ADAR_PRKRA_breast_CHRONOS.eps", height = 2, width = 2, units = "in", dpi = 600, device=cairo_ps)

ggplot(CHRONOS, aes(ADAR, PRKRA)) + geom_point(size = 0.5, alpha = 0.5) + 
  geom_smooth(method = "lm", colour = cbPalette[2], linewidth = 0.2) + theme_science() + 
  stat_cor(method = "pearson", size = 3, label.x.npc = 0, label.y.npc = 0.07) + 
  labs(x = "ADAR1 CHRONOS Score", y = "PACT CHRONOS Score")
ggsave("PACT/ADAR_PRKRA_CHRONOS.eps", height = 2, width = 2, units = "in", dpi = 600, device=cairo_ps)

CHRONOS_breast <- merge(CHRONOS_breast, cell_line_ids, by.x = "V1", by.y = "ModelID")

ggplot(CHRONOS_breast, aes(ADAR, PRKRA)) + geom_point(size = 0.5, alpha = 0.5) + 
  geom_smooth(method = "lm", colour = cbPalette[2], linewidth = 0.2) + theme_science() + 
  labs(x = "ADAR1 CHRONOS Score", y = "PACT CHRONOS Score") +
  geom_point(data= subset(CHRONOS_breast, V1 %in% c("ACH-000624", "ACH-000849", "ACH-000910", "ACH-000288")), 
             aes(ADAR, PRKRA,colour = CellLineName), size = 1) + 
  geom_text_repel(data= subset(CHRONOS_breast, V1 %in% c("ACH-000624", "ACH-000849", "ACH-000910", "ACH-000288")), 
                               aes(ADAR, PRKRA,label= CellLineName, colour = CellLineName), family = 'Arial', size = 2.5) +
  scale_colour_manual(values = cbPalette) + guides(colour = "none")
ggsave("PACT/ADAR_PRKRA_CHRONOS_labeled.eps", height = 2, width = 2, units = "in", dpi = 600, device=cairo_ps)




#merge dependency data and cell line information
demeter2_breast <- merge(demeter2_breast, cell_line_info, by.x = "V1", by.y = "ModelID")
demeter2_breast$LegacySubSubtype <- factor(demeter2_breast$LegacySubSubtype, c("ERneg_HER2neg", "ERpos_HER2neg", "ERneg_HER2pos", "ERpos_HER2pos", ""))


#make scatter plots
ggplot(subset(demeter2_breast, !LegacySubSubtype == ""), aes(ADAR, PRKRA, colour = LegacySubSubtype)) + geom_point(size = 0.5) + 
  geom_smooth(method = "lm", linewidth = 0.2) + theme_science() + 
  stat_cor(method = "pearson", label.x.npc = 0.17, label.y.npc = 0.2, size = 2) + 
  scale_colour_manual(breaks = c("ERneg_HER2neg", "ERpos_HER2neg", "ERneg_HER2pos", "ERpos_HER2pos"), 
                      labels = c("TNBC", "ER+", "HER2+", "HER2+/ER+"), values =cbPalette) + 
  labs(x = "ADAR1 DEMETER2 Score", y = "PACT DEMETER2 Score", colour = "Subtype")
ggsave("PACT/ADAR_PRKRA_breast_subtype_DEMETER2.eps", height = 2, width = 3, units = "in", dpi = 600, device=cairo_ps)

ggplot(demeter2_breast, aes(ADAR, PRKRA)) + geom_point(size = 0.5) + 
  geom_smooth(method = "lm", colour = cbPalette[2], linewidth = 0.2) + theme_science() + 
  stat_cor(method = "pearson", size = 3, label.x.npc = 0, label.y.npc = 0.07) + 
  labs(x = "ADAR1 DEMETER2 Score", y = "PACT DEMETER2 Score")
ggsave("PACT/ADAR_PRKRA_breast_DEMETER2.eps", height = 2, width = 2, units = "in", dpi = 600, device=cairo_ps)

ggplot(demeter2, aes(ADAR, PRKRA)) + geom_point(size = 0.5, alpha = 0.5) + 
  geom_smooth(method = "lm", colour = cbPalette[2], linewidth = 0.2) + theme_science() + 
  stat_cor(method = "pearson", size = 3, label.x.npc = 0, label.y.npc = 0.07) + 
  labs(x = "ADAR1 DEMETER2 Score", y = "PACT DEMETER2 Score")
ggsave("PACT/ADAR_PRKRA_DEMETER2.eps", height = 2, width = 2, units = "in", dpi = 600, device=cairo_ps)

#merge dependency with cell line information data frames

CHRONOS_long <- tidyr::gather(CHRONOS, "Gene", "CHRONOS", -V1)
demeter2_long <- tidyr::gather(demeter2, "Gene", "DEMETER2", -V1)

CHRONOS_long$ID <- paste(CHRONOS_long$V1, CHRONOS_long$Gene, sep = "_")
demeter2_long$ID <- paste(demeter2_long$V1, demeter2_long$Gene, sep = "_")

Dependency_all <- merge(CHRONOS_long, demeter2_long, by = "ID", all = TRUE)

#split row.names by '_' and make new columns with cell_line and Site
out <- strsplit(Dependency_all$ID, "_", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)

Dependency_all$DepMap_ID <- out$V1

Dependency_all$Gene <- out$V2

Dependency_all$ID <- NULL

Dependency_all$V1.x <- NULL
Dependency_all$V1.y <- NULL
Dependency_all$Gene.x <- NULL
Dependency_all$Gene.y <- NULL

Dependency_pact <- subset(Dependency_all, Gene == "PRKRA")

Dependency_pact$Gene <- NULL

Dependency_pact <- merge(Dependency_pact, cell_line_info, by.x = "DepMap_ID", by.y = 'ModelID')

cell_line_subtypes$cell_line <- toupper(cell_line_subtypes$cell_line)

Dependency_pact <- merge(Dependency_pact, cell_line_subtypes, by.x = 'StrippedCellLineName', by.y = 'cell_line', all = TRUE)

Dependency_pact <- subset(Dependency_pact, !OncotreePrimaryDisease == "Non-Cancerous")

Dependency_pact_nounknown <- subset(Dependency_pact,!LegacySubSubtype == "")

#make boxplots
aov_PACT_CHRONOS <- aov(CHRONOS ~ LegacySubSubtype, subset(Dependency_pact_nounknown, !is.na(CHRONOS) & OncotreeLineage == "Breast"))
summary(aov_PACT_CHRONOS)

ggplot(subset(Dependency_pact_nounknown, !is.na(CHRONOS) & OncotreeLineage == "Breast"), aes(reorder(LegacySubSubtype, CHRONOS), CHRONOS)) +
  geom_boxplot(outliers = FALSE) + geom_jitter(width = 0.2, alpha = 0.5) + theme_science() + 
  labs(y = "PACT CHRONOS Score") +
  scale_x_discrete(limits = c("ERneg_HER2neg","ERpos_HER2neg", "ERneg_HER2pos", "ERpos_HER2pos"), labels = c("ERneg_HER2neg" = "TNBC", "ERpos_HER2neg" = "ER+", "ERneg_HER2pos" = "HER2+", "ERpos_HER2pos" = "HER2+\nER+")) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_hline(yintercept = -0.5, linetype = "dashed", colour = 'grey', linewidth = 0.3) +
  theme(axis.title.x = element_blank())
ggsave("PACT/PACT_CHRONOS_subtype.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)


aov_PACT_demeter2 <- aov(DEMETER2 ~ LegacySubSubtype, subset(Dependency_pact_nounknown, !is.na(DEMETER2) & OncotreeLineage == "Breast"))
summary(aov_PACT_demeter2)

PACT_demeter2 <- TukeyHSD(aov_PACT_demeter2)

PACT_demeter2_df <- as.data.frame(PACT_demeter2$LegacySubSubtype[,1:4])

ggplot(subset(Dependency_pact_nounknown, !is.na(DEMETER2) & OncotreeLineage == "Breast"), aes(reorder(LegacySubSubtype, DEMETER2), DEMETER2)) +
  geom_boxplot(outliers = FALSE) + geom_jitter(width = 0.2, alpha = 0.5) + theme_science() + 
  labs(y = "PACT DEMETER2 Score") +
  scale_x_discrete(limits = c("ERneg_HER2neg","ERpos_HER2neg", "ERneg_HER2pos", "ERpos_HER2pos"), labels = c("ERneg_HER2neg" = "TNBC", "ERpos_HER2neg" = "ER+", "ERneg_HER2pos" = "HER2+", "ERpos_HER2pos" = "HER2+\nER+")) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_hline(yintercept = -0.5, linetype = "dashed", colour = 'grey', linewidth = 0.3) +
  theme(axis.title.x = element_blank())
ggsave("PACT/PACT_DEMETER2_subtype.eps", width = 2, height = 2, units = "in", dpi =300)

#for correlation between ADAR and PACT by lineage----
Dependency_adar <- subset(Dependency_all, Gene == "ADAR")

Dependency_adar$Gene <- NULL

Dependency_pact_adar <- merge(Dependency_pact, Dependency_adar, by = "DepMap_ID")

Dependency_pact_adar <- data.frame(PACT_chronos = Dependency_pact_adar$CHRONOS.x, ADAR1_chronos = Dependency_pact_adar$CHRONOS.y,
                                   lineage = Dependency_pact_adar$OncotreeLineage)

Dependency_pact_adar <- na.omit(Dependency_pact_adar)

func <- function(Dependency_pact_adar)
{
  return(data.frame(COR = cor(Dependency_pact_adar$PACT_chronos, Dependency_pact_adar$ADAR1_chronos)))
}
library(plyr)
cor_lineage <- ddply(Dependency_pact_adar, .(lineage), func)

Dependency_all$is_pact <- ifelse(Dependency_all$Gene == "PRKRA", "PACT", "Other Genes")


#Lineage plots ----

p1 <- ggplot(subset(Dependency_all, !is.na(CHRONOS)), aes(CHRONOS, colour = is_pact)) +
  geom_density() + theme_science(base_size = 8) + labs(y = "Density") + scale_color_manual(values = cbPalette) +
  geom_vline(xintercept = 0, linewidth = 0.3) + geom_vline(xintercept = -0.5, linetype = "dashed", colour = 'grey', linewidth = 0.3) +
  scale_x_continuous(limits = quantile(Dependency_all$CHRONOS, na.rm = T, probs = c(0.0001, 0.9999))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.35,0.8), legend.title = element_blank(), legend.text = element_text(size = 5), 
        legend.key.size = unit(0.1, "in")) 

p2 <- ggplot(subset(Dependency_pact, !is.na(CHRONOS)), aes(x = fct_reorder(OncotreeLineage, -CHRONOS, median), CHRONOS)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) + theme_science(base_size = 8) + 
  labs(y = "PACT CHRONOS Score") + scale_y_continuous(limits = quantile(Dependency_all$CHRONOS, na.rm = T, probs = c(0.0001, 0.9999))) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_hline(yintercept = -0.5, linetype = "dashed", colour = 'grey', linewidth = 0.3) +
  theme(axis.title.y = element_blank(),panel.background = element_blank()) + coord_flip()

free(p1, type = "label") + p2 + plot_layout(nrow = 2, heights = c(2,7), axes = "collect")
ggsave("PACT/PACT_CHRONOS_lineage.eps", height = 4.4, width = 3.15, units = "in", dpi = 600, device=cairo_ps)


p1 <- NULL
p2 <- NULL

p1 <- ggplot(subset(Dependency_all, !is.na(DEMETER2)), aes(DEMETER2, colour = is_pact)) +
  geom_density() + theme_science(base_size = 7) + labs(y = "Density") + scale_color_manual(values = cbPalette) +
  geom_vline(xintercept = 0, linewidth = 0.3) + geom_vline(xintercept = -0.5, linetype = "dashed", colour = 'grey', linewidth = 0.3) +
  scale_x_continuous(limits = quantile(Dependency_all$DEMETER2, na.rm = T, probs = c(0.0001, 0.9999))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.35,0.8), legend.title = element_blank(), legend.text = element_text(size = 5), 
        legend.key.size = unit(0.1, "in")) 

p2 <- ggplot(subset(Dependency_pact, !is.na(DEMETER2)), aes(x = fct_reorder(OncotreeLineage, -DEMETER2, median), DEMETER2)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) + theme_science(base_size = 7) + 
  labs(y = "PACT DEMETER2 Score") + scale_y_continuous(limits = quantile(Dependency_all$DEMETER2, na.rm = T, probs = c(0.0001, 0.9999))) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_hline(yintercept = -0.5, linetype = "dashed", colour = 'grey', linewidth = 0.3) +
  theme(axis.title.y = element_blank(),panel.background = element_blank()) + coord_flip()

free(p1, type = "label") + p2 + plot_layout(nrow = 2, heights = c(2,7), axes = "collect")
ggsave("PACT/PACT_DEMETER2_lineage.eps", height = 4.12, width = 3, units = "in", dpi = 600, device=cairo_ps)

p1 <- NULL
p2 <- NULL

Dependency_all$is_pact <- NULL

ggplot(subset(Dependency_pact, !is.na(CHRONOS) & OncotreeLineage == "Breast"), aes(reorder(subtype_neve, CHRONOS), CHRONOS)) +
  geom_boxplot(outliers = FALSE) + geom_jitter(width = 0.2, alpha = 0.5) + theme_science() + 
  labs(y = "PACT CHRONOS Score", colour = "Subtype") +
  scale_x_discrete(labels = c("luminal" = "Luminal", "her2" = "HER2+", "basala" = "Basal A", "basalb" = "Basal B")) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_hline(yintercept = -0.5, linetype = "dashed", colour = 'grey') +
  theme(axis.title.x = element_blank())
ggsave("PACT/PACT_CHRONOS_neve.eps", width = 2.2, height = 2, units= "in", dpi = 600, device=cairo_ps)


ggplot(subset(Dependency_pact, !is.na(DEMETER2) & OncotreeLineage == "Breast"), aes(reorder(subtype_neve, DEMETER2), DEMETER2)) +
  geom_boxplot(outliers = FALSE) + geom_jitter(width = 0.2, alpha = 0.5) + theme_science() + 
  labs(y = "PACT DEMETER2 Score", colour = "Subtype") +
  scale_x_discrete(labels = c("luminal" = "Luminal", "her2" = "HER2+", "basala" = "Basal A", "basalb" = "Basal B")) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_hline(yintercept = -0.5, linetype = "dashed", colour = 'grey') +
  theme(axis.title.x = element_blank())
ggsave("PACT/PACT_DEMETER2_neve.eps", width = 2.2, height = 2, units= "in", dpi = 600, device=cairo_ps)


##### determine if PACT is best at discriminating TNBC from other subtypes, 
#merge to identify breast cancer lines
Dependency_info_breast <- merge(Dependency_all, cell_line_info_breast, by.x = "DepMap_ID", by.y = "ModelID")

#classify as TNBC or not
Dependency_info_breast$TNBC <- ifelse(grepl("ERneg_HER2neg", Dependency_info_breast$LegacySubSubtype), "TNBC", "nonTNBC")

#remove non-cancerous models
Dependency_info_breast <- subset(Dependency_info_breast, !OncotreePrimaryDisease == "Non-Cancerous")

#group by gene and TNBC classification
Dependency_info_breast_sum <- dplyr::group_by(Dependency_info_breast, Gene, TNBC)

#summarize and determine average dependency scores
Dependency_info_breast_sum <- dplyr::summarise(Dependency_info_breast_sum, chronos_mean = mean(na.omit(CHRONOS)),
                                              demeter2_mean = mean(na.omit(DEMETER2)))


#spread data for each dependency score
Dependency_info_breast_sum_wide_chronos <- tidyr::spread(Dependency_info_breast_sum[,c("Gene", "TNBC", "chronos_mean")], 
                                                         TNBC, chronos_mean)
setnames(Dependency_info_breast_sum_wide_chronos, old = c('TNBC','nonTNBC'), new = c('TNBC_chronos','nonTNBC_chronos'))

Dependency_info_breast_sum_wide_demeter2 <- tidyr::spread(Dependency_info_breast_sum[,c("Gene", "TNBC", "demeter2_mean")], 
                                                          TNBC, demeter2_mean)
setnames(Dependency_info_breast_sum_wide_demeter2, old = c('TNBC','nonTNBC'), new = c('TNBC_demeter2','nonTNBC_demeter2'))

#merge all of the data.frames back together
Dependency_info_breast_sum_wide <- Reduce(function(x,y) merge(x,y,by="Gene",all=TRUE),
                                          list(Dependency_info_breast_sum_wide_chronos,
                                               Dependency_info_breast_sum_wide_demeter2))

#get difference in dependency scores
Dependency_info_breast_sum_wide$c_diff_chronos <- Dependency_info_breast_sum_wide$TNBC_chronos - Dependency_info_breast_sum_wide$nonTNBC_chronos
Dependency_info_breast_sum_wide$c_diff_demeter2 <- Dependency_info_breast_sum_wide$TNBC_demeter2 - Dependency_info_breast_sum_wide$nonTNBC_demeter2


Dependency_info_breast_sum_wide <- Dependency_info_breast_sum_wide[with(Dependency_info_breast_sum_wide,order(c_diff_chronos)),]
Dependency_info_breast_sum_wide_top_chronos <- Dependency_info_breast_sum_wide[1:5,]

Dependency_info_breast_sum_wide <- Dependency_info_breast_sum_wide[with(Dependency_info_breast_sum_wide,order(c_diff_demeter2)),]
Dependency_info_breast_sum_wide_top_demeter2 <- Dependency_info_breast_sum_wide[1:15,]


#CHRONOS
Dependency_info_breast_sum_wide_top_chronos$Gene <- factor(Dependency_info_breast_sum_wide_top_chronos$Gene, levels = Dependency_info_breast_sum_wide_top_chronos$Gene[order(Dependency_info_breast_sum_wide_top_chronos$c_diff_chronos)], ordered = TRUE)

#group by gene and TNBC classification
Dependency_info_breast_sum <- dplyr::group_by(Dependency_info_breast, Gene)

Dependency_info_breast_sum_CHRONOS <- subset(Dependency_info_breast_sum, !is.na(CHRONOS))

Dependency_info_breast_sum_CHRONOS_keep <- dplyr::summarise(Dependency_info_breast_sum_CHRONOS, tnbc = length(unique(TNBC)))

Dependency_info_breast_sum_CHRONOS_keep <- subset(Dependency_info_breast_sum_CHRONOS_keep, tnbc == 2)

Dependency_info_breast_sum_CHRONOS <- subset(Dependency_info_breast_sum_CHRONOS, Gene %in% Dependency_info_breast_sum_CHRONOS_keep$Gene)

#summarize and determine average dependency scores
Dependency_info_breast_sum_CHRONOS <- dplyr::summarise(Dependency_info_breast_sum_CHRONOS, CHRONOS_pvalue = ks.test(CHRONOS ~ TNBC)$p.value)

Dependency_info_breast_sum_CHRONOS$fdr <- p.adjust(Dependency_info_breast_sum_CHRONOS$CHRONOS_pvalue, method = "fdr")

Dependency_info_breast_sum_wide_top_chronos <- merge(Dependency_info_breast_sum_wide_top_chronos, Dependency_info_breast_sum_CHRONOS, by = "Gene")

Dependency_info_breast_CHRONOS <- subset(Dependency_info_breast, Gene %in% Dependency_info_breast_sum_wide_top_chronos$Gene)

Dependency_info_breast_sum_wide_top_chronos$Gene <- gsub("PRKRA", "PACT", Dependency_info_breast_sum_wide_top_chronos$Gene)
Dependency_info_breast_CHRONOS$Gene <- gsub("PRKRA", "PACT", Dependency_info_breast_CHRONOS$Gene)


p1 <- ggplot(Dependency_info_breast_sum_wide_top_chronos, aes(fct_rev(Gene), c_diff_chronos)) + geom_col(colour = "black", linewidth = 0.3, fill = "white") + 
  theme_science() + coord_flip() + scale_y_continuous(limits = c(-0.7, 0)) + geom_hline(linewidth = 0.3, yintercept = 0) +
  theme(axis.title.y = element_blank()) + 
  labs(y = "CHRONOS Score\nDifference\n(TNBC - nonTNBC)")


p2 <- ggplot(subset(Dependency_info_breast_CHRONOS, !is.na(CHRONOS)), aes(x = Gene, CHRONOS, colour = TNBC)) +
  geom_beeswarm(size = 0.5, dodge.width = 1) + geom_crossbar(stat = "summary", fun = "mean", position = position_dodge(width = 1)) + 
  theme_science() + scale_x_discrete(limits = sort(Dependency_info_breast_sum_wide_top_chronos$Gene, decreasing = TRUE)) +
  labs(y = "CHRONOS Score") + scale_colour_manual(values = cbPalette) + 
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_hline(yintercept = -0.5, linetype = "dashed", colour = 'grey', linewidth = 0.3) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
        axis.line.y = element_blank(), legend.position = c(0.15,0.1), axis.ticks.y = element_blank(), 
        legend.key.size = unit(0.1, "cm"), legend.title = element_blank(), legend.background = element_blank()) + coord_flip()

(p1 + theme(plot.margin = unit(c(0,5,0,0), "pt"))) + p2 + plot_layout(nrow = 1, widths = c(4,6))
ggsave("PACT/TNBC_top_chronos.eps", height = 2.1, width = 4.3, units = "in", dpi = 600, device=cairo_ps)

p1 <- NULL
p2 <- NULL




rm(CHRONOS)
rm(CHRONOS_long)
rm(demeter2)
rm(demeter2_long)


gc()

####RNA expression
RNA <- fread("DepMap_24Q4/Batch_corrected_Expression_Public_24Q4_subsetted.csv")
str(RNA)

#merge RNA expression with dependency for PACT
Dependency_RNA <- merge(Dependency_pact, RNA, by.x = "DepMap_ID", by.y = "V1", all = TRUE)

#remove non cancer cell lines
Dependency_RNA <- subset(Dependency_RNA, !OncotreePrimaryDisease == "Non-Cancerous")

#subset for BRCA only
Dependency_RNA_breast <- subset(Dependency_RNA, OncotreeLineage == "Breast" & !LegacySubSubtype == "")

#make boxplots
#make plots of interest and perform hypothesis testing
aov_PRKRA_tnbc <- aov(PRKRA ~ LegacySubSubtype, Dependency_RNA_breast)
summary(aov_PRKRA_tnbc)

ggplot(subset(Dependency_RNA_breast, !is.na(PRKRA)), aes(reorder(LegacySubSubtype, - PRKRA), PRKRA)) +
  geom_boxplot(outliers = FALSE, linewidth = 0.3) + geom_jitter(width = 0.2, size= 0.5, alpha = 0.5) + theme_science() + 
  labs(y = "PACT RNA (log2 TPM)") + theme(axis.title.x = element_blank()) +
  scale_x_discrete(limits = c("ERneg_HER2neg", "ERpos_HER2neg", "ERneg_HER2pos", "ERpos_HER2pos"), labels = c("ERneg_HER2neg" = "TNBC", "ERpos_HER2neg" = "ER+", "ERneg_HER2pos" = "HER2+", "ERpos_HER2pos" = "HER2+\nER+"))
ggsave("PACT/PACT_subtype.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)

aov_PRKRA_neve <- aov(PRKRA ~ subtype_neve, Dependency_RNA_breast)
summary(aov_PRKRA_neve)

ggplot(subset(Dependency_RNA_breast, !is.na(PRKRA) & !is.na(subtype_neve)), aes(reorder(subtype_neve, PRKRA), PRKRA)) +
  geom_boxplot(outliers = FALSE, linewidth = 0.3) + geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) + theme_science() + 
  labs(y = "PACT RNA (log2 TPM)") +
  theme(axis.title.x = element_blank()) 
ggsave("PACT/PACT_neve.eps", width = 2, height = 1.5, units = "in", dpi = 600, device=cairo_ps)


#Add row.name labels to RNA
RNA <- data.frame(RNA, row.names = RNA$V1)

RNA$V1 <- NULL

str(RNA)

#determine RNA expression z-scores
RNA_z <- as.data.frame(scale(as.matrix(RNA)))

#split col.names
out <- strsplit(colnames(RNA_z), "..", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)

colnames(RNA_z) <- out$V1

RNA_z$V1 <- row.names(RNA_z)

#merge RNA z-scores with dependency scores
Dependency_RNA_z <- merge(Dependency_pact, RNA_z, by.x = "DepMap_ID", by.y = "V1", all = TRUE)

#remove non-cancer cell lines
Dependency_RNA_z <- subset(Dependency_RNA_z, !OncotreePrimaryDisease == "Non-Cancerous")


#read core_isgs 
Core_ISGs <- read.delim("Core_ISGs.txt")

ISGS <- as.character(Core_ISGs$ISGs)

#get core_isg_score, select only core_isgs from rnaseq data
CCLE_Z_ISG_Core <- dplyr::select(Dependency_RNA_z, one_of(ISGS)) 

#get median z-score
Core_ISG_scores <- data.frame(Core_ISG_median = rowMedians(as.matrix(CCLE_Z_ISG_Core)), DepMap_ID = Dependency_RNA_z$DepMap_ID)

#merge PACT and ADAR1 dependency data frames
Dependency_pact_adar <- merge(Dependency_pact, Dependency_adar, by = "DepMap_ID")

#merge core isg scores with dependency data
Dependency_pact_adar <- merge(Dependency_pact_adar, Core_ISG_scores, by.x = "DepMap_ID", by.y = "DepMap_ID")


#plots for ISG vs PACT or ADAR1 dependency
ggplot(Dependency_pact_adar, aes(CHRONOS.x, Core_ISG_median)) + geom_point(size = 0.5, alpha = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "top") + theme_science() +
  labs(x = "PACT CHRONOS Score", y = "Core ISG Score") 
ggsave("PACT/PACT_CHRONOS_ISG.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)

ggplot(Dependency_pact_adar, aes(CHRONOS.y, Core_ISG_median)) + geom_point(size = 0.5, alpha = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "top") + theme_science() +
  labs(x = "ADAR CHRONOS Score", y = "Core ISG Score") 
ggsave("PACT/ADAR_CHRONOS_ISG.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)


ggplot(subset(Dependency_pact_adar, OncotreeLineage == "Breast"), aes(CHRONOS.x, Core_ISG_median)) + geom_point(size = 0.5, alpha = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "bottom") + theme_science() +
  labs(x = "PACT CHRONOS Score", y = "Core ISG Score") 
ggsave("PACT/PACT_CHRONOS_ISG_breast.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)

ggplot(subset(Dependency_pact_adar, OncotreeLineage == "Breast"), aes(CHRONOS.y, Core_ISG_median)) + geom_point(size = 0.5, alpha = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "bottom") + theme_science() +
  labs(x = "ADAR CHRONOS Score", y = "Core ISG Score") 
ggsave("PACT/ADAR_CHRONOS_ISG_breast.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)



ggplot(subset(Dependency_pact_adar, OncotreeLineage == "Breast"), aes(CHRONOS.x, CHRONOS.y, colour = Core_ISG_median)) + geom_point(size = 0.5, alpha = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "top") + theme_science() +
  labs(x = "PACT CHRONOS Score", y = "ADAR1 CHRONOS Score") + scale_color_gradient2(breaks = c(min(Dependency_pact_adar$Core_ISG_median),
                                                                                    high = min(Dependency_pact_adar$Core_ISG_median),
                                                                                    mid = 0))
ggsave("PACT/PACT_ADAR_ISG.eps", width = 3.5, height = 2, units = "in", dpi = 600, device=cairo_ps)


#read depmap proteomics data
proteomics <- fread("Proteomics/Harmonized_MS_CCLE_Gygi_subsetted.csv")
str(proteomics)

#split row.names by ' '
out <- strsplit(colnames(proteomics), " ", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)

out$V2 <- gsub("\\(|\\)", "", out$V2)

#make protein ids data.frame
protein_ids <- data.frame(IDs = colnames(proteomics), Symbol = out$V2)
protein_ids <- protein_ids[-1,]

#merge proteomics data with dependency data
Dependency_protein <- merge(Dependency_pact, proteomics, by.x = "DepMap_ID", by.y = "V1", all = TRUE)

#determine correlation between protein abundance (for all proteins) and PACT chronos score
CORS_PACT_CHRONOS_protein<- matrix(nrow=length(colnames(Dependency_protein)), ncol=3) 
for (i in 63:length(Dependency_protein)) {
  a <- cor.test(Dependency_protein$CHRONOS, Dependency_protein[,i], method = "pearson", na.action = "na.exclude")
  CORS_PACT_CHRONOS_protein[i,] <- c(colnames(Dependency_protein)[i], a$estimate, a$p.value)
}


CORS_PACT_CHRONOS_protein <- as.data.frame(na.omit(CORS_PACT_CHRONOS_protein))

#CORS_PACT_CHRONOS_protein[CORS_PACT_CHRONOS_protein == 0]  <- 2E-20

CORS_PACT_CHRONOS_protein$Pearson <- as.numeric(as.character(CORS_PACT_CHRONOS_protein$V2))
CORS_PACT_CHRONOS_protein$p_value <- as.numeric(as.character(CORS_PACT_CHRONOS_protein$V3))
CORS_PACT_CHRONOS_protein$FDR <- p.adjust(CORS_PACT_CHRONOS_protein$p_value, method = "fdr")
CORS_PACT_CHRONOS_protein <- merge(CORS_PACT_CHRONOS_protein, protein_ids, by.x = "V1", by.y = "IDs")

ggplot(CORS_PACT_CHRONOS_protein, aes(Pearson, -log10(FDR))) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "PACT CHRONOS Score vs All Proteins") +
  geom_point(alpha = 0.07, size = 0.5) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_PACT_CHRONOS_protein[order(CORS_PACT_CHRONOS_protein$FDR),],5), aes(Pearson, -log10(FDR), colour = Symbol), size = 0.5) + 
  geom_text_repel(data=head(CORS_PACT_CHRONOS_protein[order(CORS_PACT_CHRONOS_protein$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=Symbol, colour = Symbol), family = 'Arial', size = 2, segment.size = 0.3) +
  theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + scale_color_manual(values = cbPalette, guide = "none")
ggsave("PACT/CORS_PACT_CHRONOS_protein.eps", height = 2, width = 2.5, units = "in", dpi = 600, device=cairo_ps)


Dependency_protein <- merge(Dependency_protein, Core_ISG_scores, by.x = "DepMap_ID", by.y = "DepMap_ID")


ggplot(Dependency_protein, aes(CHRONOS, `P19525 (EIF2AK2)`, colour = Core_ISG_median)) + geom_point(size = 0.5, alpha = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "bottom") + theme_science() +
  labs(x = "PACT CHRONOS Score", y = "PKR Protein (z-score)") + scale_colour_gradient2(mid = 0)
ggsave("PACT/PACT_CHRONOS_PKR_protein_ISG.eps", width = 3, height = 2, units = "in", dpi = 600, device=cairo_ps)


ggplot(subset(Dependency_protein, OncotreeLineage == "Breast"), aes(CHRONOS, `P19525 (EIF2AK2)`, colour = Core_ISG_median)) + 
  geom_point(size = 0.5, alpha = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "bottom") + theme_science() +
  labs(x = "PACT CHRONOS Score", y = "PKR Protein (z-score)") + scale_colour_gradient2(mid = 0)
ggsave("PACT/PACT_CHRONOS_PKR_protein_ISG_breast.eps", width = 3, height = 2, units = "in", dpi = 600, device=cairo_ps)

########Protein correlation by lineage

protein_lineage <- Dependency_protein[,c("CHRONOS", "OncotreeLineage", "P19525 (EIF2AK2)")]

protein_lineage <- na.omit(protein_lineage)

protein_lineage <- protein_lineage[protein_lineage$OncotreeLineage %in% names(which(table(protein_lineage$OncotreeLineage) > 3)), ]



protein_lineage_dfl <- split(protein_lineage[, c(1,3)], protein_lineage$OncotreeLineage)
results.lst <- lapply(protein_lineage_dfl, function(x) cor.test(x[, 1], x[, 2], method="pearson"))
results.stats <- lapply(results.lst, "[", c("estimate", "conf.int", "p.value"))
stats <- as.data.frame(do.call(rbind, lapply(results.stats, unlist)))
stats$lineage <- row.names(stats)

stats$FDR <- p.adjust(stats$p.value, method = "fdr")

library(scales)
ggplot(stats, aes(reorder(lineage, -estimate.cor), estimate.cor, ymin = conf.int1, ymax = conf.int2, colour = p.value)) +
  geom_pointrange(size = 0.3, linewidth = 0.3) + theme_science() + coord_flip() + scale_colour_gradientn(
    colours = c(cbPalette[5], cbPalette[2], cbPalette[1]),
    values = c(0, rescale(0.05, from = range(stats$p.value)), 1), breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = "", y = "Pearson Correlation, \nPACT CHRONOS vs PKR protein expression", colour = "p-value") +
  theme(legend.key.width = unit(0.25, "cm"))
ggsave("PACT/PACT_chronos_PKR_protein_lineage.eps", height = 2.2, width = 3.2, units = "in", dpi = 600, device=cairo_ps)



#repeat as above for DEMTER2 scores
Dependency_protein_DEMETER2 <- subset(Dependency_protein, !is.na(Dependency_protein$DEMETER2))

CORS_PACT_DEMETER2_protein<- matrix(nrow=length(colnames(Dependency_protein_DEMETER2)), ncol=3) 
for (i in 63:length(Dependency_protein_DEMETER2)) {
  a <- cor.test(Dependency_protein_DEMETER2$DEMETER2, Dependency_protein_DEMETER2[,i], method = "pearson", na.action = "na.exclude")
  CORS_PACT_DEMETER2_protein[i,] <- c(colnames(Dependency_protein_DEMETER2)[i], a$estimate, a$p.value)
}

CORS_PACT_DEMETER2_protein <- as.data.frame(na.omit(CORS_PACT_DEMETER2_protein))

#CORS_PACT_DEMETER2_protein[CORS_PACT_DEMETER2_protein == 0]  <- 2E-20

CORS_PACT_DEMETER2_protein$Pearson <- as.numeric(as.character(CORS_PACT_DEMETER2_protein$V2))
CORS_PACT_DEMETER2_protein$p_value <- as.numeric(as.character(CORS_PACT_DEMETER2_protein$V3))
CORS_PACT_DEMETER2_protein$FDR <- p.adjust(CORS_PACT_DEMETER2_protein$p_value, method = "fdr")
CORS_PACT_DEMETER2_protein <- merge(CORS_PACT_DEMETER2_protein, protein_ids, by.x = "V1", by.y = "IDs")

ggplot(CORS_PACT_DEMETER2_protein, aes(Pearson, -log10(FDR))) + ylab("-log10 FDR") + labs(x = "Pearson Correlation", title = "PACT DEMETER2 Score vs All Proteins") +
  geom_point(alpha = 0.07, size = 0.5) + geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", colour = "grey40") + 
  theme_science() + scale_x_continuous(limits = c(-1,1)) +
  geom_point(data=head(CORS_PACT_DEMETER2_protein[order(CORS_PACT_DEMETER2_protein$FDR),],5), aes(Pearson, -log10(FDR), colour = Symbol), size = 0.5) + 
  geom_text_repel(data=head(CORS_PACT_DEMETER2_protein[order(CORS_PACT_DEMETER2_protein$FDR),],5), 
                  aes(Pearson, -log10(FDR),label=Symbol, colour = Symbol), family = 'Arial', size = 2, segment.size = 0.3) +
  theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + scale_color_manual(values = cbPalette, guide = "none")
ggsave("PACT/CORS_PACT_DEMETER2_protein.eps", height = 2, width = 2.5, units = "in", dpi = 600, device=cairo_ps)



#make scatter plots for comparisons between PACT depednency score and PACT or PKR expression
ggplot(Dependency_protein, aes(CHRONOS, `P19525 (EIF2AK2)`)) + geom_point(size = 0.5, alpha = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "bottom") + theme_science() +
  labs(x = "PACT CHRONOS Score", y = "PKR Protein (z-score)") 
ggsave("PACT/PACT_CHRONOS_PKR_protein.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)

ggplot(subset(Dependency_protein, OncotreeLineage == "Breast"), aes(CHRONOS, `P19525 (EIF2AK2)`)) + geom_point(size = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "bottom") + theme_science() +
  labs(x = "PACT CHRONOS Score", y = "PKR Protein (z-score)")
ggsave("PACT/PACT_CHRONOS_PKR_protein_breast.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)

ggplot(Dependency_protein, aes(CHRONOS, `O75569 (PRKRA)`)) + geom_point(size = 0.5, alpha = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3) + theme_science() +
  labs(x = "PACT CHRONOS Score", y = "PACT Protein (z-score)")
ggsave("PACT/PACT_CHRONOS_PACT_protein.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)

ggplot(subset(Dependency_protein, OncotreeLineage == "Breast"), aes(CHRONOS, `O75569 (PRKRA)`)) + geom_point(size = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3) + theme_science() +
  labs(x = "PACT CHRONOS Score", y = "PACT Protein (z-score)")
ggsave("PACT/PACT_CHRONOS_PACT_protein_breast.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)




ggplot(subset(Dependency_protein, OncotreeLineage == "Breast"), aes(CHRONOS, `P19525 (EIF2AK2)`)) + geom_point(size = 0.5, alpha = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = 0) + theme_science() +
  geom_point(data= subset(Dependency_protein, DepMap_ID %in% c("ACH-000624", "ACH-000849", "ACH-000910", "ACH-000288")), 
             aes(CHRONOS, `P19525 (EIF2AK2)`,colour = CellLineName), size = 1) +
  geom_text_repel(data= subset(Dependency_protein, DepMap_ID %in% c("ACH-000624", "ACH-000849", "ACH-000910", "ACH-000288")), 
                  aes(CHRONOS, `P19525 (EIF2AK2)`,label= CellLineName, colour = CellLineName), family = 'Arial', size = 2.5) +
  labs(x = "PACT CHRONOS Score", y = "PKR Protein (z-score)") +
  scale_colour_manual(values = cbPalette) + guides(colour = "none")
ggsave("PACT/PACT_CHRONOS_PKR_protein_breast_labeled.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)

ggplot(subset(Dependency_protein, OncotreeLineage == "Breast"), aes(CHRONOS, `O75569 (PRKRA)`)) + geom_point(size = 0.5, alpha = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "bottom") + theme_science() +
  geom_point(data= subset(Dependency_protein, DepMap_ID %in% c("ACH-000624", "ACH-000849", "ACH-000910", "ACH-000288")), 
             aes(CHRONOS, `O75569 (PRKRA)`,colour = CellLineName), size = 1) +
  geom_text_repel(data= subset(Dependency_protein, DepMap_ID %in% c("ACH-000624", "ACH-000849", "ACH-000910", "ACH-000288")), 
                  aes(CHRONOS, `O75569 (PRKRA)`,label= CellLineName, colour = CellLineName), family = 'Arial', size = 2) +
  labs(x = "PACT CHRONOS Score", y = "PACT Protein (z-score)") +
  scale_colour_manual(values = cbPalette) + guides(colour = "none")
ggsave("PACT/PACT_CHRONOS_PACT_protein_breast_labeled.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)



ggplot(Dependency_protein, aes(DEMETER2, `P19525 (EIF2AK2)`)) + geom_point(size = 0.5, alpha = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "bottom") + theme_science() +
  labs(x = "PACT DEMETER2 Score", y = "PKR Protein (z-score)") 
ggsave("PACT/PACT_DEMETER2_PKR_protein.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)

ggplot(subset(Dependency_protein, OncotreeLineage == "Breast"), aes(DEMETER2, `P19525 (EIF2AK2)`)) + geom_point(size = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "bottom") + theme_science() +
  labs(x = "PACT DEMETER2 Score", y = "PKR Protein (z-score)")
ggsave("PACT/PACT_DEMETER2_PKR_protein_breast.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)

ggplot(Dependency_protein, aes(DEMETER2, `O75569 (PRKRA)`)) + geom_point(size = 0.5, alpha = 0.5) + 
  geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "top") + theme_science() +
  labs(x = "PACT DEMETER2 Score", y = "PACT Protein (z-score)")
ggsave("PACT/PACT_DEMETER2_PACT_protein.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)

ggplot(subset(Dependency_protein, OncotreeLineage == "Breast"), aes(DEMETER2, `O75569 (PRKRA)`)) + 
  geom_point(size = 0.5) + geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "top") + theme_science() +
  labs(x = "PACT DEMETER2 Score", y = "PACT Protein (z-score)")
ggsave("PACT/PACT_DEMETER2_PACT_protein_breast.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)


#scatter plots for ATF4, ATF4 and PACT

ggplot(Dependency_protein, aes(`P18848 (ATF4)`, `O75569 (PRKRA)`)) + geom_point(size = 0.5) + 
  geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "top") + theme_science() +
  labs(x = "ATF4 Protein (z-score)", y = "PACT Protein (z-score)")
ggsave("PACT/PACT_ATF4_protein.eps", width = 1.7, height = 1.7, units = "in", dpi = 600, device=cairo_ps)


ggplot(Dependency_protein, aes(`P18848 (ATF4)`, `P18847 (ATF3)`)) + geom_point(size = 0.5) + 
  geom_smooth(method = "lm", linewidth = 0.3, colour = cbPalette[2]) + 
  stat_cor(size = 3, label.x.npc = "left", label.y.npc = "top") + theme_science() +
  labs(x = "ATF4 Protein (z-score)", y = "ATF3 Protein (z-score)")
ggsave("PACT/ATF3_ATF4_protein.eps", width = 1.7, height = 1.7, units = "in", dpi = 600, device=cairo_ps)



Dependency_protein_breast <- subset(Dependency_protein, OncotreeLineage == "Breast")


#make boxplots
aov_PRKRA_tnbc <- aov(`O75569 (PRKRA)` ~ LegacySubSubtype, Dependency_protein_breast)
summary(aov_PRKRA_tnbc)

ggplot(subset(Dependency_protein_breast, !is.na(`O75569 (PRKRA)`)), aes(reorder(LegacySubSubtype, - `O75569 (PRKRA)`), `O75569 (PRKRA)`)) +
  geom_boxplot(outliers = FALSE, linewidth = 0.3) + geom_jitter(width = 0.2, size= 0.5, alpha = 0.5) + theme_science() + 
  labs(y = "PACT protein (z-score)") + theme(axis.title.x = element_blank()) +
  scale_x_discrete(limits = c("ERneg_HER2neg", "ERpos_HER2neg", "ERneg_HER2pos", "ERpos_HER2pos"), labels = c("ERneg_HER2neg" = "TNBC", "ERpos_HER2neg" = "ER+", "ERneg_HER2pos" = "HER2+", "ERpos_HER2pos" = "HER2+\nER+"))
ggsave("PACT/PACT_subtype_protein.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)


aov_PRKRA_neve <- aov(`O75569 (PRKRA)` ~ subtype_neve, Dependency_protein_breast)
summary(aov_PRKRA_neve)
ggplot(subset(Dependency_protein_breast, !is.na(`O75569 (PRKRA)`) & !is.na(subtype_neve)), aes(reorder(subtype_neve, - `O75569 (PRKRA)`), `O75569 (PRKRA)`)) +
  geom_boxplot(outliers = FALSE, linewidth = 0.3) + geom_jitter(width = 0.2, size= 0.5, alpha = 0.5) + theme_science() + 
  labs(y = "PACT protein (z-score)") + theme(axis.title.x = element_blank()) 
ggsave("PACT/PACT_neve_protein.eps", width = 2, height = 1.5, units = "in", dpi = 600, device=cairo_ps)


aov_PKR_tnbc <- aov(`P19525 (EIF2AK2)` ~ LegacySubSubtype, Dependency_protein_breast)
summary(aov_PKR_tnbc)

ggplot(subset(Dependency_protein_breast, !is.na(`P19525 (EIF2AK2)`)), aes(reorder(LegacySubSubtype, - `P19525 (EIF2AK2)`), `P19525 (EIF2AK2)`)) +
  geom_boxplot(outliers = FALSE, linewidth = 0.3) + geom_jitter(width = 0.2, size= 0.5, alpha = 0.5) + theme_science() + 
  labs(y = "PKR protein (z-score)") + theme(axis.title.x = element_blank()) +
  scale_x_discrete(limits = c("ERneg_HER2neg", "ERpos_HER2neg", "ERneg_HER2pos", "ERpos_HER2pos"), labels = c("ERneg_HER2neg" = "TNBC", "ERpos_HER2neg" = "ER+", "ERneg_HER2pos" = "HER2+", "ERpos_HER2pos" = "HER2+\nER+"))
ggsave("PACT/PKR_subtype_protein.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)


aov_PKR_neve <- aov(`P19525 (EIF2AK2)` ~ subtype_neve, Dependency_protein_breast)
summary(aov_PKR_neve)


ggplot(subset(Dependency_protein_breast, !is.na(`P19525 (EIF2AK2)`) & !is.na(subtype_neve)), aes(reorder(subtype_neve, - `P19525 (EIF2AK2)`), `P19525 (EIF2AK2)`)) +
  geom_boxplot(outliers = FALSE, linewidth = 0.3) + geom_jitter(width = 0.2, size= 0.5, alpha = 0.5) + theme_science() + 
  labs(y = "PKR protein (z-score)") + theme(axis.title.x = element_blank()) 
ggsave("PACT/PKR_neve_protein.eps", width = 2, height = 1.5, units = "in", dpi = 600, device=cairo_ps)






