library(RTCGA)
library(Biobase)
library(affy)
library(genefilter)
library(XML)
library(reshape2)
library(survival)
library(survminer)
require("survival")
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(TCGAretriever)
#library(extrafont)
library(DESeq2)
library(tidyverse)
library(edgeR)
library(DescTools)
library(MetBrewer)
library(ggpubr)
library(ggrepel)

# Load fonts
#extrafont::loadfonts(device="win")
#fonts()

#library(systemfonts)

library(showtext)

font_add(family = "Arial", regular = "C:\\WINDOWS\\Fonts\\arial.ttf")

showtext.auto()
showtext_opts(dpi = 600, device=cairo_ps)




cbPalette <- c("#440154FF", "#31688EFF", "#35B779FF", "#E69F00")
palette2 <- c(cbPalette[1], cbPalette[2])



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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# read in RNAseq data for breast cancer samples from TCGA - 
RNAseq_BRCA <- fread("BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt")

#first row has relevant column names, make col.names match the first row
colnames(RNAseq_BRCA) <- paste(colnames(RNAseq_BRCA), RNAseq_BRCA[1,], sep= "-")

#hybridization reference has gene name and id as a string, split it and make a data.frame of the split
out <- strsplit(as.character(RNAseq_BRCA$`Hybridization REF`), "\\|")
test <- as.data.frame(do.call(rbind, out))


#get just the columns with raw counts
RNAseq_BRCA <- dplyr::select(RNAseq_BRCA, contains("raw_count"))

#Make a new column for gene name using the split from above
RNAseq_BRCA$Gene  <- test$V1

#remove the first row
RNAseq_BRCA <- RNAseq_BRCA[-1,]

#get rid of genes with ? or duplicate names
RNAseq_BRCA <- subset(RNAseq_BRCA, Gene != "?")
RNAseq_BRCA <- subset(RNAseq_BRCA, Gene != "SLC35E2")

#make a new data.frame with RNAseq_BRCA data with row.names = gene name
RNAseq_BRCA <- data.frame(RNAseq_BRCA, row.names = RNAseq_BRCA$Gene)

#get rid of the the gene name column
RNAseq_BRCA$Gene <- NULL

#make a new data.frame for RNAseq_BRCA where all of the raw counts are numeric
RNAseq_BRCA_df <- as.data.frame(row.names = row.names(RNAseq_BRCA), lapply(RNAseq_BRCA, function(x) as.numeric(as.character(x))))


str(RNAseq_BRCA_df)

row.names(RNAseq_BRCA_df)

#the next few lines are based off of this example: https://www.biostars.org/p/153013/
#quick check to see how many normal and primary samples there are. 
#In TCGA barcode Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29.
#the 14th character will be 0 for tumor and 1 for normal
table(substr(colnames(RNAseq_BRCA_df),14,14))

#index tumor (t) or normal (n) based on barcode value
n_index <- which(substr(colnames(RNAseq_BRCA_df),14,14) == '1')
t_index <- which(substr(colnames(RNAseq_BRCA_df),14,14) == '0')

#get CPM values for RNAseq data
logCPM <- cpm(RNAseq_BRCA_df, log=TRUE)

#define scal function, calculates a modified z-score, see link above for more details
scal <- function(x,y){
  mean_n <- rowMeans(y)  # mean of normal
  sd_n <- apply(y,1,sd)  # SD of normal
  # z score as (value - mean normal)/SD normal
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}

makeStars <- function(x){
  stars <- c("***", "**", "*", "ns")
  vec <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

#apply scal function to the logCPM data 
z_rna <- scal(logCPM,logCPM[,n_index])

#transpose z_rna and make new RNAseq_BRCA data.frame
RNAseq_BRCA <- as.data.frame(t(z_rna))


#make a column for sample barcode based on row.names
RNAseq_BRCA$Sample <- row.names(RNAseq_BRCA)

#split the sample barcode
out <- strsplit(as.character(RNAseq_BRCA$Sample), "\\.")
test <- as.data.frame(do.call(rbind, out))


#make the patient barcode using the first three values of the split
test$bcr_patient_barcode  <- do.call(paste, c(test[c("V1", "V2", "V3")], sep = "-")) 

#add the patient_barcode to the data.frame
RNAseq_BRCA$patient_barcode <- test$bcr_patient_barcode

#make the sample barcode using the first three values of the split
test$Sample_Barcode <- do.call(paste, c(test[c("V1", "V2", "V3","V4")], sep = "-")) 

#add the sample_barcode to the data.frame
RNAseq_BRCA$Sample_Barcode <- test$Sample_Barcode

#add the sample type
RNAseq_BRCA$Sample_Type <- test$V4

#rename the sample types, add DROP for the B samples for each normal, primary or metastatic
RNAseq_BRCA$Sample_Type <- gsub("11A", "Normal", RNAseq_BRCA$Sample_Type)
RNAseq_BRCA$Sample_Type <- gsub("11B", "DROP", RNAseq_BRCA$Sample_Type)

RNAseq_BRCA$Sample_Type <- gsub("01A", "Primary", RNAseq_BRCA$Sample_Type)
RNAseq_BRCA$Sample_Type <- gsub("01B", "DROP", RNAseq_BRCA$Sample_Type)

RNAseq_BRCA$Sample_Type <- gsub("06A", "Metastatic", RNAseq_BRCA$Sample_Type)
RNAseq_BRCA$Sample_Type <- gsub("06B", "DROP", RNAseq_BRCA$Sample_Type)

#drop the b samples
RNAseq_BRCA <- subset(RNAseq_BRCA, !Sample_Type == "DROP")

#subset only the primary tumor samples
RNAseq_BRCA_primary <- subset(RNAseq_BRCA, Sample_Type == "Primary")

RNAseq_BRCA_primary$Sample_Barcode

#load annotation file from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4911051/ 
Annotation <- fread("journal.pone.0157368.s008.txt")

#keep only columns 2-4 of annotation file
Annotation <- Annotation[,c(2:4,30,36,41)]

#merge with RNAseq data, patient barcode works here because there are no normal or metastic samples
RNAseq_BRCA_primary <- merge(RNAseq_BRCA_primary, Annotation, by.y = "BARCODE", by.x = "patient_barcode", all.x = TRUE)

#subset only the normal samples
RNAseq_BRCA_normal <- subset(RNAseq_BRCA, Sample_Type == "Normal")

#defne TNBC and PAM50 as Normal or Normal Breast
RNAseq_BRCA_normal$TNBC <- "Normal"
RNAseq_BRCA_normal$PAM50 <- "Normal Breast"
RNAseq_BRCA_normal$ER_CALL <- "Normal"
RNAseq_BRCA_normal$PGR_CALL <- "Normal"
RNAseq_BRCA_normal$HER2_CALL <- "Normal"

#bind primary and normal data
RNAseq_BRCA_primary_normal <- rbind(RNAseq_BRCA_normal, RNAseq_BRCA_primary)

#make plots of interest and perform hypothesis testing
aov_PRKRA <- aov(PRKRA ~ TNBC, RNAseq_BRCA_primary_normal)

summary(aov_PRKRA)

RNAseq_BRCA_primary_normal$TNBC <- factor(as.factor(RNAseq_BRCA_primary_normal$TNBC), levels = c("Normal", "NO", "YES"))

DT_tnbc <- DunnettTest(RNAseq_BRCA_primary_normal$PRKRA, RNAseq_BRCA_primary_normal$TNBC)

DT_tnbc <- as.data.frame(DT_tnbc$Normal[,1:4])

out <- strsplit(as.character(row.names(DT_tnbc)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

DT_tnbc$group1 <- out$V1
DT_tnbc$group2 <- out$V2
DT_tnbc$labels <- makeStars(DT_tnbc$pval)
DT_tnbc$y.position <- max(RNAseq_BRCA_primary_normal$PRKRA*1.05)
DT_tnbc <- subset(DT_tnbc, pval < 0.05)


ggplot(subset(RNAseq_BRCA_primary_normal, !is.na(TNBC)), aes(as.factor(TNBC), PRKRA)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "NO", "YES"), labels=c("NO" = "Non TNBC", "YES" = "TNBC")) + ylab("PACT RNA (z-Score)") +
  geom_bracket(xmin = DT_tnbc$group2, 
               xmax = DT_tnbc$group1, 
               y.position = DT_tnbc$y.position, label = DT_tnbc$labels, 
               tip.length = 0.05, step.increase = 0.1) +
  geom_boxplot(notch = TRUE, outliers = FALSE) + geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) + theme_science() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +   geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey")
ggsave("TNBC_v_PRKRA.eps", width = 2.2, height = 2, units = "in", dpi = 600, device=cairo_ps)


#make plots of interest and perform hypothesis testing
aov_EIF2AK2 <- aov(EIF2AK2 ~ TNBC, RNAseq_BRCA_primary_normal)

summary(aov_EIF2AK2)

RNAseq_BRCA_primary_normal$TNBC <- factor(as.factor(RNAseq_BRCA_primary_normal$TNBC), levels = c("Normal", "NO", "YES"))

DT_tnbc <- DunnettTest(RNAseq_BRCA_primary_normal$EIF2AK2, RNAseq_BRCA_primary_normal$TNBC)

DT_tnbc <- as.data.frame(DT_tnbc$Normal[,1:4])

out <- strsplit(as.character(row.names(DT_tnbc)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

DT_tnbc$group1 <- out$V1
DT_tnbc$group2 <- out$V2
DT_tnbc$labels <- makeStars(DT_tnbc$pval)
DT_tnbc$y.position <- max(RNAseq_BRCA_primary_normal$EIF2AK2*1.05)
DT_tnbc <- subset(DT_tnbc, pval < 0.05)

ggplot(subset(RNAseq_BRCA_primary_normal, !is.na(TNBC)), aes(as.factor(TNBC), EIF2AK2)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "NO", "YES"), labels=c("NO" = "Non TNBC", "YES" = "TNBC")) + ylab("PKR RNA (z-Score)") +
  geom_bracket(xmin = DT_tnbc$group2, 
               xmax = DT_tnbc$group1, 
               y.position = DT_tnbc$y.position, label = DT_tnbc$labels, 
               tip.length = 0.05, step.increase = 0.1) +
  geom_boxplot(notch = TRUE, outliers = FALSE) + geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) + theme_science() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +   geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey40")
ggsave("TNBC_v_EIF2AK2.eps", width = 2.2, height = 2, units = "in", dpi = 600, device=cairo_ps)




#make plots of interest and perform hypothesis testing
aov_PRKRA <- aov(PRKRA ~ PAM50, RNAseq_BRCA_primary_normal)

summary(aov_PRKRA)

RNAseq_BRCA_primary_normal$PAM50 <- factor(as.factor(RNAseq_BRCA_primary_normal$PAM50), 
                                           levels = c("Normal Breast", "Normal", "Basal", "Her2", "LumA", "LumB"))

DTpam50 <- DunnettTest(RNAseq_BRCA_primary_normal$PRKRA, RNAseq_BRCA_primary_normal$PAM50)



DT_pam50 <- DunnettTest(RNAseq_BRCA_primary_normal$PRKRA, RNAseq_BRCA_primary_normal$PAM50)

DT_pam50 <- as.data.frame(DT_pam50$Normal[,1:4])

out <- strsplit(as.character(row.names(DT_pam50)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

DT_pam50$group1 <- out$V1
DT_pam50$group2 <- out$V2
DT_pam50$labels <- makeStars(DT_pam50$pval)
DT_pam50$y.position <- max(RNAseq_BRCA_primary_normal$PRKRA*1.05)
DT_pam50 <- subset(DT_pam50, pval < 0.05)

ggplot(subset(RNAseq_BRCA_primary_normal, !is.na(PAM50)), aes(as.factor(PAM50), PRKRA)) + xlab("") + 
  scale_x_discrete(limits= c("Normal Breast", "Normal", "Basal", "Her2", "LumA", "LumB"), labels=c("Normal\nBreast", "PAM50\nNormal", "Basal", "Her2", "LumA", "LumB")) +
  ylab("PACT RNA (z-Score)") +
  geom_bracket(xmin = DT_pam50$group2, 
               xmax = DT_pam50$group1, 
               y.position = DT_pam50$y.position, label = DT_pam50$labels, 
               tip.length = 0.05, step.increase = 0.1) +
  geom_boxplot(notch = TRUE, outliers = FALSE) + geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) + theme_science() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +   geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey40")
ggsave("pam50_v_PRKRA.eps", width = 3, height = 2, units = "in", dpi = 600, device=cairo_ps)




#make plots of interest and perform hypothesis testing
aov_EIF2AK2 <- aov(EIF2AK2 ~ PAM50, RNAseq_BRCA_primary_normal)

summary(aov_EIF2AK2)

RNAseq_BRCA_primary_normal$PAM50 <- factor(as.factor(RNAseq_BRCA_primary_normal$PAM50), 
                                           levels = c("Normal Breast", "Normal", "Basal", "Her2", "LumA", "LumB"))

DTpam50 <- DunnettTest(RNAseq_BRCA_primary_normal$EIF2AK2, RNAseq_BRCA_primary_normal$PAM50)



DT_pam50 <- DunnettTest(RNAseq_BRCA_primary_normal$EIF2AK2, RNAseq_BRCA_primary_normal$PAM50)

DT_pam50 <- as.data.frame(DT_pam50$Normal[,1:4])

out <- strsplit(as.character(row.names(DT_pam50)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

DT_pam50$group1 <- out$V1
DT_pam50$group2 <- out$V2
DT_pam50$labels <- makeStars(DT_pam50$pval)
DT_pam50$y.position <- max(RNAseq_BRCA_primary_normal$EIF2AK2*1.05)
DT_pam50 <- subset(DT_pam50, pval < 0.05)


ggplot(subset(RNAseq_BRCA_primary_normal, !is.na(PAM50)), aes(as.factor(PAM50), EIF2AK2)) + xlab("") + 
  scale_x_discrete(limits= c("Normal Breast", "Normal", "Basal", "Her2", "LumA", "LumB"), labels=c("Normal\nBreast", "PAM50\nNormal", "Basal", "Her2", "LumA", "LumB")) +
  ylab("PKR RNA (z-Score)") +
  geom_bracket(xmin = DT_pam50$group2, 
               xmax = DT_pam50$group1, 
               y.position = DT_pam50$y.position, label = DT_pam50$labels, 
               tip.length = 0.05, step.increase = 0.1) +
  geom_boxplot(notch = TRUE, outliers = FALSE) + geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) + theme_science() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +   geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey40")
ggsave("pam50_v_EIF2AK2.eps", width = 3, height = 2, units = "in", dpi = 600, device=cairo_ps)





#make plots of interest and perform hypothesis testing
aov_PRKRA <- aov(PRKRA ~ Sample_Type, RNAseq_BRCA_primary_normal)

summary(aov_PRKRA)

RNAseq_BRCA$Sample_Type <- factor(as.factor(RNAseq_BRCA$Sample_Type), levels = c("Normal", "Primary", "Metastatic"))

DT_sample <- DunnettTest(RNAseq_BRCA$PRKRA, RNAseq_BRCA$Sample_Type)

DT_sample <- as.data.frame(DT_sample$Normal)

out <- strsplit(as.character(row.names(DT_sample)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

DT_sample$group1 <- out$V1
DT_sample$group2 <- out$V2
DT_sample$labels <- makeStars(DT_sample$pval)
DT_sample$y.position <- max(RNAseq_BRCA_primary_normal$PRKRA*1.05)
DT_sample <- subset(DT_sample, pval < 0.05)


ggplot(subset(RNAseq_BRCA, !is.na(Sample_Type)), aes(as.factor(Sample_Type), PRKRA)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "Primary", "Metastatic")) + ylab("PACT RNA (z-Score)") +
  geom_boxplot(notch = TRUE, outliers = FALSE) + geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) + theme_science() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +   geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey40")
ggsave("sample_v_PRKRA.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)



aov_EIF2AK2 <- aov(EIF2AK2 ~ Sample_Type, RNAseq_BRCA_primary_normal)

summary(aov_EIF2AK2)

RNAseq_BRCA$Sample_Type <- factor(as.factor(RNAseq_BRCA$Sample_Type), levels = c("Normal", "Primary", "Metastatic"))

DT_sample <- DunnettTest(RNAseq_BRCA$EIF2AK2, RNAseq_BRCA$Sample_Type)

DT_sample <- as.data.frame(DT_sample$Normal)

out <- strsplit(as.character(row.names(DT_sample)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

DT_sample$group1 <- out$V1
DT_sample$group2 <- out$V2
DT_sample$labels <- makeStars(DT_sample$pval)
DT_sample$y.position <- max(RNAseq_BRCA_primary_normal$EIF2AK2*1.05)
DT_sample <- subset(DT_sample, pval < 0.05)


ggplot(subset(RNAseq_BRCA, !is.na(Sample_Type)), aes(as.factor(Sample_Type), EIF2AK2)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "Primary", "Metastatic")) + ylab("PKR RNA (z-Score)") +
  geom_bracket(xmin = DT_sample$group2, 
               xmax = DT_sample$group1, 
               y.position = DT_sample$y.position, label = DT_sample$labels, 
               tip.length = 0.05, step.increase = 0.1) +
  geom_boxplot(notch = TRUE, outliers = FALSE) + geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) + theme_science() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +   geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey40")
ggsave("sample_v_EIF2AK2.eps", width = 2, height = 2, units = "in", dpi = 600, device=cairo_ps)



#change value for call to be subtype names


RNAseq_BRCA_primary_normal$call <- apply(RNAseq_BRCA_primary_normal[,c("ER_CALL","PGR_CALL","HER2_CALL")], 1, paste, collapse="-")

RNAseq_BRCA_primary_normal$call <- gsub("0-0-0", "TNBC", RNAseq_BRCA_primary_normal$call)

RNAseq_BRCA_primary_normal$call <- gsub("1-0-0", "ER+", RNAseq_BRCA_primary_normal$call)

RNAseq_BRCA_primary_normal$call <- gsub("0-1-0", "PR+", RNAseq_BRCA_primary_normal$call)

RNAseq_BRCA_primary_normal$call <- gsub("0-0-1", "HER2+", RNAseq_BRCA_primary_normal$call)

RNAseq_BRCA_primary_normal$call <- gsub("1-1-0", "ER+/PR+", RNAseq_BRCA_primary_normal$call)

RNAseq_BRCA_primary_normal$call <- gsub("1-0-1", "ER+/HER2+", RNAseq_BRCA_primary_normal$call)

RNAseq_BRCA_primary_normal$call <- gsub("0-1-1", "PR+/HER2+", RNAseq_BRCA_primary_normal$call)

RNAseq_BRCA_primary_normal$call <- gsub("1-1-1", "ER+/PR+/HER2+", RNAseq_BRCA_primary_normal$call)

RNAseq_BRCA_primary_normal$call <- gsub("Normal-Normal-Normal", "Normal", RNAseq_BRCA_primary_normal$call)

RNAseq_BRCA_primary_normal_type <- subset(RNAseq_BRCA_primary_normal, call %in% c("TNBC", "ER+", "PR+", "HER2+", "ER+/PR+", "ER+/HER2+", "PR+/HER2+", "ER+/PR+/HER2+", "Normal"))


aov_PRKRA <- aov(PRKRA ~ call, RNAseq_BRCA_primary_normal_type)

summary(aov_PRKRA)

RNAseq_BRCA_primary_normal_type$call <- factor(as.factor(RNAseq_BRCA_primary_normal_type$call), levels = c("Normal", "TNBC", "ER+", "PR+", "HER2+", "ER+/PR+", "ER+/HER2+", "PR+/HER2+", "ER+/PR+/HER2+"))

DT_type <- DunnettTest(RNAseq_BRCA_primary_normal_type$PRKRA, RNAseq_BRCA_primary_normal_type$call)

DT_type <- as.data.frame(DT_type$Normal)

out <- strsplit(as.character(row.names(DT_type)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

DT_type$group1 <- out$V1
DT_type$group2 <- out$V2
DT_type$labels <- makeStars(DT_type$pval)
DT_type$y.position <- max(RNAseq_BRCA_primary_normal_type$PRKRA*1.05)
DT_type <- subset(DT_type, pval < 0.05)


ggplot(subset(RNAseq_BRCA_primary_normal_type, !is.na(call)), aes(as.factor(call), PRKRA)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "TNBC", "ER+", "PR+", "HER2+", "ER+/PR+", "ER+/HER2+", "PR+/HER2+", "ER+/PR+/HER2+"),
                   labels= c("Normal", "TNBC", "ER+", "PR+", "HER2+", "ER+\nPR+", "ER+\nHER2+", "PR+\nHER2+", "ER+/PR+\nHER2+")) +
  ylab("PACT RNA (z-Score)") +
  geom_bracket(xmin = DT_type$group2, 
               xmax = DT_type$group1, 
               y.position = DT_type$y.position, label = DT_type$labels, 
               tip.length = 0.05, step.increase = 0.1) +
  geom_boxplot(notch = TRUE, outliers = FALSE) + geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) + theme_science() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +   geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey40")
ggsave("type_v_PRKRA.eps", width = 4.4, height = 2, units = "in", dpi = 600, device=cairo_ps)





aov_EIF2AK2 <- aov(EIF2AK2 ~ call, RNAseq_BRCA_primary_normal_type)

summary(aov_EIF2AK2)

RNAseq_BRCA_primary_normal_type$call <- factor(as.factor(RNAseq_BRCA_primary_normal_type$call), levels = c("Normal", "TNBC", "ER+", "PR+", "HER2+", "ER+/PR+", "ER+/HER2+", "PR+/HER2+", "ER+/PR+/HER2+"))

DT_type <- DunnettTest(RNAseq_BRCA_primary_normal_type$EIF2AK2, RNAseq_BRCA_primary_normal_type$call)

DT_type <- as.data.frame(DT_type$Normal)

out <- strsplit(as.character(row.names(DT_type)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

DT_type$group1 <- out$V1
DT_type$group2 <- out$V2
DT_type$labels <- makeStars(DT_type$pval)
DT_type$y.position <- max(RNAseq_BRCA_primary_normal_type$EIF2AK2*1.1)

DT_type <- subset(DT_type, !labels == "ns")


ggplot(subset(RNAseq_BRCA_primary_normal_type, !is.na(call)), aes(as.factor(call), EIF2AK2)) + xlab("") + 
  scale_x_discrete(limits= c("Normal", "TNBC", "ER+", "PR+", "HER2+", "ER+/PR+", "ER+/HER2+", "PR+/HER2+", "ER+/PR+/HER2+"),
                   labels= c("Normal", "TNBC", "ER+", "PR+", "HER2+", "ER+\nPR+", "ER+\nHER2+", "PR+\nHER2+", "ER+/PR+\nHER2+")) +
  ylab("PKR RNA (z-Score)") +
  geom_bracket(xmin = DT_type$group2, 
               xmax = DT_type$group1, 
               y.position = DT_type$y.position, label = DT_type$labels, 
               tip.length = 0.05, step.increase = 0.1) +
  geom_boxplot(notch = TRUE, outliers = FALSE) + geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) + theme_science() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +   geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey40")
ggsave("type_v_EIF2AK2.eps", width = 4.4, height = 2, units = "in", dpi = 600, device=cairo_ps)



#read breast cancer clinical data
Clinical_BRCA <- readTCGA("/Users/kacottre/Desktop/TCGA/gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/BRCA.merged_only_clinical_clin_format.txt", dataType = 'clinical')

Clinical_BRCA <- readTCGA("/Users/kacottre/Desktop/TCGA/gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/BRCA.merged_only_clinical_clin_format.txt", dataType = 'clinical')

Clinical_BRCA$patient_barcode <- toupper(Clinical_BRCA$patient.bcr_patient_barcode)


#use the survivalTCGA function to extract survival data
BRCA_survival <- survivalTCGA(Clinical_BRCA)

BRCA_survival$times <- BRCA_survival$times/365

#merge survival data with primary tumor sample RNAseq
BRCA <- merge(BRCA_survival, RNAseq_BRCA_primary, by.x = "bcr_patient_barcode", by.y="patient_barcode")

#remove data with times <0
BRCA <- subset(BRCA, !times < 0)

#determine PRKRA expression cutoff using surv_cutpoint function

BRCA.surv_rnaseq.cut <- surv_cutpoint(
  BRCA,
  time = "times",
  event = "patient.vital_status",
  variables = "PRKRA"
)
summary(BRCA.surv_rnaseq.cut)
plot(BRCA.surv_rnaseq.cut, "PRKRA")

#score patients by PRKRA expression above or below cutoff
BRCA$PRKRA_score <- ifelse(BRCA$PRKRA < BRCA.surv_rnaseq.cut$cutpoint$cutpoint, "low", "high")

#fit the survival curve and plot
fit_PRKRA <- survfit(Surv(times, patient.vital_status)
                    ~ PRKRA_score , data = BRCA)


ggsurvplot(fit_PRKRA, data = BRCA, 
           size = 0.3,                 # change line size
           palette = palette2,# custom color palettes
           conf.int = TRUE,          # Add confidence interval
           pval = TRUE,              # Add p-value
           legend.title = "PACT\nExpression",
           legend = "right",
           legend.labs = c("High", "Low"),
           legend.size = 2,
           font.main = 8,
           font.x = 8,
           font.y = 8,
           font.tickslab = 8,
           title = "All BRCA",
           xlab = "Years",
           pval.size = 2.5,
           risk.table = FALSE,        # Add risk table
           font.family = "Arial",
           censor.size = 1.5,
           ggtheme = theme_science(), # Change ggplot2 theme
)
ggsave("PRKRA_cutoff.eps", width = 2.5, height = 1.8, units = "in", dpi = 600, device=cairo_ps)


tnbc_tumors <- subset(Annotation, TNBC == "YES")

#subset for TNBC
BRCA_tnbc <- subset(BRCA, bcr_patient_barcode %in% tnbc_tumors$BARCODE)

#determine PRKRA expression cutoff using surv_cutpoint function

BRCA_tnbc.surv_rnaseq.cut <- surv_cutpoint(
  BRCA_tnbc,
  time = "times",
  event = "patient.vital_status",
  variables = "PRKRA"
)
summary(BRCA_tnbc.surv_rnaseq.cut)
plot(BRCA_tnbc.surv_rnaseq.cut, "PRKRA")

#score patients by PRKRA expression above or below cutoff
BRCA_tnbc$PRKRA_score <- ifelse(BRCA_tnbc$PRKRA < BRCA_tnbc.surv_rnaseq.cut$cutpoint$cutpoint, "low", "high")

#fit the survival curve and plot
fit_PRKRA <- survfit(Surv(times, patient.vital_status)
                     ~ PRKRA_score , data = BRCA_tnbc)


ggsurvplot(fit_PRKRA, data = BRCA_tnbc, 
           size = 0.3,                 # change line size
           palette = palette2,# custom color palettes
           conf.int = TRUE,          # Add confidence interval
           pval = TRUE,              # Add p-value
           legend.title = "PACT\nExpression",
           legend = "right",
           legend.labs = c("High", "Low"),
           title = "TNBC",
           legend.size = 3,
           xlab = "Years",
           font.main = 8,
           font.x = 8,
           font.y = 8,
           font.tickslab = 8,
           pval.size = 2.5,
           risk.table = FALSE,        # Add risk table
           font.family = "Arial",
           censor.size = 1.5,
           ggtheme = theme_science() # Change ggplot2 theme
)
ggsave("PRKRA_cutoff_tnbc.eps", width = 2.5, height = 1.8, units = "in", dpi = 600, device=cairo_ps)





