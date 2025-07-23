library(ggplot2)
library(dplyr)
library(extrafont)
library(gridExtra)
library(RGraphics)
library(data.table)
library(ggpubr)
library(viridis)
library(scales)
library(gmodels)
library(DescTools)
library(MetBrewer)
library(svglite)
library(ggbeeswarm)
library(ggrepel)
library(forcats)
library(cowplot)
library(patchwork)
library(ggh4x)
library(gginnards)


# Load fonts
#extrafont::loadfonts(device="win")
#fonts()

#library(systemfonts)

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
          strip.background = element_rect(colour = "black", size = 0.5), legend.key = element_blank())
}



scaleFUN <- function(x) sprintf("%.2f", x)

makeStars <- function(x){
  stars <- c("***", "**", "*", "ns")
  vec <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}


cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbviridis <- c("#440154FF", "#31688EFF", "#35B779FF", "#E69F00")

palette <- met.brewer(name="Veronese", n=7, type="discrete")

met.brewer(name="Veronese", n=7, type="discrete")

palette2 <- c(palette[3], palette[7])
pallete4 <- c(palette[1], palette[3], palette[5],palette[7])
pallete5 <- c(palette[1], palette[4], palette[2],palette[5],palette[7])
pallete6 <- c(palette[1], palette[7], palette[2],palette[3], palette[6], palette[4])

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Define dataframe and omit blank values  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

filenames <- gsub("\\.csv$","", list.files(pattern="\\.csv$"))

for(i in filenames){
  assign(i, read.csv(paste(i, ".csv", sep="")))
}


####Immunoblot pact ko ----
#rename cell lines
immunoblot_ko$Cell_line <- gsub("BT549", "BT-549", immunoblot_ko$Cell_line)
immunoblot_ko$Cell_line <- gsub("MB453", "MDA-MB-453", immunoblot_ko$Cell_line)
immunoblot_ko$Cell_line <- gsub("MB468", "MDA-MB-468", immunoblot_ko$Cell_line)
immunoblot_ko$Cell_line <- factor(immunoblot_ko$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))


#group and summarise data
immunoblot_ko_sum <- dplyr::group_by(immunoblot_ko, Cell_line, Protein, sgRNA)

immunoblot_ko_sum<- dplyr::summarise(immunoblot_ko_sum, mean_expression = mean(Fold_change), sd_expression = sd(Fold_change))

#statistics for P-PKR in each cell line
#subset by cell line and P-PKR
#perform ANOVA and post-hoc Tukey
#create data.frame of Tukey results with additional information needed to annotate plots
ppkr_ko_hcc1806 <-subset(immunoblot_ko, Protein == "pPKR/PKR" & Cell_line == "HCC1806")

hcc1806 <- aov(Fold_change ~ sgRNA, data = ppkr_ko_hcc1806)
summary(hcc1806)

hcc1806ppkr <- TukeyHSD(hcc1806)
hcc1806ppkr <- as.data.frame(hcc1806ppkr$sgRNA[,1:4])

out <- strsplit(as.character(row.names(hcc1806ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

hcc1806ppkr$stars <- makeStars(hcc1806ppkr$`p adj`)

hcc1806ppkr$group1 <- out$V1
hcc1806ppkr$group2 <- out$V2
hcc1806ppkr$y <- max(ppkr_ko_hcc1806$Fold_change*1.1)
hcc1806ppkr$Cell_line <- "HCC1806"



ppkr_ko_bt549 <-subset(immunoblot_ko, Protein == "pPKR/PKR" & Cell_line == "BT-549")

bt549 <- aov(Fold_change ~ sgRNA, data = ppkr_ko_bt549)
summary(bt549)

bt549ppkr <- TukeyHSD(bt549)
bt549ppkr <- as.data.frame(bt549ppkr$sgRNA[,1:4])

out <- strsplit(as.character(row.names(bt549ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549ppkr$stars <- makeStars(bt549ppkr$`p adj`)

bt549ppkr$group1 <- out$V1
bt549ppkr$group2 <- out$V2
bt549ppkr$y <- max(ppkr_ko_bt549$Fold_change*1.1)
bt549ppkr$Cell_line = "BT-549"



ppkr_ko_mb453 <- subset(immunoblot_ko, Protein == "pPKR/PKR" & Cell_line == "MDA-MB-453")

mb453 <- aov(Fold_change ~ sgRNA, data = ppkr_ko_mb453)
summary(mb453)

mb453ppkr <- TukeyHSD(mb453)
mb453ppkr <- as.data.frame(mb453ppkr$sgRNA[,1:4])

out <- strsplit(as.character(row.names(mb453ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb453ppkr$stars <- makeStars(mb453ppkr$`p adj`)

mb453ppkr$group1 <- out$V1
mb453ppkr$group2 <- out$V2
mb453ppkr$y <- max(ppkr_ko_mb453$Fold_change*1.1)
mb453ppkr$Cell_line = "MDA-MB-453"


ppkr_ko_mb468 <- subset(immunoblot_ko, Protein == "pPKR/PKR" & Cell_line == "MDA-MB-468")

mb468 <- aov(Fold_change ~ sgRNA, data = ppkr_ko_mb468)
summary(mb468)

mb468ppkr <- TukeyHSD(mb468)
mb468ppkr <- as.data.frame(mb468ppkr$sgRNA[,1:4])

out <- strsplit(as.character(row.names(mb468ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb468ppkr$stars <- makeStars(mb468ppkr$`p adj`)

mb468ppkr$group1 <- out$V1
mb468ppkr$group2 <- out$V2
mb468ppkr$y <- max(ppkr_ko_mb468$Fold_change*1.1)
mb468ppkr$Cell_line = "MDA-MB-468"

#combine statistical analysis for all cell lines

tnbcppkr <- rbind(mb468ppkr, hcc1806ppkr, bt549ppkr, mb453ppkr)
tnbcppkr$y <- tnbcppkr$y
tnbcppkr <- subset(tnbcppkr, !stars == "ns")
tnbcppkr$Cell_line <- factor(tnbcppkr$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))

#colour of the facet labels
strip <- strip_themed(text_x = elem_list_text(colour = pallete4))

ggplot(subset(immunoblot_ko, Protein == "pPKR/PKR"), aes(sgRNA, Fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.3) + theme_science() + 
  scale_x_discrete(limits = c("C1", "NTA", "PACT1", "PACT2"), labels = c("sgC1", "sgNTA", "sgPACT-1", "sgPACT-2")) + labs(x = "", y = "Fold Change P-PKR/PKR", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 0, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none", fill = "none") + facet_grid2(cols = vars(Cell_line), strip = strip) +
  scale_colour_manual(values = pallete4) +
  geom_col(data = subset(immunoblot_ko_sum, Protein == "pPKR/PKR"), aes(sgRNA, mean_expression, fill = Cell_line), width = 0.5, alpha = 0.2, linewidth = 0.3) +
  geom_errorbar(data = subset(immunoblot_ko_sum, Protein == "pPKR/PKR"), aes(x = sgRNA, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2, linewidth = 0.3) + 
  stat_pvalue_manual(data = tnbcppkr, step.increase = 0.1, color = "Cell_line", label = "stars", y.position = "y", label.size = 3, linewidth = 0.3) +
  scale_fill_manual(values = pallete4)
ggsave("/Manuscripts/PACT/Western_plots/ppkr_ko_facet.tiff", height =2, width = 2.9, units = "in", dpi = 300)



####immunoblot ko ISR-NFKB ----


#group and summarise data
immunoblot_ko_effectors$log2FoldChange <- log2(immunoblot_ko_effectors$FoldChange)

immunoblot_ko_effectors$Gene <- factor(immunoblot_ko_effectors$Gene, levels = c("GADD34", "ATF3", "P-eIF2a", "C-PARP", "P-p65"))

immunoblot_ko_effectors_sum <- dplyr::group_by(immunoblot_ko_effectors, Gene, Sample)

immunoblot_ko_effectors_sum<- dplyr::summarise(immunoblot_ko_effectors_sum, mean_expression = mean(log2FoldChange), sd_expression = sd(log2FoldChange))

#t-tests for Fold Change and Sample
#adjust for multiple comparisons with FDR

isrnfkb <- immunoblot_ko_effectors %>% group_by(Gene) %>% 
  do(w = t.test(log2FoldChange ~ Sample, data=.)) %>% 
  summarise(Gene, p.value = w$p.value)
isrnfkb$padj <- p.adjust(isrnfkb$p.value, method = "fdr")


#add additional information needed for annotating plots
isrnfkb$stars <- makeStars(isrnfkb$padj)
isrnfkb$group1 <- "sgNTA"
isrnfkb$group2 <- "sgPACT-2"
isrnfkb$Gene <- factor(isrnfkb$Gene, levels = c("GADD34", "ATF3", "P-eIF2a", "C-PARP", "P-p65"))
isrnfkb <- subset(isrnfkb, padj < 0.05)
isrnfkb <-  isrnfkb[order(isrnfkb$Gene),]
isrnfkb$y <- c(9.7, 4.2, 6.4, 3.3)

strip <- strip_themed(text_x = elem_list_text(colour = pallete5))

ggplot(immunoblot_ko_effectors, aes(Sample, log2FoldChange, colour = Gene)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.3) + theme_science() + 
  scale_x_discrete(limits = c("sgNTA", "sgPACT-2")) + labs(x = "", y = "log2(Fold Change)", colour = "", fill = "") +
  geom_hline(yintercept = 0, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none", fill = "none") + facet_grid2(cols = vars(Gene), strip = strip) +
  scale_colour_manual(values = pallete5) +
  geom_col(data = immunoblot_ko_effectors_sum, aes(Sample, mean_expression, fill = Gene), width = 0.5, alpha = 0.2, linewidth = 0.3) +
  geom_errorbar(data = immunoblot_ko_effectors_sum, aes(x = Sample, ymax = mean_expression + sd_expression, ymin = ifelse(mean_expression - sd_expression <0, 0, mean_expression - sd_expression), y = mean_expression), width=0.2, linewidth = 0.3) + 
  stat_pvalue_manual(data = isrnfkb, color = "Gene", label = "stars", y.position = "y", label.size = 3, linewidth = 0.3) +
  scale_fill_manual(values = pallete5) + coord_cartesian(ylim = c(0, 10))
ggsave("/Manuscripts/PACT/Western_plots/isrnfkb_ko_facet.tiff", height =2, width = 2.6, units = "in", dpi = 300)



####Immunoblot rescue PACT oe ----
#immunoblot_rescue_pact <- subset(immunoblot_rescue_pact_pact, !Replicate == 2)
immunoblot_rescue_pact <- immunoblot_rescue_pact_pact

immunoblot_rescue_pact_sum <- dplyr::group_by(immunoblot_rescue_pact, Sample)

immunoblot_rescue_pact_sum<- dplyr::summarise(immunoblot_rescue_pact_sum, mean_area = mean(Fold_change), sd_area = sd(Fold_change))

immunoblot_rescue_pact_aov <- aov(Fold_change ~ Sample, data = immunoblot_rescue_pact)
summary(immunoblot_rescue_pact_aov)

immunoblot_rescue_pact$Sample <- factor(as.factor(immunoblot_rescue_pact$Sample), 
                                 levels = c("EV_PACT2", "EV_NTA", "PACT_NTA", "EAA_NTA", "AA_NTA", "PACT_PACT2", "EAA_PACT2", "AA_PACT2"))

str(immunoblot_rescue_pact)

rescue_pact_immunoblot <- DunnettTest(Fold_change ~ Sample, data = immunoblot_rescue_pact)
rescue_pact_immunoblot <- as.data.frame(rescue_pact_immunoblot$EV_PACT2[,1:4])

out <- strsplit(as.character(row.names(rescue_pact_immunoblot)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

rescue_pact_immunoblot$group1 <- out$V1
rescue_pact_immunoblot$group2 <- out$V2
rescue_pact_immunoblot$pval <- makeStars(rescue_pact_immunoblot$pval)
rescue_pact_immunoblot$y.position <- max(immunoblot_rescue_pact$Fold_change*1.1)

rescue_pact_immunoblot <- subset(rescue_pact_immunoblot, !pval == "ns")

ggplot(immunoblot_rescue_pact, aes(Sample, Fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.5) + theme_science() +
  scale_x_discrete(limits = c("EV_NTA", "PACT_NTA", "EAA_NTA", "AA_NTA", "EV_PACT2", "PACT_PACT2", "EAA_PACT2", "AA_PACT2"), 
                   labels = c("EV/sgNTA", "PACT/sgNTA", "PACT-EAA/sgNTA", "PACT-AA/sgNTA", "EV/sgPACT-2", "PACT/sgPACT-2", "PACT-EAA/sgPACT-2", "PACT-AA/sgPACT-2")) + 
  labs(x = "", y = "Fold Change P-PKR/PKR", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none") +
  geom_col(data = immunoblot_rescue_pact_sum, aes(Sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") + 
  geom_errorbar(data = immunoblot_rescue_pact_sum, aes(x = Sample, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area), width=0.2) +
  stat_pvalue_manual(data = rescue_pact_immunoblot, step.increase = 0.1, label = "pval", y.position = "y.position", label.size = 3, linewidth = 0.3)
ggsave("/Manuscripts/PACT/Western_plots/immunoblot_pact_pact_dunnetEV-PACT2.tiff", height = 2.5, width = 2, units = "in", dpi = 300)



####Immunoblot rescue d3 ----
immunoblot_rescue_d3 <- subset(immunoblot_rescue_pact_d3)

immunoblot_rescue_d3_sum <- dplyr::group_by(immunoblot_rescue_d3, Sample)

immunoblot_rescue_d3_sum<- dplyr::summarise(immunoblot_rescue_d3_sum, mean_area = mean(Fold_change), sd_area = sd(Fold_change))


#statistical analysis for P-PKR
#perform ANOVA

immunoblot_rescue_d3_aov <- aov(Fold_change ~ Sample, data = immunoblot_rescue_d3)
summary(immunoblot_rescue_d3_aov)

#set factor levels for samples
immunoblot_rescue_d3$Sample <- factor(as.factor(immunoblot_rescue_d3$Sample), 
                                        levels = c("PACT2_EV", "NTA_EV", "NTA_PACT", "NTA_PACTd3", "NTA_PACTd3GST",
                                                   "PACT2_PACT", "PACT2_PACTd3", "PACT2_PACTd3GST"))

str(immunoblot_rescue_d3)

#perform Dunnett's Test
#add additional information for annotating plots
rescue_pact_immunoblot <- DunnettTest(Fold_change ~ Sample, data = immunoblot_rescue_d3)
rescue_pact_immunoblot <- as.data.frame(rescue_pact_immunoblot$PACT2_EV[,1:4])

out <- strsplit(as.character(row.names(rescue_pact_immunoblot)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

rescue_pact_immunoblot$group1 <- out$V1
rescue_pact_immunoblot$group2 <- out$V2
rescue_pact_immunoblot$pval <- makeStars(rescue_pact_immunoblot$pval)
rescue_pact_immunoblot$y.position <- max(immunoblot_rescue_d3$Fold_change*1.1)

rescue_pact_immunoblot <- subset(rescue_pact_immunoblot, !pval == "ns")

ggplot(immunoblot_rescue_d3, aes(Sample, Fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.5) + theme_science() +
  scale_x_discrete(limits = c("NTA_EV", "NTA_PACT", "NTA_PACTd3", "NTA_PACTd3GST",
                              "PACT2_EV", "PACT2_PACT", "PACT2_PACTd3", "PACT2_PACTd3GST"), 
                   labels = c("EV/sgNTA", "PACT/sgNTA", "PACTd3/sgNTA", "PACTd3-GST/sgNTA",
                              "EV/sgPACT-2", "PACT/sgPACT-2", "PACTd3/sgPACT-2", "PACTd3-GST/sgPACT-2")) + 
  labs(x = "", y = "Fold Change P-PKR/PKR", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none") +
  geom_col(data = immunoblot_rescue_d3_sum, aes(Sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") + 
  geom_errorbar(data = immunoblot_rescue_d3_sum, aes(x = Sample, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area), width=0.2) +
  stat_pvalue_manual(data = rescue_pact_immunoblot, step.increase = 0.1, label = "pval", y.position = "y.position", label.size = 3, linewidth = 0.3)
ggsave("/Manuscripts/PACT/Western_plots/immunoblot_pact_d3_dunnetEV-PACT2.tiff", height = 2.5, width = 2, units = "in", dpi = 300)


####same as above for GST only control

immunoblot_gst_control_sum <- dplyr::group_by(immunoblot_gst_control, Sample)

immunoblot_gst_control_sum<- dplyr::summarise(immunoblot_gst_control_sum, mean_area = mean(Fold_Change), sd_area = sd(Fold_Change))

immunoblot_gst_control_aov <- aov(Fold_Change ~ Sample, data = immunoblot_gst_control)
summary(immunoblot_gst_control_aov)

immunoblot_gst_control$Sample <- factor(as.factor(immunoblot_gst_control$Sample), 
                                      levels = c("P2_EV", "NTA_EV", "NTA_GST", "P2_GST"))

str(immunoblot_gst_control)

gst_pact_immunoblot <- DunnettTest(Fold_Change ~ Sample, data = immunoblot_gst_control)
gst_pact_immunoblot <- as.data.frame(gst_pact_immunoblot$P2_EV[,1:4])

out <- strsplit(as.character(row.names(gst_pact_immunoblot)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

gst_pact_immunoblot$group1 <- out$V1
gst_pact_immunoblot$group2 <- out$V2
gst_pact_immunoblot$pval <- makeStars(gst_pact_immunoblot$pval)
gst_pact_immunoblot$y.position <- max(immunoblot_gst_control$Fold_Change*1.1)

gst_pact_immunoblot <- subset(gst_pact_immunoblot, !pval == "ns")

ggplot(immunoblot_gst_control, aes(Sample, Fold_Change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.5) + theme_science() +
  scale_x_discrete(limits = c("NTA_EV", "NTA_GST", "P2_EV", "P2_GST"), 
                   labels = c("EV/sgNTA", "GST/sgNTA", "EV/sgPACT-2", "GST/sgPACT-2")) + 
  labs(x = "", y = "Fold Change P-PKR/PKR", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none") +
  geom_col(data = immunoblot_gst_control_sum, aes(Sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") + 
  geom_errorbar(data = immunoblot_gst_control_sum, aes(x = Sample, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area), width=0.2) +
  stat_pvalue_manual(data = gst_pact_immunoblot, step.increase = 0.1, label = "pval", y.position = "y.position", label.size = 3, linewidth = 0.3)
ggsave("/Manuscripts/PACT/Western_plots/immunoblot_pact_gst_dunnetEV-PACT2.tiff", height = 2, width = 1.5, units = "in", dpi = 300)




####immunoblot_baseline ----
immunoblot_baseline$Protein <- gsub("pPKR", "P-PKR/PKR", immunoblot_baseline$Protein)

immunoblot_baseline_sum <- dplyr::group_by(immunoblot_baseline, Cell_line, Protein)

immunoblot_baseline_sum<- dplyr::summarise(immunoblot_baseline_sum, mean_area = mean(Fold_Change), sd_area = sd(Fold_Change))

ggplot(immunoblot_baseline, aes(Cell_line, Fold_Change, colour = Protein)) +
  geom_jitter(width = 0.12, size = 0.5) + theme_science() +
  scale_x_discrete(limits = c("MCF7", "SKBR3", "BT549", "MB453", "MB468", "HCC1806"), 
                   labels = c("MCF-7", "SK-BR-3", "BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806")) + 
  labs(x = "", y = "Relative Abundance", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none") + facet_grid2(cols = vars(Protein)) + scale_colour_manual(values = pallete4) +
  scale_fill_manual(values = pallete4) +
  geom_col(data = immunoblot_baseline_sum, aes(Cell_line, mean_area, colour = Protein, fill = Protein), width = 0.5, alpha = 0.3) + 
  geom_errorbar(data = immunoblot_baseline_sum, aes(x = Cell_line, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area, colour = Protein), width=0.2) 
ggsave("/Manuscripts/PACT/Western_plots/immunoblot_baseline.tiff", height = 2.2, width = 4, units = "in", dpi = 300)




###rescue PACT - ADAR ----

#immunoblot_rescue_adar <- subset(immunoblot_rescue_pact_adar, !Replicate == 2)
immunoblot_rescue_adar <- immunoblot_rescue_pact_adar

immunoblot_rescue_adar_sum <- dplyr::group_by(immunoblot_rescue_adar, Sample)

immunoblot_rescue_adar_sum<- dplyr::summarise(immunoblot_rescue_adar_sum, mean_area = mean(Fold_change), sd_area = sd(Fold_change))

immunoblot_rescue_adar_aov <- aov(Fold_change ~ Sample, data = immunoblot_rescue_adar)
summary(immunoblot_rescue_adar_aov)

immunoblot_rescue_adar$Sample <- factor(as.factor(immunoblot_rescue_adar$Sample), 
                                        levels = c("EV_PACT2", "EV_NTA", "p110_NTA", "p150_NTA", "p110_PACT2", "p150_PACT2"))

str(immunoblot_rescue_adar)

rescue_adar_immunoblot <- DunnettTest(Fold_change ~ Sample, data = immunoblot_rescue_adar)
rescue_adar_immunoblot <- as.data.frame(rescue_adar_immunoblot$EV_PACT2[,1:4])

out <- strsplit(as.character(row.names(rescue_adar_immunoblot)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

rescue_adar_immunoblot$group1 <- out$V1
rescue_adar_immunoblot$group2 <- out$V2
rescue_adar_immunoblot$pval <- makeStars(rescue_adar_immunoblot$pval)
rescue_adar_immunoblot$y.position <- max(immunoblot_rescue_adar$Fold_change*1.1)

rescue_adar_immunoblot <- subset(rescue_adar_immunoblot, !pval == "ns")

ggplot(immunoblot_rescue_adar, aes(Sample, Fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.5) + theme_science() +
  scale_x_discrete(limits = c("EV_NTA", "p110_NTA", "p150_NTA", "EV_PACT2", "p110_PACT2", "p150_PACT2"), 
                   labels = c("EV/sgNTA", "p110/sgNTA", "p150/sgNTA", "EV/sgPACT-2", "p110/sgPACT-2", "p150/sgPACT-2")) + 
  labs(x = "", y = "Fold Change pPKR/PKR", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none") +
  geom_col(data = immunoblot_rescue_adar_sum, aes(Sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") + 
  geom_errorbar(data = immunoblot_rescue_adar_sum, aes(x = Sample, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area), width=0.2) +
stat_pvalue_manual(data = rescue_adar_immunoblot, step.increase = 0.1, label = "pval", y.position = "y.position", label.size = 3, linewidth = 0.3)
ggsave("/Manuscripts/PACT/Western_plots/rescue_pact_adar_dunnetEV-PACT2.tiff", height = 2.5, width = 1.8, units = "in", dpi = 300)




####ctg PACT ko ----
ctg_ko$Cell_line <- gsub("BT549", "BT-549", ctg_ko$Cell_line)
ctg_ko$Cell_line <- gsub("MB453", "MDA-MB-453", ctg_ko$Cell_line)
ctg_ko$Cell_line <- gsub("MB468", "MDA-MB-468", ctg_ko$Cell_line)
ctg_ko$Cell_line <- factor(ctg_ko$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))



#group and summarise data
ctg_ko_sum <- dplyr::group_by(ctg_ko, Cell_line, sgRNA)

ctg_ko_sum<- dplyr::summarise(ctg_ko_sum, mean_area = mean(Fold_change), sd_area = sd(Fold_change))

#statistics for CTG in each cell line
#subset by cell line
#perform ANOVA and post-hoc Tukey
#create data.frame of Tukey results with additional information needed to annotate plots

ctg_ko_hcc1806 <-subset(ctg_ko, Cell_line == "HCC1806")

hcc1806 <- aov(Fold_change ~ sgRNA, data = ctg_ko_hcc1806)
summary(hcc1806)

hcc1806ctg <- TukeyHSD(hcc1806)
hcc1806ctg <- as.data.frame(hcc1806ctg$sgRNA[,1:4])

out <- strsplit(as.character(row.names(hcc1806ctg)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

hcc1806ctg$stars <- makeStars(hcc1806ctg$`p adj`)

hcc1806ctg$group1 <- out$V1
hcc1806ctg$group2 <- out$V2
hcc1806ctg$y <- max(ctg_ko_hcc1806$Fold_change*1.1)
hcc1806ctg$Cell_line <- "HCC1806"



ctg_ko_bt549 <-subset(ctg_ko, Cell_line == "BT-549")

bt549 <- aov(Fold_change ~ sgRNA, data = ctg_ko_bt549)
summary(bt549)

bt549ctg <- TukeyHSD(bt549)
bt549ctg <- as.data.frame(bt549ctg$sgRNA[,1:4])

out <- strsplit(as.character(row.names(bt549ctg)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549ctg$stars <- makeStars(bt549ctg$`p adj`)

bt549ctg$group1 <- out$V1
bt549ctg$group2 <- out$V2
bt549ctg$y <- max(ctg_ko_bt549$Fold_change*1.1)
bt549ctg$Cell_line = "BT-549"



ctg_ko_mb453 <-subset(ctg_ko, Cell_line == "MDA-MB-453")

mb453 <- aov(Fold_change ~ sgRNA, data = ctg_ko_mb453)
summary(mb453)

mb453ctg <- TukeyHSD(mb453)
mb453ctg <- as.data.frame(mb453ctg$sgRNA[,1:4])

out <- strsplit(as.character(row.names(mb453ctg)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb453ctg$stars <- makeStars(mb453ctg$`p adj`)

mb453ctg$group1 <- out$V1
mb453ctg$group2 <- out$V2
mb453ctg$y <- max(ctg_ko_mb453$Fold_change*1.1)
mb453ctg$Cell_line = "MDA-MB-453"


ctg_ko_mb468 <-subset(ctg_ko, Cell_line == "MDA-MB-468")

mb468 <- aov(Fold_change ~ sgRNA, data = ctg_ko_mb468)
summary(mb468)

mb468ctg <- TukeyHSD(mb468)
mb468ctg <- as.data.frame(mb468ctg$sgRNA[,1:4])

out <- strsplit(as.character(row.names(mb468ctg)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb468ctg$stars <- makeStars(mb468ctg$`p adj`)

mb468ctg$group1 <- out$V1
mb468ctg$group2 <- out$V2
mb468ctg$y <- max(ctg_ko_mb468$Fold_change*1.1)
mb468ctg$Cell_line = "MDA-MB-468"

#combine stats from each cell line
tnbcctg <- rbind(mb468ctg, hcc1806ctg, bt549ctg, mb453ctg)
tnbcctg <- subset(tnbcctg, !stars == "ns")
tnbcctg$Cell_line <- factor(tnbcctg$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))

strip <- strip_themed(text_x = elem_list_text(colour = pallete4))

ggplot(ctg_ko, aes(sgRNA, Fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.3) + theme_science() + 
  scale_x_discrete(limits = c("C1", "NTA", "PACT1", "PACT2"), labels = c("sgC1", "sgNTA", "sgPACT-1", "sgPACT-2")) + labs(x = "", y = "Relative Viability", colour = "", fill = "") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75), expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none", fill = "none") + facet_grid2(cols = vars(Cell_line), strip = strip) + 
  scale_colour_manual(values = pallete4) +
  geom_col(data = ctg_ko_sum, aes(sgRNA, mean_area, fill = Cell_line), width = 0.5, alpha = 0.2, linewidth = 0.3) +
  geom_errorbar(data = ctg_ko_sum, aes(x = sgRNA, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area), width=0.2, linewidth = 0.3) +
  stat_pvalue_manual(data = tnbcctg, step.increase = 0.1, step.group.by = "Cell_line", color = "Cell_line", label = "stars", y.position = "y", label.size = 3, linewidth = 0.3) +
  scale_fill_manual(values = pallete4)
ggsave("/Manuscripts/PACT/CTG_plots/ctg_ko_facet.tiff", height =2, width = 2.9, units = "in", dpi = 300)



####ctg PACT and ADAR1 depletion ----
ctg_dko <- ctg_pact_adar_dko
ctg_dko$Cell_line <- gsub("BT549", "BT-549", ctg_dko$Cell_line)
ctg_dko$Cell_line <- gsub("MB453", "MDA-MB-453", ctg_dko$Cell_line)
ctg_dko$Cell_line <- factor(ctg_dko$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))

ctg_dko$shRNA <- grepl("A6", ctg_dko$sgRNA)

#group and summarise data
ctg_dko_sum <- dplyr::group_by(ctg_dko, Cell_line, sgRNA)

ctg_dko_sum<- dplyr::summarise(ctg_dko_sum, mean_area = mean(Fold_change), sd_area = sd(Fold_change))

ctg_dko_sum$shRNA <- grepl("A6", ctg_dko_sum$sgRNA)

#statistics for CTG in each cell line
#subset by cell line
#perform ANOVA and post-hoc Tukey
#create data.frame of Tukey results with additional information needed to annotate plots

ctg_dko_bt549 <-subset(ctg_dko, Cell_line == "BT-549")

bt549 <- aov(Fold_change ~ sgRNA, data = ctg_dko_bt549)
summary(bt549)

bt549ctg <- TukeyHSD(bt549)
bt549ctg <- as.data.frame(bt549ctg$sgRNA[,1:4])

out <- strsplit(as.character(row.names(bt549ctg)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549ctg$stars <- makeStars(bt549ctg$`p adj`)

bt549ctg$group1 <- out$V1
bt549ctg$group2 <- out$V2
bt549ctg$y <- max(ctg_dko_bt549$Fold_change*1.1)
bt549ctg$Cell_line = "BT-549"

ctg_dko_mb453 <-subset(ctg_dko, Cell_line == "MDA-MB-453")

mb453 <- aov(Fold_change ~ sgRNA, data = ctg_dko_mb453)
summary(mb453)

mb453ctg <- TukeyHSD(mb453)
mb453ctg <- as.data.frame(mb453ctg$sgRNA[,1:4])

out <- strsplit(as.character(row.names(mb453ctg)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb453ctg$stars <- makeStars(mb453ctg$`p adj`)

mb453ctg$group1 <- out$V1
mb453ctg$group2 <- out$V2
mb453ctg$y <- max(ctg_dko_mb453$Fold_change*1.1)
mb453ctg$Cell_line = "MDA-MB-453"

#combine stats for each cell line

dkoctg <- rbind(bt549ctg, mb453ctg)
dkoctg <- subset(dkoctg, !stars == "ns")

dkoctg$Cell_line <- factor(dkoctg$Cell_line, levels = c("BT-549", "MDA-MB-453"))

p <- ggplot(ctg_dko, aes(sgRNA, Fold_change, colour = shRNA)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.3) + theme_science() + 
  scale_x_discrete(limits = c("C1_SCR", "NTA_SCR", "PACT1_SCR", "PACT2_SCR", "C1_A6", "NTA_A6", "PACT1_A6", "PACT2_A6"), 
                   labels = c("sgC1", "sgNTA", "sgPACT1", "sgPACT2", "sgC1", "sgNTA", "sgPACT1", "sgPACT2")) + 
  labs(x = "", y = "Relative Viability", colour = "", fill = "") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75), expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + facet_grid2(cols = vars(Cell_line)) + 
  geom_col(data = ctg_dko_sum, aes(sgRNA, mean_area, colour = shRNA, fill = shRNA), width = 0.5, linewidth = 0.3, alpha = 0.2) +
  geom_errorbar(data = ctg_dko_sum, aes(x = sgRNA, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area, colour = shRNA), width = 0.2, linewidth = 0.3) +
  stat_pvalue_manual(data = dkoctg, step.increase = 0.1, step.group.by = "Cell_line", label = "stars", y.position = "y") +
  scale_colour_manual(values = palette2, labels = c("shSCR", "shADAR1")) + 
  scale_fill_manual(values = palette2, labels = c("shSCR", "shADAR1"))
move_layers(p, "GeomPoint", position = "top")
ggsave("/Manuscripts/PACT/CTG_plots/ctg_pact_adar_dko_facet_all.tiff", height = 2.1, width = 2.6, units = "in", dpi = 300)

#remove comparisons between the two control sgRNAs or the two PACT sgRNAs to make the plot less cluttered 

dkoctg <- subset(dkoctg, !grepl("PACT1.*PACT2|PACT2.*PACT1", paste(dkoctg$group1, dkoctg$group2)))
dkoctg <- subset(dkoctg, !grepl("C1.*NTA|NTA.*C1", paste(dkoctg$group1, dkoctg$group2)))

p <- ggplot(ctg_dko, aes(sgRNA, Fold_change, colour = shRNA)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.3) + theme_science() + 
  scale_x_discrete(limits = c("C1_SCR", "NTA_SCR", "PACT1_SCR", "PACT2_SCR", "C1_A6", "NTA_A6", "PACT1_A6", "PACT2_A6"), 
                   labels = c("sgC1", "sgNTA", "sgPACT1", "sgPACT2", "sgC1", "sgNTA", "sgPACT1", "sgPACT2")) + 
  labs(x = "", y = "Relative Viability", colour = "", fill = "") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75), expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + facet_grid2(cols = vars(Cell_line)) + 
  geom_col(data = ctg_dko_sum, aes(sgRNA, mean_area, colour = shRNA, fill = shRNA), width = 0.5, linewidth = 0.3, alpha = 0.2) +
  geom_errorbar(data = ctg_dko_sum, aes(x = sgRNA, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area, colour = shRNA), width = 0.2, linewidth = 0.3) +
  stat_pvalue_manual(data = dkoctg, step.increase = 0.1, step.group.by = "Cell_line", label = "stars", y.position = "y") +
  scale_colour_manual(values = palette2, labels = c("shSCR", "shADAR1")) + 
  scale_fill_manual(values = palette2, labels = c("shSCR", "shADAR1"))
move_layers(p, "GeomPoint", position = "top")
ggsave("/Manuscripts/PACT/CTG_plots/ctg_pact_adar_dko_facet.tiff", height = 2.1, width = 2.6, units = "in", dpi = 300)


####CTG PKR PACT dko ----

ctg_dko_pkr <- ctg_pact_pkr_dko

#group and summarise data
ctg_dko_pkr_sum <- dplyr::group_by(ctg_dko_pkr, Cell_line, Sample)

ctg_dko_pkr_sum<- dplyr::summarise(ctg_dko_pkr_sum, mean_area = mean(Fold_change), sd_area = sd(Fold_change))

ctg_dko_pkr_aov <- aov(Fold_change ~ Sample, data = ctg_dko_pkr)
summary(ctg_dko_pkr_aov)

ctg_dko_pkr$Sample <- factor(as.factor(ctg_dko_pkr$Sample), 
                             levels = c("sgC1_sgPACT","sgC1_sgNTA", "sgPKR_sgNTA", "sgPKR_sgPACT2"))

str(ctg_dko_pkr)

pkr_pact_ctg <- DunnettTest(Fold_change ~ Sample, data = ctg_dko_pkr)
pkr_pact_ctg <- as.data.frame(pkr_pact_ctg$sgC1_sgPACT[,1:4])

out <- strsplit(as.character(row.names(pkr_pact_ctg)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

pkr_pact_ctg$group1 <- out$V1
pkr_pact_ctg$group2 <- out$V2
pkr_pact_ctg$pval <- makeStars(pkr_pact_ctg$pval)
pkr_pact_ctg$y.position <- max(ctg_dko_pkr$Fold_change*1.1)

pkr_pact_ctg <- subset(pkr_pact_ctg, !pval == "ns")

ggplot(ctg_dko_pkr, aes(Sample, Fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.5) + theme_science() + 
  scale_x_discrete(limits = c("sgC1_sgNTA", "sgPKR_sgNTA", "sgC1_sgPACT", "sgPKR_sgPACT2"), 
                   labels = c("sgC1/sgNTA", "sgPKR/sgNTA", "sgC1/sgPACT-2", "sgPKR/sgPACT-2")) + 
  labs(x = "", y = "Relative Viability", colour = "", fill = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5), expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none") +
  geom_col(data = ctg_dko_pkr_sum, aes(Sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") + 
  geom_errorbar(data = ctg_dko_pkr_sum, aes(x = Sample, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area), width=0.2) +
  stat_pvalue_manual(data = pkr_pact_ctg, step.increase = 0.1, label = "pval", y.position = "y.position", label.size = 3, linewidth = 0.3)
ggsave("/Manuscripts/PACT/CTG_plots/ctg_pact_pkr_dko_dunnet.tiff", height = 2.4, width = 1.5, units = "in")


#########CTG rescue ----

#group and summarise data
ctg_rescue_pact <- subset(ctg_rescue, Sample %in% c("EV_sgNTA", "PACT_sgNTA", "PACTeaa_sgNTA", "PACTaa_sgNTA", "EV_sgPACT2", "PACT_sgPACT2", "PACTeaa_sgPACT2", "PACTaa_sgPACT2"))

ctg_rescue_pact_sum <- dplyr::group_by(ctg_rescue_pact, Cell_line, Sample)

ctg_rescue_pact_sum<- dplyr::summarise(ctg_rescue_pact_sum, mean_area = mean(Fold_change), sd_area = sd(Fold_change))

ctg_rescue_adar <- subset(ctg_rescue, Sample %in% c("EV_sgNTA", "p110_sgNTA", "p150_sgNTA", "EV_sgPACT2", "p110_sgPACT2", "p150_sgPACT2"))

ctg_rescue_adar_sum <- dplyr::group_by(ctg_rescue_adar, Cell_line, Sample)

ctg_rescue_adar_sum<- dplyr::summarise(ctg_rescue_adar_sum, mean_area = mean(Fold_change), sd_area = sd(Fold_change))

ctg_rescue_pact_aov <- aov(Fold_change ~ Sample, data = ctg_rescue_pact)
summary(ctg_rescue_pact_aov)

ctg_rescue_pact$Sample <- factor(as.factor(ctg_rescue_pact$Sample), 
                                 levels = c("EV_sgPACT2", "EV_sgNTA", "PACT_sgNTA", "PACTeaa_sgNTA", "PACTaa_sgNTA", "PACT_sgPACT2", "PACTeaa_sgPACT2", "PACTaa_sgPACT2"))

str(ctg_rescue_pact)

rescue_pact_ctg <- DunnettTest(Fold_change ~ Sample, data = ctg_rescue_pact)
rescue_pact_ctg <- as.data.frame(rescue_pact_ctg$EV_sgPACT2[,1:4])

out <- strsplit(as.character(row.names(rescue_pact_ctg)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

rescue_pact_ctg$group1 <- out$V1
rescue_pact_ctg$group2 <- out$V2
rescue_pact_ctg$pval <- makeStars(rescue_pact_ctg$pval)
rescue_pact_ctg$y.position <- max(ctg_rescue_pact$Fold_change*1.1)

rescue_pact_ctg <- subset(rescue_pact_ctg, !pval == "ns")

ggplot(ctg_rescue_pact, aes(Sample, Fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.5) + theme_science() +
  scale_x_discrete(limits = c("EV_sgNTA", "PACT_sgNTA", "PACTeaa_sgNTA", "PACTaa_sgNTA", "EV_sgPACT2", "PACT_sgPACT2", "PACTeaa_sgPACT2", "PACTaa_sgPACT2"), 
                   labels = c("EV/sgNTA", "PACT/sgNTA", "PACT-EAA/sgNTA", "PACT-AA/sgNTA", "EV/sgPACT-2", "PACT/sgPACT-2", "PACT-EAA/sgPACT-2", "PACT-AA/sgPACT-2")) + 
  labs(x = "", y = "Relative Viability", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none") +
  geom_col(data = ctg_rescue_pact_sum, aes(Sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") + 
  geom_errorbar(data = ctg_rescue_pact_sum, aes(x = Sample, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area), width=0.2) +
  stat_pvalue_manual(data = rescue_pact_ctg, step.increase = 0.1, label = "pval", y.position = "y.position", label.size = 3, linewidth = 0.3)
ggsave("/Manuscripts/PACT/CTG_plots/ctg_pact_pact_dunnetEV-PACT2.tiff", height = 2.5, width = 2, units = "in", dpi = 300)

rescue_pact_ctg <- NULL


ctg_rescue_d3_sum <- dplyr::group_by(ctg_rescue_d3, Sample)

ctg_rescue_d3_sum<- dplyr::summarise(ctg_rescue_d3_sum, mean_area = mean(Fold_change), sd_area = sd(Fold_change))

ctg_rescue_d3_aov <- aov(Fold_change ~ Sample, data = ctg_rescue_d3)
summary(ctg_rescue_d3_aov)

ctg_rescue_d3$Sample <- factor(as.factor(ctg_rescue_d3$Sample), 
                               levels = c("PACT2_EV", "NTA_EV", "NTA_PACT", "NTA_PACTd3", "NTA_PACTd3GST",
                                          "PACT2_PACT", "PACT2_PACTd3", "PACT2_PACTd3GST"))
str(ctg_rescue_d3)

rescue_pact_ctg <- DunnettTest(Fold_change ~ Sample, data = ctg_rescue_d3)
rescue_pact_ctg <- as.data.frame(rescue_pact_ctg$PACT2_EV[,1:4])

out <- strsplit(as.character(row.names(rescue_pact_ctg)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

rescue_pact_ctg$group1 <- out$V1
rescue_pact_ctg$group2 <- out$V2
rescue_pact_ctg$pval <- makeStars(rescue_pact_ctg$pval)
rescue_pact_ctg$y.position <- max(ctg_rescue_d3$Fold_change*1.1)

rescue_pact_ctg <- subset(rescue_pact_ctg, !pval == "ns")

ggplot(ctg_rescue_d3, aes(Sample, Fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.5) + theme_science() +
  scale_x_discrete(limits = c("NTA_EV", "NTA_PACT", "NTA_PACTd3", "NTA_PACTd3GST",
                              "PACT2_EV", "PACT2_PACT", "PACT2_PACTd3", "PACT2_PACTd3GST"), 
                   labels = c("EV/sgNTA", "PACT/sgNTA", "PACTd3/sgNTA", "PACTd3-GST/sgNTA",
                              "EV/sgPACT-2", "PACT/sgPACT-2", "PACTd3/sgPACT-2", "PACTd3-GST/sgPACT-2")) +  labs(x = "", y = "Relative Viability", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none") +
  geom_col(data = ctg_rescue_d3_sum, aes(Sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") + 
  geom_errorbar(data = ctg_rescue_d3_sum, aes(x = Sample, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area), width=0.2) +
  stat_pvalue_manual(data = rescue_pact_ctg, step.increase = 0.1, label = "pval", y.position = "y.position", label.size = 3, linewidth = 0.3)
ggsave("/Manuscripts/PACT/CTG_plots/ctg_pact_d3_dunnetEV-PACT2.tiff", height = 2.5, width = 2, units = "in", dpi = 300)

rescue_pact_ctg <- NULL



ctg_rescue_adar_aov <- aov(Fold_change ~ Sample, data = ctg_rescue_adar)
summary(ctg_rescue_adar_aov)

ctg_rescue_adar$Sample <- factor(as.factor(ctg_rescue_adar$Sample), 
                                 levels = c("EV_sgPACT2", "EV_sgNTA", "p110_sgNTA", "p150_sgNTA", "p110_sgPACT2", "p150_sgPACT2"))

str(ctg_rescue_adar)

rescue_adar_ctg <- DunnettTest(Fold_change ~ Sample, data = ctg_rescue_adar)
rescue_adar_ctg <- as.data.frame(rescue_adar_ctg$EV_sgPACT2[,1:4])

out <- strsplit(as.character(row.names(rescue_adar_ctg)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

rescue_adar_ctg$group1 <- out$V1
rescue_adar_ctg$group2 <- out$V2
rescue_adar_ctg$pval <- makeStars(rescue_adar_ctg$pval)
rescue_adar_ctg$y.position <- max(ctg_rescue_adar$Fold_change*1.1)

rescue_adar_ctg <- subset(rescue_adar_ctg, !pval == "ns")

ggplot(ctg_rescue_adar, aes(Sample, Fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.5) + theme_science() + 
  scale_x_discrete(limits = c("EV_sgNTA", "p110_sgNTA", "p150_sgNTA", "EV_sgPACT2", "p110_sgPACT2", "p150_sgPACT2"), 
                   labels = c("EV/sgNTA", "p110/sgNTA", "p150/sgNTA", "EV/sgPACT-2", "p110/sgPACT-2", "p150/sgPACT-2")) + 
  labs(x = "", y = "Relative Viability", colour = "", fill = "") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75), expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none") +
  geom_col(data = ctg_rescue_adar_sum, aes(Sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") + 
  geom_errorbar(data = ctg_rescue_adar_sum, aes(x = Sample, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area), width=0.2) +
  stat_pvalue_manual(data = rescue_adar_ctg, step.increase = 0.1, label = "pval", y.position = "y.position", label.size = 3, linewidth = 0.3)
ggsave("/Manuscripts/PACT/CTG_plots/ctg_pact_ADAR_dunnetEV-PACT2.tiff", height = 2.5, width = 1.8, units = "in")



ctg_gst_control_sum <- dplyr::group_by(ctg_gst_control, Sample)

ctg_gst_control_sum<- dplyr::summarise(ctg_gst_control_sum, mean_area = mean(Fold_change), sd_area = sd(Fold_change))

ctg_gst_control_aov <- aov(Fold_change ~ Sample, data = ctg_gst_control)
summary(ctg_gst_control_aov)

ctg_gst_control$Sample <- factor(as.factor(ctg_gst_control$Sample), 
                               levels = c("P2_EV", "NTA_EV", "NTA_GST", "P2_GST"))
str(ctg_gst_control)

rescue_pact_ctg <- DunnettTest(Fold_change ~ Sample, data = ctg_gst_control)
rescue_pact_ctg <- as.data.frame(rescue_pact_ctg$P2_EV[,1:4])

out <- strsplit(as.character(row.names(rescue_pact_ctg)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

rescue_pact_ctg$group1 <- out$V1
rescue_pact_ctg$group2 <- out$V2
rescue_pact_ctg$pval <- makeStars(rescue_pact_ctg$pval)
rescue_pact_ctg$y.position <- max(ctg_gst_control$Fold_change*1.1)

rescue_pact_ctg <- subset(rescue_pact_ctg, !pval == "ns")

ggplot(ctg_gst_control, aes(Sample, Fold_change)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.5) + theme_science() +
  scale_x_discrete(limits = c("NTA_EV", "NTA_GST", "P2_EV", "P2_GST"), 
                   labels = c("EV/sgNTA", "GST/sgNTA", "EV/sgPACT-2", "GST/sgPACT-2")) +  labs(x = "", y = "Relative Viability", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none") +
  geom_col(data = ctg_gst_control_sum, aes(Sample, mean_area), width = 0.5, fill = "white", alpha = 0, colour = "black") + 
  geom_errorbar(data = ctg_gst_control_sum, aes(x = Sample, ymax = mean_area + sd_area, ymin = mean_area - sd_area, y = mean_area), width=0.2) +
  stat_pvalue_manual(data = rescue_pact_ctg, step.increase = 0.1, label = "pval", y.position = "y.position", label.size = 3, linewidth = 0.3)
ggsave("/Manuscripts/PACT/CTG_plots/ctg_pact_dunnetEV-PACT2.tiff", height = 2, width = 1.5, units = "in", dpi = 300)

rescue_pact_ctg <- NULL





####Immunoblot PACT/ADAR1 dko ----
immunoblot_dko$Cell_line <- gsub("BT549", "BT-549", immunoblot_dko$Cell_line)
immunoblot_dko$Cell_line <- gsub("MB453", "MDA-MB-453", immunoblot_dko$Cell_line)
immunoblot_dko$Cell_line <- factor(immunoblot_dko$Cell_line, levels = c("BT-549", "MDA-MB-453"))

immunoblot_dko$shRNA <- grepl("A6", immunoblot_dko$Sample)

#group and summarise data
immunoblot_dko_sum <- dplyr::group_by(immunoblot_dko, Cell_line, Protein, Sample)

immunoblot_dko_sum<- dplyr::summarise(immunoblot_dko_sum, mean_expression = mean(Fold_change), sd_expression = sd(Fold_change))

immunoblot_dko_sum$shRNA <- grepl("A6", immunoblot_dko_sum$Sample)

ppkr_ko_bt549 <-subset(immunoblot_dko, Protein == "pPKR/PKR" & Cell_line == "BT-549")

bt549 <- aov(Fold_change ~ Sample, data = ppkr_ko_bt549)
summary(bt549)

bt549ppkr <- TukeyHSD(bt549)
bt549ppkr <- as.data.frame(bt549ppkr$Sample[,1:4])

out <- strsplit(as.character(row.names(bt549ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549ppkr$stars <- makeStars(bt549ppkr$`p adj`)

bt549ppkr$group1 <- out$V1
bt549ppkr$group2 <- out$V2
bt549ppkr$y <- max(ppkr_ko_bt549$Fold_change*1.05)
bt549ppkr$Cell_line = "BT-549"



ppkr_ko_mb453 <- subset(immunoblot_dko, Protein == "pPKR/PKR" & Cell_line == "MDA-MB-453")

mb453 <- aov(Fold_change ~ Sample, data = ppkr_ko_mb453)
summary(mb453)

mb453ppkr <- TukeyHSD(mb453)
mb453ppkr <- as.data.frame(mb453ppkr$Sample[,1:4])

out <- strsplit(as.character(row.names(mb453ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb453ppkr$stars <- makeStars(mb453ppkr$`p adj`)

mb453ppkr$group1 <- out$V1
mb453ppkr$group2 <- out$V2
mb453ppkr$y <- max(ppkr_ko_mb453$Fold_change*1.05)
mb453ppkr$Cell_line = "MDA-MB-453"


tnbcppkr <- rbind(bt549ppkr, mb453ppkr)
tnbcppkr$y <- tnbcppkr$y
tnbcppkr <- subset(tnbcppkr, !stars == "ns")
tnbcppkr$Cell_line <- factor(tnbcppkr$Cell_line, levels = c("BT-549", "MDA-MB-453"))


p <- ggplot(subset(immunoblot_dko, Protein == "pPKR/PKR"), aes(Sample, Fold_change, colour = shRNA)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.3) + theme_science() + 
  scale_x_discrete(limits = c("C1_SCR", "NTA_SCR", "PACT1_SCR", "PACT2_SCR", "C1_A6", "NTA_A6", "PACT1_A6", "PACT2_A6"), 
                   labels = c("sgC1", "sgNTA", "sgPACT1", "sgPACT2", "sgC1", "sgNTA", "sgPACT1", "sgPACT2")) +
  labs(x = "", y = "Fold Change P-PKR/PKR", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 0, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + facet_grid2(cols = vars(Cell_line), scales = "free_y", independent = "y") +
  geom_col(data = subset(immunoblot_dko_sum, Protein == "pPKR/PKR"), aes(Sample, mean_expression, colour = shRNA, fill = shRNA), width = 0.5, linewidth = 0.3, alpha = 0.2) +
  geom_errorbar(data = subset(immunoblot_dko_sum, Protein == "pPKR/PKR"), aes(x = Sample, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression, colour = shRNA), width=0.2, linewidth = 0.3) + 
  stat_pvalue_manual(data = tnbcppkr, step.increase = 0.05, step.group.by = "Cell_line", label = "stars", y.position = "y", label.size = 3, linewidth = 0.3) +
  scale_colour_manual(values = palette2, labels = c("shSCR", "shADAR1")) + 
  scale_fill_manual(values = palette2, labels = c("shSCR", "shADAR1"))
move_layers(p, "GeomPoint", position = "top")
ggsave("/Manuscripts/PACT/Western_plots/ppkr_pact_adar_dko_facet_all.tiff", height = 2.1, width = 2.6, units = "in")

#remove comparisons between the two control sgRNAs or the two PACT sgRNAs to make the plot less cluttered 

tnbcppkr <- subset(tnbcppkr, !grepl("PACT1.*PACT2|PACT2.*PACT1", paste(tnbcppkr$group1, tnbcppkr$group2)))
tnbcppkr <- subset(tnbcppkr, !grepl("C1.*NTA|NTA.*C1", paste(tnbcppkr$group1, tnbcppkr$group2)))

p <- ggplot(subset(immunoblot_dko, Protein == "pPKR/PKR"), aes(Sample, Fold_change, colour = shRNA)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.3) + theme_science() + 
  scale_x_discrete(limits = c("C1_SCR", "NTA_SCR", "PACT1_SCR", "PACT2_SCR", "C1_A6", "NTA_A6", "PACT1_A6", "PACT2_A6"), 
                   labels = c("sgC1", "sgNTA", "sgPACT1", "sgPACT2", "sgC1", "sgNTA", "sgPACT1", "sgPACT2")) +
  labs(x = "", y = "Fold Change P-PKR/PKR", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 0, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + facet_grid2(cols = vars(Cell_line), scales = "free_y", independent = "y") +
  geom_col(data = subset(immunoblot_dko_sum, Protein == "pPKR/PKR"), aes(Sample, mean_expression, colour = shRNA, fill = shRNA), width = 0.5, linewidth = 0.3, alpha = 0.2) +
  geom_errorbar(data = subset(immunoblot_dko_sum, Protein == "pPKR/PKR"), aes(x = Sample, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression, colour = shRNA), width=0.2, linewidth = 0.3) + 
  stat_pvalue_manual(data = tnbcppkr, step.increase = 0.05, step.group.by = "Cell_line", label = "stars", y.position = "y", label.size = 3, linewidth = 0.3) +
  scale_colour_manual(values = palette2, labels = c("shSCR", "shADAR1")) + 
  scale_fill_manual(values = palette2, labels = c("shSCR", "shADAR1"))
move_layers(p, "GeomPoint", position = "top")
ggsave("/Manuscripts/PACT/Western_plots/ppkr_pact_adar_dko_facet.tiff", height = 2.1, width = 2.6, units = "in")


####Immunoblot PACT overexpression ----

immunoblot_oe <- rbind(western_pact_oe_bt549, western_pact_oe_mb453, western_pact_oe_hcc1806, western_pact_oe_mb468)


###Immunoblot TNBC lines
immunoblot_oe$Cell_line <- gsub("BT549", "BT-549", immunoblot_oe$Cell_line)
immunoblot_oe$Cell_line <- gsub("MB453", "MDA-MB-453", immunoblot_oe$Cell_line)
immunoblot_oe$Cell_line <- gsub("MB468", "MDA-MB-468", immunoblot_oe$Cell_line)
immunoblot_oe$Cell_line <- factor(immunoblot_oe$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))


#group and summarise data
immunoblot_oe_sum <- dplyr::group_by(immunoblot_oe, Cell_line, Protein, Sample)

immunoblot_oe_sum<- dplyr::summarise(immunoblot_oe_sum, mean_expression = mean(Fold_change), sd_expression = sd(Fold_change))

#stats as above for knockout of PACT
ppkr_oe_hcc1806 <-subset(immunoblot_oe, Protein == "pPKR" & Cell_line == "HCC1806")

hcc1806 <- aov(Fold_change ~ Sample, data = ppkr_oe_hcc1806)
summary(hcc1806)

hcc1806ppkr <- TukeyHSD(hcc1806)
hcc1806ppkr <- as.data.frame(hcc1806ppkr$Sample[,1:4])

out <- strsplit(as.character(row.names(hcc1806ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

hcc1806ppkr$stars <- makeStars(hcc1806ppkr$`p adj`)

hcc1806ppkr$group1 <- out$V1
hcc1806ppkr$group2 <- out$V2
hcc1806ppkr$y <- max(ppkr_oe_hcc1806$Fold_change*1.1)
hcc1806ppkr$Cell_line <- "HCC1806"



ppkr_oe_bt549 <-subset(immunoblot_oe, Protein == "pPKR" & Cell_line == "BT-549")

bt549 <- aov(Fold_change ~ Sample, data = ppkr_oe_bt549)
summary(bt549)

bt549ppkr <- TukeyHSD(bt549)
bt549ppkr <- as.data.frame(bt549ppkr$Sample[,1:4])

out <- strsplit(as.character(row.names(bt549ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

bt549ppkr$stars <- makeStars(bt549ppkr$`p adj`)

bt549ppkr$group1 <- out$V1
bt549ppkr$group2 <- out$V2
bt549ppkr$y <- max(ppkr_oe_bt549$Fold_change*1.1)
bt549ppkr$Cell_line = "BT-549"



ppkr_oe_mb453 <- subset(immunoblot_oe, Protein == "pPKR" & Cell_line == "MDA-MB-453")

mb453 <- aov(Fold_change ~ Sample, data = ppkr_oe_mb453)
summary(mb453)

mb453ppkr <- TukeyHSD(mb453)
mb453ppkr <- as.data.frame(mb453ppkr$Sample[,1:4])

out <- strsplit(as.character(row.names(mb453ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb453ppkr$stars <- makeStars(mb453ppkr$`p adj`)

mb453ppkr$group1 <- out$V1
mb453ppkr$group2 <- out$V2
mb453ppkr$y <- max(ppkr_oe_mb453$Fold_change*1.1)
mb453ppkr$Cell_line = "MDA-MB-453"


ppkr_oe_mb468 <- subset(immunoblot_oe, Protein == "pPKR" & Cell_line == "MDA-MB-468")

mb468 <- aov(Fold_change ~ Sample, data = ppkr_oe_mb468)
summary(mb468)

mb468ppkr <- TukeyHSD(mb468)
mb468ppkr <- as.data.frame(mb468ppkr$Sample[,1:4])

out <- strsplit(as.character(row.names(mb468ppkr)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

mb468ppkr$stars <- makeStars(mb468ppkr$`p adj`)

mb468ppkr$group1 <- out$V1
mb468ppkr$group2 <- out$V2
mb468ppkr$y <- max(ppkr_oe_mb468$Fold_change*1.1)
mb468ppkr$Cell_line = "MDA-MB-468"


tnbcppkr <- rbind(mb468ppkr, hcc1806ppkr, bt549ppkr, mb453ppkr)
tnbcppkr$y <- tnbcppkr$y
tnbcppkr <- subset(tnbcppkr, !stars == "ns")
tnbcppkr$Cell_line <- factor(tnbcppkr$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))

strip <- strip_themed(text_x = elem_list_text(colour = pallete4))

ggplot(subset(immunoblot_oe, Protein == "pPKR"), aes(Sample, Fold_change, colour = Cell_line)) +
  geom_jitter(aes(shape = as.factor(Replicate)), width = 0.12, size = 0.3) + theme_science() + 
  scale_x_discrete(limits = c("EV", "PACT", "PACT S287A", "PACT S287D")) + labs(x = "", y = "Fold Change P-PKR/PKR", colour = "", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 0, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none", colour = "none", fill = "none") + facet_grid2(cols = vars(Cell_line), strip = strip) +
  scale_colour_manual(values = pallete4) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  geom_col(data = subset(immunoblot_oe_sum, Protein == "pPKR"), aes(Sample, mean_expression, fill = Cell_line), width = 0.5, alpha = 0.2, linewidth = 0.3) +
  geom_errorbar(data = subset(immunoblot_oe_sum, Protein == "pPKR"), aes(x = Sample, ymax = mean_expression + sd_expression, ymin = mean_expression - sd_expression, y = mean_expression), width=0.2, linewidth = 0.3) + 
  scale_fill_manual(values = pallete4)
ggsave("/Manuscripts/PACT/Western_plots/ppkr_oe_facet.tiff", height =2, width = 2.9, units = "in", dpi = 300)




####qPCR PACT ko ----
qPCR_ko <- na.omit(qpcr_relative_pact_ko)
#convert to linear scale
qPCR_ko$Fold_Change <- 2^-(qPCR_ko$Fold_Change)
qPCR_ko$Cell_line <- gsub("549", "BT-549", qPCR_ko$Cell_line)
qPCR_ko$Cell_line <- gsub("453", "MDA-MB-453", qPCR_ko$Cell_line)
qPCR_ko$Cell_line <- gsub("468", "MDA-MB-468", qPCR_ko$Cell_line)
qPCR_ko$Cell_line <- gsub("1806", "HCC1806", qPCR_ko$Cell_line)
qPCR_ko$Cell_line <- factor(qPCR_ko$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))


#group and summarise data
qPCR_ko_sum <- dplyr::group_by(qPCR_ko, Cell_line, Gene, Sample)

qPCR_ko_sum<- dplyr::summarise(qPCR_ko_sum, mean_expression = mean(Fold_Change), sd_expression = sd(Fold_Change))

qPCR_ko$condition <- ifelse(grepl("PACT", qPCR_ko$Sample), "WT", "KO")

#stats as above for knockout of PACT, grouped by gene set

atf4_1806 <- aov(Fold_Change ~ Sample, data = subset(qPCR_ko, Gene %in% c("ATF3", "ASNS", "SESN2") & Cell_line == "HCC1806"))
summary(atf4_1806)

atf4_1806 <- TukeyHSD(atf4_1806)
atf4_1806 <- as.data.frame(atf4_1806$Sample[,1:4])

out <- strsplit(as.character(row.names(atf4_1806)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

atf4_1806$stars <- makeStars(atf4_1806$`p adj`)

atf4_1806$group1 <- out$V1
atf4_1806$group2 <- out$V2
atf4_1806$y <- max(qPCR_ko$Fold_Change[qPCR_ko$Cell_line == "HCC1806" & 
                                        qPCR_ko$Gene %in% c("ATF3", "ASNS", "SESN2") ]*1.1)
atf4_1806$Cell_line <- "HCC1806"

atf4_468 <- aov(Fold_Change ~ Sample, data = subset(qPCR_ko, Gene %in% c("ATF3", "ASNS", "SESN2") & Cell_line == "MDA-MB-468"))
summary(atf4_468)

atf4_468 <- TukeyHSD(atf4_468)
atf4_468 <- as.data.frame(atf4_468$Sample[,1:4])

out <- strsplit(as.character(row.names(atf4_468)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

atf4_468$stars <- makeStars(atf4_468$`p adj`)

atf4_468$group1 <- out$V1
atf4_468$group2 <- out$V2
atf4_468$y <- max(qPCR_ko$Fold_Change[qPCR_ko$Cell_line == "MDA-MB-468" & 
                                       qPCR_ko$Gene %in% c("ATF3", "ASNS", "SESN2") ]*1.1)
atf4_468$Cell_line <- "MDA-MB-468"

atf4_549 <- aov(Fold_Change ~ Sample, data = subset(qPCR_ko, Gene %in% c("ATF3", "ASNS", "SESN2") & Cell_line == "BT-549"))
summary(atf4_549)

atf4_549 <- TukeyHSD(atf4_549)
atf4_549 <- as.data.frame(atf4_549$Sample[,1:4])

out <- strsplit(as.character(row.names(atf4_549)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

atf4_549$stars <- makeStars(atf4_549$`p adj`)

atf4_549$group1 <- out$V1
atf4_549$group2 <- out$V2
atf4_549$y <- max(qPCR_ko$Fold_Change[qPCR_ko$Cell_line == "BT-549" & 
                                       qPCR_ko$Gene %in% c("ATF3", "ASNS", "SESN2") ]*1.1)
atf4_549$Cell_line <- "BT-549"

atf4_453 <- aov(Fold_Change ~ Sample, data = subset(qPCR_ko, Gene %in% c("ATF3", "ASNS", "SESN2") & Cell_line == "MDA-MB-453"))
summary(atf4_453)

atf4_453 <- TukeyHSD(atf4_453)
atf4_453 <- as.data.frame(atf4_453$Sample[,1:4])

out <- strsplit(as.character(row.names(atf4_453)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

atf4_453$stars <- makeStars(atf4_453$`p adj`)

atf4_453$group1 <- out$V1
atf4_453$group2 <- out$V2
atf4_453$y <- max(qPCR_ko$Fold_Change[qPCR_ko$Cell_line == "MDA-MB-453" & 
                                       qPCR_ko$Gene %in% c("ATF3", "ASNS", "SESN2") ]*1.1)
atf4_453$Cell_line <- "MDA-MB-453"


atf4_ko <- rbind(atf4_468, atf4_1806, atf4_549, atf4_453)
atf4_ko$y <- atf4_ko$y
atf4_ko <- subset(atf4_ko, !stars == "ns")
atf4_ko$Cell_line <- factor(atf4_ko$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))



#make data.frame for annotation of plot to add brackets above samples grouping the three genes

bracket_468 <- data.frame(x =  c(0.6,1.6,2.6,NA), x_end = c(1.4,2.4,3.4,NA), y= c(1.4, 1.4, 5.1, NA),
                       Cell_line = "MDA-MB-468")
bracket_1806 <- data.frame(x =  c(0.6,1.6,2.6,3.6), x_end = c(1.4,2.4,3.4,4.4), y= c(1.4, 2.4, 10, 8.2),
                          Cell_line = "HCC1806")
bracket_549 <- data.frame(x =  NA, x_end = NA, y= NA,
                          Cell_line = "BT-549")
bracket_453 <- data.frame(x =  NA, x_end = NA, y= NA,
                          Cell_line = "MDA-MB-453")
brackets <- rbind(bracket_1806, bracket_468, bracket_549, bracket_453)
                       
brackets$Cell_line <- factor(brackets$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))
brackets$Gene = NA


ggplot(subset(qPCR_ko, Gene %in% c("ATF3", "ASNS", "SESN2")), aes(Sample, Fold_Change, colour = Gene, group = Gene)) +
  geom_point(aes(shape = as.factor(Replicate)), size = 0.3, 
             position = position_jitterdodge(dodge.width = 0.7)) + theme_science() + 
  scale_x_discrete(limits = c("sgC1", "sgNTA", "sgPACT.1", "sgPACT.2"), labels = c("sgC1", "sgNTA", "sgPACT-1", "sgPACT-2")) + 
  labs(x = "", y = "Relative Expression", colour = "", fill = "", title = "ATF4 Target Genes") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = c(0.1,0.8), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + facet_grid2(cols = vars(Cell_line)) +
  scale_colour_manual(values = pallete4) + scale_fill_manual(values = pallete4) + 
  geom_col(data = subset(qPCR_ko_sum, Gene %in% c("ATF3", "ASNS", "SESN2")), aes(Sample, mean_expression, colour = Gene, fill = Gene),
           width = 0.5, alpha = 0.2, linewidth = 0.3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = subset(qPCR_ko_sum, Gene %in% c("ATF3", "ASNS", "SESN2")), 
                aes(x = Sample, colour = Gene, ymax = mean_expression + sd_expression, 
                    ymin = mean_expression - sd_expression, y = mean_expression), width=0.2, linewidth = 0.3,
                position = position_dodge(width = 0.7)) +
  stat_pvalue_manual(data = atf4_ko, label = "stars", y.position = "y", step.increase = 0.07, step.group.by = "Cell_line", label.size = 2.5, linewidth = 0.3, tip.length = 0.01) +
  geom_bracket(data = brackets, aes(xmin = x, xmax = x_end, y.position = y, label = ""), colour = "black") 
ggsave("/Manuscripts/PACT/qPCR_plots/pact_ko_atf4.tiff", height =2.2, width = 3.3, units = "in", dpi = 300)

#remove bracket annotation data.frames

bracket_468 <- NULL
bracket_1806 <- NULL
bracket_549 <- NULL
bracket_453 <- NULL
brackets <- NULL


isg_1806 <- aov(Fold_Change ~ Sample, data = subset(qPCR_ko, Gene %in% c("PARP14", "IFIT2", "IFIH1") & Cell_line == "HCC1806"))
summary(isg_1806)

isg_1806 <- TukeyHSD(isg_1806)
isg_1806 <- as.data.frame(isg_1806$Sample[,1:4])

out <- strsplit(as.character(row.names(isg_1806)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

isg_1806$stars <- makeStars(isg_1806$`p adj`)

isg_1806$group1 <- out$V1
isg_1806$group2 <- out$V2
isg_1806$y <- max(qPCR_ko$Fold_Change[qPCR_ko$Cell_line == "HCC1806" & 
                                         qPCR_ko$Gene %in% c("PARP14", "IFIT2", "IFIH1") ]*1.1)
isg_1806$Cell_line <- "HCC1806"

isg_468 <- aov(Fold_Change ~ Sample, data = subset(qPCR_ko, Gene %in% c("PARP14", "IFIT2", "IFIH1") & Cell_line == "MDA-MB-468"))
summary(isg_468)

isg_468 <- TukeyHSD(isg_468)
isg_468 <- as.data.frame(isg_468$Sample[,1:4])

out <- strsplit(as.character(row.names(isg_468)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

isg_468$stars <- makeStars(isg_468$`p adj`)

isg_468$group1 <- out$V1
isg_468$group2 <- out$V2
isg_468$y <- max(qPCR_ko$Fold_Change[qPCR_ko$Cell_line == "MDA-MB-468" & 
                                        qPCR_ko$Gene %in% c("PARP14", "IFIT2", "IFIH1") ]*1.1)
isg_468$Cell_line <- "MDA-MB-468"

isg_549 <- aov(Fold_Change ~ Sample, data = subset(qPCR_ko, Gene %in% c("PARP14", "IFIT2", "IFIH1") & Cell_line == "BT-549"))
summary(isg_549)

isg_549 <- TukeyHSD(isg_549)
isg_549 <- as.data.frame(isg_549$Sample[,1:4])

out <- strsplit(as.character(row.names(isg_549)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

isg_549$stars <- makeStars(isg_549$`p adj`)

isg_549$group1 <- out$V1
isg_549$group2 <- out$V2
isg_549$y <- max(qPCR_ko$Fold_Change[qPCR_ko$Cell_line == "BT-549" & 
                                        qPCR_ko$Gene %in% c("PARP14", "IFIT2", "IFIH1") ]*1.1)
isg_549$Cell_line <- "BT-549"

isg_453 <- aov(Fold_Change ~ Sample, data = subset(qPCR_ko, Gene %in% c("PARP14", "IFIT2", "IFIH1") & Cell_line == "MDA-MB-453"))
summary(isg_453)

isg_453 <- TukeyHSD(isg_453)
isg_453 <- as.data.frame(isg_453$Sample[,1:4])

out <- strsplit(as.character(row.names(isg_453)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

isg_453$stars <- makeStars(isg_453$`p adj`)

isg_453$group1 <- out$V1
isg_453$group2 <- out$V2
isg_453$y <- max(qPCR_ko$Fold_Change[qPCR_ko$Cell_line == "MDA-MB-453" & 
                                        qPCR_ko$Gene %in% c("PARP14", "IFIT2", "IFIH1") ]*1.1)
isg_453$Cell_line <- "MDA-MB-453"


isg_ko <- rbind(isg_468, isg_1806, isg_549, isg_453)
isg_ko$y <- isg_ko$y
isg_ko <- subset(isg_ko, !stars == "ns")
isg_ko$Cell_line <- factor(isg_ko$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))


bracket_468 <- data.frame(x =  NA, x_end = NA, y= NA,
                          Cell_line = "MDA-MB-468")
bracket_1806 <- data.frame(x =  c(0.6,1.6,2.6,3.6), x_end = c(1.4,2.4,3.4,4.4), y= c(1.5, 1.5, 2.4, 2.5),
                           Cell_line = "HCC1806")
bracket_549 <- data.frame(x =  NA, x_end = NA, y= NA,
                          Cell_line = "BT-549")
bracket_453 <- data.frame(x =  c(0.6,1.6,2.6,NA), x_end = c(1.4,2.4,3.4,NA), y= c(1.8, 1.35, 4.65, NA),
                          Cell_line = "MDA-MB-453")
brackets <- rbind(bracket_1806, bracket_468, bracket_549, bracket_453)

brackets$Cell_line <- factor(brackets$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))
brackets$Gene = NA


ggplot(subset(qPCR_ko, Gene %in% c("PARP14", "IFIT2", "IFIH1")), aes(Sample, Fold_Change, colour = Gene, group = Gene)) +
  geom_point(aes(shape = as.factor(Replicate)), size = 0.3, 
             position = position_jitterdodge(dodge.width = 0.7)) + theme_science() + 
  scale_x_discrete(limits = c("sgC1", "sgNTA", "sgPACT.1", "sgPACT.2"), labels = c("sgC1", "sgNTA", "sgPACT-1", "sgPACT-2")) + 
  labs(x = "", y = "Relative Expression", colour = "", fill = "", title  =  "IFN Stimulated Genes") + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = c(0.1,0.8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + facet_grid2(cols = vars(Cell_line)) +
  scale_colour_manual(values = pallete4) + scale_fill_manual(values = pallete4) +
  geom_col(data = subset(qPCR_ko_sum, Gene %in% c("PARP14", "IFIT2", "IFIH1")), aes(Sample, mean_expression, colour = Gene, fill = Gene),
           width = 0.5, alpha = 0.2, linewidth = 0.3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = subset(qPCR_ko_sum, Gene %in% c("PARP14", "IFIT2", "IFIH1")), 
                aes(x = Sample, colour = Gene, ymax = mean_expression + sd_expression, 
                    ymin = mean_expression - sd_expression, y = mean_expression), width=0.2, linewidth = 0.3,
                position = position_dodge(width = 0.7)) +
  stat_pvalue_manual(data = isg_ko, label = "stars", step.increase = 0.07, , step.group.by = "Cell_line", y.position = "y", label.size = 2.5, linewidth = 0.3, tip.length = 0.01) +
  geom_bracket(data = brackets, aes(xmin = x, xmax = x_end, y.position = y, label = ""), colour = "black") 
ggsave("/Manuscripts/PACT/qPCR_plots/pact_ko_isg.tiff", height =2.2, width = 3.3, units = "in", dpi = 300)


bracket_468 <- NULL
bracket_1806 <- NULL
bracket_549 <- NULL
bracket_453 <- NULL
brackets <- NULL

nfkb_1806 <- aov(Fold_Change ~ Sample, data = subset(qPCR_ko, Gene %in% c("NFKBIE", "IRF1", "BIRC3") & Cell_line == "HCC1806"))
summary(nfkb_1806)

nfkb_1806 <- TukeyHSD(nfkb_1806)
nfkb_1806 <- as.data.frame(nfkb_1806$Sample[,1:4])

out <- strsplit(as.character(row.names(nfkb_1806)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

nfkb_1806$stars <- makeStars(nfkb_1806$`p adj`)

nfkb_1806$group1 <- out$V1
nfkb_1806$group2 <- out$V2
nfkb_1806$y <- max(qPCR_ko$Fold_Change[qPCR_ko$Cell_line == "HCC1806" & 
                                         qPCR_ko$Gene %in% c("NFKBIE", "IRF1", "BIRC3") ]*1.1)
nfkb_1806$Cell_line <- "HCC1806"

nfkb_468 <- aov(Fold_Change ~ Sample, data = subset(qPCR_ko, Gene %in% c("NFKBIE", "IRF1", "BIRC3") & Cell_line == "MDA-MB-468"))
summary(nfkb_468)

nfkb_468 <- TukeyHSD(nfkb_468)
nfkb_468 <- as.data.frame(nfkb_468$Sample[,1:4])

out <- strsplit(as.character(row.names(nfkb_468)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

nfkb_468$stars <- makeStars(nfkb_468$`p adj`)

nfkb_468$group1 <- out$V1
nfkb_468$group2 <- out$V2
nfkb_468$y <- max(qPCR_ko$Fold_Change[qPCR_ko$Cell_line == "MDA-MB-468" & 
                                        qPCR_ko$Gene %in% c("NFKBIE", "IRF1", "BIRC3") ]*1.3)
nfkb_468$Cell_line <- "MDA-MB-468"

nfkb_549 <- aov(Fold_Change ~ Sample, data = subset(qPCR_ko, Gene %in% c("NFKBIE", "IRF1") & Cell_line == "BT-549"))
summary(nfkb_549)

nfkb_549 <- TukeyHSD(nfkb_549)
nfkb_549 <- as.data.frame(nfkb_549$Sample[,1:4])

out <- strsplit(as.character(row.names(nfkb_549)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

nfkb_549$stars <- makeStars(nfkb_549$`p adj`)

nfkb_549$group1 <- out$V1
nfkb_549$group2 <- out$V2
nfkb_549$y <- max(qPCR_ko$Fold_Change[qPCR_ko$Cell_line == "BT-549" & 
                                        qPCR_ko$Gene %in% c("NFKBIE", "IRF1", "BIRC3") ]*1.1)
nfkb_549$Cell_line <- "BT-549"

nfkb_453 <- aov(Fold_Change ~ Sample, data = subset(qPCR_ko, Gene %in% c("NFKBIE", "IRF1", "BIRC3") & Cell_line == "MDA-MB-453"))
summary(nfkb_453)

nfkb_453 <- TukeyHSD(nfkb_453)
nfkb_453 <- as.data.frame(nfkb_453$Sample[,1:4])

out <- strsplit(as.character(row.names(nfkb_453)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

nfkb_453$stars <- makeStars(nfkb_453$`p adj`)

nfkb_453$group1 <- out$V1
nfkb_453$group2 <- out$V2
nfkb_453$y <- max(qPCR_ko$Fold_Change[qPCR_ko$Cell_line == "MDA-MB-453" & 
                                        qPCR_ko$Gene %in% c("NFKBIE", "IRF1", "BIRC3") ]*1.1)
nfkb_453$Cell_line <- "MDA-MB-453"


nfkb_ko <- rbind(nfkb_468, nfkb_1806, nfkb_549, nfkb_453)
nfkb_ko$y <- nfkb_ko$y
nfkb_ko <- subset(nfkb_ko, !stars == "ns")
nfkb_ko$Cell_line <- factor(nfkb_ko$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))


bracket_468 <- data.frame(x =  c(0.6,1.6,2.6,3.6), x_end = c(1.4,2.4,3.4,4.4), y= c(1.4, 1.4, 2.5, 2.4),
                          Cell_line = "MDA-MB-468")
bracket_1806 <- data.frame(x =  c(0.6,1.6,2.6,3.6), x_end = c(1.4,2.4,3.4,4.4), y= c(1.6, 1.6, 8.8, 8.2),
                           Cell_line = "HCC1806")
bracket_549 <- data.frame(x =  NA, x_end = NA, y= NA,
                          Cell_line = "BT-549")
bracket_453 <- data.frame(x =  NA, x_end = NA, y= NA,
                          Cell_line = "MDA-MB-453")
brackets <- rbind(bracket_1806, bracket_468, bracket_549, bracket_453)

brackets$Cell_line <- factor(brackets$Cell_line, levels = c("BT-549", "MDA-MB-453", "MDA-MB-468", "HCC1806"))
brackets$Gene = NA


ggplot(subset(qPCR_ko, Gene %in% c("NFKBIE", "IRF1", "BIRC3")), aes(Sample, Fold_Change, colour = Gene, group = Gene)) +
  geom_point(aes(shape = as.factor(Replicate)), size = 0.3, 
             position = position_jitterdodge(dodge.width = 0.7)) + theme_science() + 
  scale_x_discrete(limits = c("sgC1", "sgNTA", "sgPACT.1", "sgPACT.2"), labels = c("sgC1", "sgNTA", "sgPACT-1", "sgPACT-2")) + 
  labs(x = "", y = "Relative Expression", colour = "", fill = "", title = "NF-kB Target Genes") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = c(0.1,0.8), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + facet_grid2(cols = vars(Cell_line)) +
  scale_colour_manual(values = pallete4) + scale_fill_manual(values = pallete4) +
  geom_col(data = subset(qPCR_ko_sum, Gene %in% c("NFKBIE", "IRF1", "BIRC3")), aes(Sample, mean_expression, colour = Gene, fill = Gene),
           width = 0.5, alpha = 0.2, linewidth = 0.3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = subset(qPCR_ko_sum, Gene %in% c("NFKBIE", "IRF1", "BIRC3")), 
                aes(x = Sample, colour = Gene, ymax = mean_expression + sd_expression, 
                    ymin = mean_expression - sd_expression, y = mean_expression), width=0.2, linewidth = 0.3,
                position = position_dodge(width = 0.7)) +
  stat_pvalue_manual(data = nfkb_ko, label = "stars", step.increase = 0.07, step.group.by = "Cell_line", y.position = "y", label.size = 2.5, linewidth = 0.3, tip.length = 0.01) +
  geom_bracket(data = brackets, aes(xmin = x, xmax = x_end, y.position = y, label = ""), colour = "black") 
ggsave("/Manuscripts/PACT/qPCR_plots/pact_ko_nfkb.tiff", height =2.2, width = 3.3, units = "in", dpi = 300)

bracket_468 <- NULL
bracket_1806 <- NULL
bracket_549 <- NULL
bracket_453 <- NULL
brackets <- NULL


####qPCR PACT/ADAR1 dko ----
qPCR_dko_pact_adar$Fold_Change <- 2^-(qPCR_dko_pact_adar$Fold_Change)

qPCR_dko_pact_adar$shRNA <- ifelse(grepl("ADAR", qPCR_dko_pact_adar$Sample), "shADAR", "shSCR")
qPCR_dko_pact_adar$shRNA <- factor(qPCR_dko_pact_adar$shRNA, levels = c("shSCR", "shADAR"))


#group and summarise data
qPCR_dko_pact_adar_sum <- dplyr::group_by(qPCR_dko_pact_adar, Gene, Sample)

qPCR_dko_pact_adar_sum<- dplyr::summarise(qPCR_dko_pact_adar_sum, mean_expression = mean(Fold_Change), sd_expression = sd(Fold_Change))

qPCR_dko_pact_adar_sum$shRNA <- ifelse(grepl("ADAR", qPCR_dko_pact_adar_sum$Sample), "shADAR", "shSCR")
qPCR_dko_pact_adar_sum$shRNA <- factor(qPCR_dko_pact_adar_sum$shRNA, levels = c("shSCR", "shADAR"))

qPCR_dko_pact_adar$sgRNA <- gsub("_shADAR", "", qPCR_dko_pact_adar$Sample)
qPCR_dko_pact_adar$sgRNA <- gsub("_shSCR", "", qPCR_dko_pact_adar$sgRNA)

qPCR_dko_pact_adar_sum$sgRNA <- gsub("_shADAR", "", qPCR_dko_pact_adar_sum$Sample)
qPCR_dko_pact_adar_sum$sgRNA <- gsub("_shSCR", "", qPCR_dko_pact_adar_sum$sgRNA)

#stats analysis as in above section
atf4_453 <- aov(Fold_Change ~ Sample, , data = subset(qPCR_dko_pact_adar, Gene %in% c("ATF3", "ASNS", "SESN2")))
summary(atf4_453)

atf4_453_t <- TukeyHSD(atf4_453)

atf4_453_t  <- as.data.frame(atf4_453_t $Sample)

atf4_453_t <- subset(atf4_453_t, `p adj` < 0.05)

atf4_453_t$stars <- makeStars(atf4_453_t$`p adj`)

out <- strsplit(as.character(row.names(atf4_453_t)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

atf4_453_t$group1 <- out$V1
atf4_453_t$group2 <- out$V2

atf4_453_t$y <- max(qPCR_dko_pact_adar$Fold_Change[qPCR_dko_pact_adar$Gene %in% c("ATF3", "ASNS", "SESN2")]*1.1)

atf4_453_t <- atf4_453_t[!(grepl("PACT.1", row.names(atf4_453_t)) & grepl("PACT.2", row.names(atf4_453_t))),]


                    

ggplot(subset(qPCR_dko_pact_adar, Gene %in% c("ATF3", "ASNS", "SESN2")), aes(sgRNA, Fold_Change, colour = Gene, group = Gene)) +
  geom_point(aes(shape = as.factor(Replicate)), size = 0.3, 
             position = position_jitterdodge(dodge.width = 0.7)) + theme_science() + 
  scale_x_discrete(limits = c("sgC1", "sgNTA", "sgPACT.1", "sgPACT.2"), labels = c("sgC1", "sgNTA", "sgPACT-1", "sgPACT-2")) + 
  labs(x = "", y = "Relative Expression", colour = "", fill = "", title = "ATF4 Target Genes (MDA-MB-453)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = c(0.2, 0.8), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + facet_grid2(cols = vars(shRNA)) +
  scale_colour_manual(values = pallete4) + scale_fill_manual(values = pallete4) +
  geom_col(data = subset(qPCR_dko_pact_adar_sum, Gene %in% c("ATF3", "ASNS", "SESN2")), aes(sgRNA, mean_expression, colour = Gene, fill = Gene),
           width = 0.5, alpha = 0.2, linewidth = 0.3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = subset(qPCR_dko_pact_adar_sum, Gene %in% c("ATF3", "ASNS", "SESN2")), 
                aes(x = sgRNA, colour = Gene, ymax = mean_expression + sd_expression, 
                    ymin = mean_expression - sd_expression, y = mean_expression), width=0.2, linewidth = 0.3,
                position = position_dodge(width = 0.7)) 
ggsave("/Manuscripts/PACT/qPCR_plots/pact_adar_dko_453_atf4.tiff", height =1.8, width = 2.6, units = "in", dpi = 300)


brackets <- data.frame(x =  c(0.6,1.6,2.6,3.6,4.6, 5.6, 6.6, 7.6), 
                       x_end = c(1.4,2.4,3.4,4.4, 5.4, 6.4, 7.4, 8.4), 
                       y= c(2.7,2.7,2.7,2.7,2.7,2.7, 40, 35))


brackets$Gene = NA

ggplot(subset(qPCR_dko_pact_adar, Gene %in% c("ATF3", "ASNS", "SESN2")), aes(Sample, Fold_Change, colour = Gene, group = Gene)) +
  geom_point(aes(shape = as.factor(Replicate)), size = 0.3, 
             position = position_jitterdodge(dodge.width = 0.7)) + theme_science() + 
  scale_x_discrete(limits = c("sgC1_shSCR", "sgNTA_shSCR", "sgPACT.1_shSCR", "sgPACT.2_shSCR",
                              "sgC1_shADAR", "sgNTA_shADAR", "sgPACT.1_shADAR", "sgPACT.2_shADAR"), 
                   labels = c("sgC1", "sgNTA", "sgPACT-1", "sgPACT-2",
                              "sgC1", "sgNTA", "sgPACT-1", "sgPACT-2")) + 
  labs(x = "", y = "Relative Expression", colour = "", fill = "", title = "ATF4 Target Genes (MDA-MB-453)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = c(0.2, 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + 
  scale_colour_manual(values = pallete4) + scale_fill_manual(values = pallete4) +
  geom_col(data = subset(qPCR_dko_pact_adar_sum, Gene %in% c("ATF3", "ASNS", "SESN2")), aes(Sample, mean_expression, colour = Gene, fill = Gene),
           width = 0.5, alpha = 0.2, linewidth = 0.3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = subset(qPCR_dko_pact_adar_sum, Gene %in% c("ATF3", "ASNS", "SESN2")), 
                aes(x = Sample, colour = Gene, ymax = mean_expression + sd_expression, 
                    ymin = mean_expression - sd_expression, y = mean_expression), width=0.2, linewidth = 0.3,
                position = position_dodge(width = 0.7))  +
  stat_pvalue_manual(data = atf4_453_t, label = "stars", step.increase = 0.02, 
                     y.position = "y", label.size = 2.5, linewidth = 0.3, tip.length = 0.01) +
  annotate("segment", x = 1, xend = 4, y = 57, yend = 57, colour = "black") + 
  annotate("segment", x = 5, xend = 8, y = 57, yend = 57, colour = "black") +
  annotate("text", x = 2.5, y = 60, colour = "black", label = "shSCR", size = 3) + 
  annotate("text", x = 6.5, y = 60, colour = "black", label = "shADAR", size = 3) +
  geom_bracket(data = brackets, aes(xmin = x, xmax = x_end, y.position = y, label = ""), colour = "black") 
ggsave("/Manuscripts/PACT/qPCR_plots/pact_adar_dko_453_atf4_stats.tiff", height =2.2, width = 2.5, units = "in", dpi = 300)

isg_453 <- aov(Fold_Change ~ Sample, , data = subset(qPCR_dko_pact_adar, Gene %in% c("IFIT2", "IFIH1", "PARP14")))
summary(isg_453)

isg_453_t <- TukeyHSD(isg_453)

isg_453_t  <- as.data.frame(isg_453_t$Sample)

isg_453_t <- subset(isg_453_t, `p adj` < 0.05)

isg_453_t$stars <- makeStars(isg_453_t$`p adj`)

out <- strsplit(as.character(row.names(isg_453_t)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

isg_453_t$group1 <- out$V1
isg_453_t$group2 <- out$V2

isg_453_t$y <- max(qPCR_dko_pact_adar$Fold_Change[qPCR_dko_pact_adar$Gene %in% c("IFIT2", "IFIH1", "PARP14")]*1.1)

isg_453_t <- isg_453_t[!(grepl("PACT.1", row.names(isg_453_t)) & grepl("PACT.2", row.names(isg_453_t))),]


ggplot(subset(qPCR_dko_pact_adar, Gene %in% c("IFIT2", "IFIH1", "PARP14")), aes(sgRNA, Fold_Change, colour = Gene, group = Gene)) +
  geom_point(aes(shape = as.factor(Replicate)), size = 0.3, 
             position = position_jitterdodge(dodge.width = 0.7)) + theme_science() + 
  scale_x_discrete(limits = c("sgC1", "sgNTA", "sgPACT.1", "sgPACT.2"), labels = c("sgC1", "sgNTA", "sgPACT-1", "sgPACT-2")) + 
  labs(x = "", y = "Relative Expression", colour = "", fill = "", title = "ISGs (MDA-MB-453)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = c(0.2, 0.8), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + facet_grid2(cols = vars(shRNA)) +
  scale_colour_manual(values = pallete4) + scale_fill_manual(values = pallete4) +
  geom_col(data = subset(qPCR_dko_pact_adar_sum, Gene %in% c("IFIT2", "IFIH1", "PARP14")), aes(sgRNA, mean_expression, colour = Gene, fill = Gene),
           width = 0.5, alpha = 0.2, linewidth = 0.3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = subset(qPCR_dko_pact_adar_sum, Gene %in% c("IFIT2", "IFIH1", "PARP14")), 
                aes(x = sgRNA, colour = Gene, ymax = mean_expression + sd_expression, 
                    ymin = mean_expression - sd_expression, y = mean_expression), width=0.2, linewidth = 0.3,
                position = position_dodge(width = 0.7)) 
ggsave("/Manuscripts/PACT/qPCR_plots/pact_adar_dko_453_isg.tiff", height =1.8, width = 2.6, units = "in", dpi = 300)



brackets <- data.frame(x =  c(NA,NA,NA,NA,NA, 5.6, 6.6, NA), 
                       x_end = c(NA,NA,NA,NA,NA,  6.4, 7.4, NA), 
                       y= c(NA,NA,NA,NA,NA, 2, 5.8, NA))


brackets$Gene = NA


ggplot(subset(qPCR_dko_pact_adar, Gene %in% c("IFIT2", "IFIH1", "PARP14")), aes(Sample, Fold_Change, colour = Gene, group = Gene)) +
  geom_point(aes(shape = as.factor(Replicate)), size = 0.3, 
             position = position_jitterdodge(dodge.width = 0.7)) + theme_science() + 
  scale_x_discrete(limits = c("sgC1_shSCR", "sgNTA_shSCR", "sgPACT.1_shSCR", "sgPACT.2_shSCR",
                              "sgC1_shADAR", "sgNTA_shADAR", "sgPACT.1_shADAR", "sgPACT.2_shADAR"), 
                   labels = c("sgC1", "sgNTA", "sgPACT-1", "sgPACT-2",
                              "sgC1", "sgNTA", "sgPACT-1", "sgPACT-2")) + 
  labs(x = "", y = "Relative Expression", colour = "", fill = "", title = "ISGs (MDA-MB-453)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = c(0.2, 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + 
  scale_colour_manual(values = pallete4) + scale_fill_manual(values = pallete4) +
  geom_col(data = subset(qPCR_dko_pact_adar_sum, Gene %in% c("IFIT2", "IFIH1", "PARP14")), aes(Sample, mean_expression, colour = Gene, fill = Gene),
           width = 0.5, alpha = 0.2, linewidth = 0.3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = subset(qPCR_dko_pact_adar_sum, Gene %in% c("IFIT2", "IFIH1", "PARP14")), 
                aes(x = Sample, colour = Gene, ymax = mean_expression + sd_expression, 
                    ymin = mean_expression - sd_expression, y = mean_expression), width=0.2, linewidth = 0.3,
                position = position_dodge(width = 0.7))  +
  stat_pvalue_manual(data = isg_453_t, label = "stars", step.increase = 0.02, 
                     y.position = "y", label.size = 2.5, linewidth = 0.3, tip.length = 0.01) +
  annotate("segment", x = 1, xend = 4, y = 7.5, yend = 7.5, colour = "black") + 
  annotate("segment", x = 5, xend = 8, y = 7.5, yend = 7.5, colour = "black") +
  annotate("text", x = 2.5, y = 8, colour = "black", label = "shSCR", size = 3) + 
  annotate("text", x = 6.5, y = 8, colour = "black", label = "shADAR", size = 3) +
  geom_bracket(data = brackets, aes(xmin = x, xmax = x_end, y.position = y, label = ""), colour = "black")
ggsave("/Manuscripts/PACT/qPCR_plots/pact_adar_dko_453_isg_stats.tiff", height =2.2, width = 2.5, units = "in", dpi = 300)

nfkb_453 <- aov(Fold_Change ~ Sample, , data = subset(qPCR_dko_pact_adar, Gene %in% c("NFKBIE", "IRF1", "BIRC3")))
summary(nfkb_453)

nfkb_453_t <- TukeyHSD(nfkb_453)

nfkb_453_t  <- as.data.frame(nfkb_453_t$Sample)

nfkb_453_t <- subset(nfkb_453_t, `p adj` < 0.05)

nfkb_453_t$stars <- makeStars(nfkb_453_t$`p adj`)

out <- strsplit(as.character(row.names(nfkb_453_t)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

nfkb_453_t$group1 <- out$V1
nfkb_453_t$group2 <- out$V2

nfkb_453_t$y <- max(qPCR_dko_pact_adar$Fold_Change[qPCR_dko_pact_adar$Gene %in% c("NFKBIE", "IRF1", "BIRC3")]*1.1)

nfkb_453_t <- nfkb_453_t[!(grepl("PACT.1", row.names(nfkb_453_t)) & grepl("PACT.2", row.names(nfkb_453_t))),]


ggplot(subset(qPCR_dko_pact_adar, Gene %in% c("NFKBIE", "IRF1", "BIRC3")), aes(sgRNA, Fold_Change, colour = Gene, group = Gene)) +
  geom_point(aes(shape = as.factor(Replicate)), size = 0.3, 
             position = position_jitterdodge(dodge.width = 0.7)) + theme_science() + 
  scale_x_discrete(limits = c("sgC1", "sgNTA", "sgPACT.1", "sgPACT.2"), labels = c("sgC1", "sgNTA", "sgPACT-1", "sgPACT-2")) + 
  labs(x = "", y = "Relative Expression", colour = "", fill = "", title = "NF-kB Target Genes (MDA-MB-453)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = c(0.2, 0.9), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + facet_grid2(cols = vars(shRNA)) +
  scale_colour_manual(values = pallete4) + scale_fill_manual(values = pallete4) +
  geom_col(data = subset(qPCR_dko_pact_adar_sum, Gene %in% c("NFKBIE", "IRF1", "BIRC3")), aes(sgRNA, mean_expression, colour = Gene, fill = Gene),
           width = 0.5, alpha = 0.2, linewidth = 0.3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = subset(qPCR_dko_pact_adar_sum, Gene %in% c("NFKBIE", "IRF1", "BIRC3")), 
                aes(x = sgRNA, colour = Gene, ymax = mean_expression + sd_expression, 
                    ymin = mean_expression - sd_expression, y = mean_expression), width=0.2, linewidth = 0.3,
                position = position_dodge(width = 0.7)) 
ggsave("/Manuscripts/PACT/qPCR_plots/pact_adar_dko_453_nfkb.tiff", height =1.8, width = 2.6, units = "in", dpi = 300)


brackets <- data.frame(x =  c(0.6,1.6,2.6,3.6,4.6, 5.6, 6.6, 7.6), 
                       x_end = c(1.4,2.4,3.4,4.4, 5.4, 6.4, 7.4, 8.4), 
                       y= c(1.4,1.4,1.4,1.4,1.1,1.1, 2, 2))


brackets$Gene = NA


ggplot(subset(qPCR_dko_pact_adar, Gene %in% c("NFKBIE", "IRF1", "BIRC3")), aes(Sample, Fold_Change, colour = Gene, group = Gene)) +
  geom_point(aes(shape = as.factor(Replicate)), size = 0.3, 
             position = position_jitterdodge(dodge.width = 0.7)) + theme_science() + 
  scale_x_discrete(limits = c("sgC1_shSCR", "sgNTA_shSCR", "sgPACT.1_shSCR", "sgPACT.2_shSCR",
                              "sgC1_shADAR", "sgNTA_shADAR", "sgPACT.1_shADAR", "sgPACT.2_shADAR"), 
                   labels = c("sgC1", "sgNTA", "sgPACT-1", "sgPACT-2",
                              "sgC1", "sgNTA", "sgPACT-1", "sgPACT-2")) + 
  labs(x = "", y = "Relative Expression", colour = "", fill = "", title = "NF-kB Target Genes (MDA-MB-453)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = c(0.2, 0.78), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + 
  scale_colour_manual(values = pallete4) + scale_fill_manual(values = pallete4) +
  geom_col(data = subset(qPCR_dko_pact_adar_sum, Gene %in% c("NFKBIE", "IRF1", "BIRC3")), aes(Sample, mean_expression, colour = Gene, fill = Gene),
           width = 0.5, alpha = 0.2, linewidth = 0.3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = subset(qPCR_dko_pact_adar_sum, Gene %in% c("NFKBIE", "IRF1", "BIRC3")), 
                aes(x = Sample, colour = Gene, ymax = mean_expression + sd_expression, 
                    ymin = mean_expression - sd_expression, y = mean_expression), width=0.2, linewidth = 0.3,
                position = position_dodge(width = 0.7))  +
  stat_pvalue_manual(data = nfkb_453_t, label = "stars", step.increase = 0.02, 
                     y.position = "y", label.size = 2.5, linewidth = 0.3, tip.length = 0.01) +
  annotate("segment", x = 1, xend = 4, y = 4.2, yend = 4.2, colour = "black") + 
  annotate("segment", x = 5, xend = 8, y = 4.2, yend = 4.2, colour = "black") +
  annotate("text", x = 2.5, y = 4.4, colour = "black", label = "shSCR", size = 3) + 
  annotate("text", x = 6.5, y = 4.4, colour = "black", label = "shADAR", size = 3) +
  geom_bracket(data = brackets, aes(xmin = x, xmax = x_end, y.position = y, label = ""), colour = "black")
ggsave("/Manuscripts/PACT/qPCR_plots/pact_adar_dko_453_nfkb_stats.tiff", height =2.2, width = 2.5, units = "in", dpi = 300)


#qPCR PKR-PACT dko ----
qpcr_pkr_pact_dko <- na.omit(qpcr_pkr_pact_dko)
qpcr_pkr_pact_dko$Fold_Change <- 2^-(qpcr_pkr_pact_dko$Fold_Change)
qpcr_pkr_pact_dko$Cell_line <- NULL

#group and summarise data
qpcr_pkr_pact_dko_sum <- dplyr::group_by(qpcr_pkr_pact_dko, Gene, Sample)

qpcr_pkr_pact_dko_sum<- dplyr::summarise(qpcr_pkr_pact_dko_sum, mean_expression = mean(Fold_Change), sd_expression = sd(Fold_Change))


qpcr_pkr_pact_dko$pathway <- ifelse(grepl("ATF3|ASNS|SESN2", qpcr_pkr_pact_dko$Gene), "ATF4 Targets", "NF-kB Targets")
qpcr_pkr_pact_dko_sum$pathway <- ifelse(grepl("ATF3|ASNS|SESN2", qpcr_pkr_pact_dko_sum$Gene), "ATF4 Targets", "NF-kB Targets")

qpcr_pkr_pact_dko$Gene <- factor(qpcr_pkr_pact_dko$Gene, levels = c("ASNS", "ATF3", "SESN2", "BIRC3", "IRF1", "NFKBIE"))
qpcr_pkr_pact_dko_sum$Gene <- factor(qpcr_pkr_pact_dko_sum$Gene, levels = c("ASNS", "ATF3", "SESN2", "BIRC3", "IRF1", "NFKBIE"))


qpcr_pkr_pact_dko$Sample <- factor(as.factor(qpcr_pkr_pact_dko$Sample), 
                             levels = c("sgC1_sgPACT","sgC1_sgNTA", "sgPKR_sgNTA", "sgPKR_sgPACT"))

str(qpcr_pkr_pact_dko)

#stats analysis as in the immunoblot PACT/PKR dko data

pkr_pact_qpcr_atf4 <- DunnettTest(Fold_Change ~ Sample, data = qpcr_pkr_pact_dko[which(qpcr_pkr_pact_dko$pathway == "ATF4 Targets"),])
pkr_pact_qpcr_atf4 <- as.data.frame(pkr_pact_qpcr_atf4$sgC1_sgPACT[,1:4])

out <- strsplit(as.character(row.names(pkr_pact_qpcr_atf4)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

pkr_pact_qpcr_atf4$group1 <- out$V1
pkr_pact_qpcr_atf4$group2 <- out$V2
pkr_pact_qpcr_atf4$pval <- makeStars(pkr_pact_qpcr_atf4$pval)
pkr_pact_qpcr_atf4$y <- max(qpcr_pkr_pact_dko$Fold_Change[qpcr_pkr_pact_dko$pathway == "ATF4 Targets"]*1.1)

pkr_pact_qpcr_atf4 <- subset(pkr_pact_qpcr_atf4, !pval == "ns")
pkr_pact_qpcr_atf4$pathway <- "ATF4 Targets"


pkr_pact_qpcr_nfkb <- DunnettTest(Fold_Change ~ Sample, data = qpcr_pkr_pact_dko[which(qpcr_pkr_pact_dko$pathway == "NF-kB Targets"),])
pkr_pact_qpcr_nfkb <- as.data.frame(pkr_pact_qpcr_nfkb$sgC1_sgPACT[,1:4])

out <- strsplit(as.character(row.names(pkr_pact_qpcr_nfkb)), "-", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

pkr_pact_qpcr_nfkb$group1 <- out$V1
pkr_pact_qpcr_nfkb$group2 <- out$V2
pkr_pact_qpcr_nfkb$pval <- makeStars(pkr_pact_qpcr_nfkb$pval)
pkr_pact_qpcr_nfkb$y <- max(qpcr_pkr_pact_dko$Fold_Change[qpcr_pkr_pact_dko$pathway == "NF-kB Targets"]*1.1)

pkr_pact_qpcr_nfkb <- subset(pkr_pact_qpcr_nfkb, !pval == "ns")
pkr_pact_qpcr_nfkb$pathway <- "NF-kB Targets"

pkr_pact_qpcr_p <- rbind(pkr_pact_qpcr_nfkb, pkr_pact_qpcr_atf4)
str(pkr_pact_qpcr_p)


bracket_atf4 <- data.frame(x =  c(0.6,1.6,2.6,3.6), x_end = c(1.4,2.4,3.4,4.4), y= c(2.2, 2.2, 26, 2.2),
                          pathway = "ATF4 Targets")
bracket_nfkb <- data.frame(x =  c(0.6,1.6,2.6,3.6), x_end = c(1.4,2.4,3.4,4.4), y= c(2.2, 2.2, 16, 2.2),
                           pathway = "NF-kB Targets")

brackets <- rbind(bracket_atf4, bracket_nfkb)

brackets$Gene = NA



ggplot(qpcr_pkr_pact_dko, aes(Sample, Fold_Change, colour = Gene, group = Gene)) +
  geom_point(aes(shape = as.factor(Replicate)), size = 0.3, 
             position = position_jitterdodge(dodge.width = 0.7)) + theme_science() + 
  scale_x_discrete(limits = c("sgC1_sgNTA", "sgPKR_sgNTA", "sgC1_sgPACT", "sgPKR_sgPACT"), labels = c("sgC1/sgNTA", "sgPKR/sgNTA", "sgC1/sgPACT", "sgPKR/sgPACT-2")) + 
  labs(x = "", y = "Relative Expression", colour = "", fill = "") + facet_grid2(cols = vars(pathway)) +
  geom_hline(yintercept = 1, size = 0.3, colour = "grey", linetype = "dashed") + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "right", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank()) +
  guides(shape = "none") + scale_colour_manual(values = pallete6) + scale_fill_manual(values = pallete6) +
  geom_col(data = qpcr_pkr_pact_dko_sum, aes(Sample, mean_expression, colour = Gene, fill = Gene),
           width = 0.5, alpha = 0.2, linewidth = 0.3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = qpcr_pkr_pact_dko_sum, 
                aes(x = Sample, colour = Gene, ymax = mean_expression + sd_expression, 
                    ymin = mean_expression - sd_expression, y = mean_expression), width=0.2, linewidth = 0.3,
                position = position_dodge(width = 0.7)) +
  stat_pvalue_manual(data = pkr_pact_qpcr_p, label = "pval", step.increase = 0.07, , step.group.by = "pathway", y.position = "y", label.size = 2.5, linewidth = 0.3, tip.length = 0.01) +
  geom_bracket(data = brackets, aes(xmin = x, xmax = x_end, y.position = y, label = ""), colour = "black") 
ggsave("/Manuscripts/PACT/qPCR_plots/pact_pkr_dko.tiff", height =2.5, width = 3, units = "in", dpi = 300)




#tumorignesis ----
#group and summarise data
tumor <- tumorigenesis_hcc1806

#group and summarise by day and sgRNA
tumor_sum <- dplyr::group_by(tumor, sgRNA, Day)

tumor_sum<- dplyr::summarise(tumor_sum, mean_volume = mean(Volume), sd_volume = sd(Volume))

tumor_sum$Mouse <- tumor_sum$sgRNA


#generate palette for plot
fc <- colorRampPalette(c(palette[4], palette[5]))
rc <- colorRampPalette(c(palette[2], palette[3]))

palette_t <- c(fc(4), rc(5), palette[5], palette[2])


p1 <- ggplot(tumor, aes(Day, Volume, colour = as.factor(Mouse), group = as.factor(Mouse))) + 
  geom_line(alpha = 0.7, linewidth = 0.3) + geom_point(alpha = 0.7, size = 0.3) + 
  scale_y_continuous(limits = c(0, 1350), expand = expansion(mult = c(0, 0.1))) + 
  scale_x_continuous(limits = c(0,27), expand = expansion(mult = c(0,0))) + guides(colour = "none") +
  scale_colour_manual(values =palette_t) + theme_science() + geom_vline(xintercept = 5, linewidth = 0.3) + 
  geom_line(data = tumor_sum, aes(Day, mean_volume, colour = Mouse), linewidth = 0.6) + 
  geom_pointrange(data = tumor_sum, aes(Day, mean_volume, ymax = mean_volume + sd_volume, ymin = mean_volume - sd_volume), size = 0.3) + 
  labs(y = "Tumor Volume (mm^3)", x = "Days post injection", colour = "") +
  geom_label(aes(x = 5, y = 750), label = "Doxycycline initiated", angle = 90, colour = "black", size = 2)
p1
ggsave("tumorigenesis.tiff", height = 3, width = 2.5, units = "in", dpi = 300)

ggplot(tumor, aes(Day, Volume, colour = as.factor(Mouse), group = as.factor(Mouse))) + 
  geom_line(alpha = 0.7, linewidth = 0.3) + geom_point(alpha = 0.7, size = 0.3) +
  scale_colour_manual(values =palette_t) + theme_science() +
  labs(y = "Tumor Volume (mm^3)", x = "Days post injection", colour = "") + theme(legend.position = "bottom")
ggsave("tumorigenesis_nosum.tiff", height = 3, width = 2.5, units = "in", dpi = 300)

#make data.frame for final day measurements

tumor_sub <- subset(tumor, Day == 26)
tumor_sum_sub <- subset(tumor_sum, Day == 26)

#perform t.test
t.test(Volume ~ sgRNA, data = tumor_sub)$p.value


p2 <- ggplot(tumor_sub, aes(sgRNA, Volume, colour = sgRNA)) +
  geom_jitter(width = 0.12, size = 0.3) + theme_science() + 
  scale_x_discrete(limits = c("sgNTA", "sgPACT-2")) + 
  labs(x = "", y = "Tumor Volume at End Point (mm^3)", colour = "", fill = "") +
  scale_y_continuous(limits = c(0, 1350), expand = expansion(mult = c(0, 0.1)), position = "right") + 
  geom_col(data = tumor_sum_sub, aes(sgRNA, mean_volume, fill = sgRNA, colour = sgRNA), width = 0.5, alpha = 0.2, linewidth = 0.3) +
  geom_errorbar(data = tumor_sum_sub, aes(x = sgRNA, ymax = mean_volume + sd_volume, ymin = mean_volume - sd_volume, y = mean_volume, colour = sgRNA), width=0.2, linewidth = 0.3) +
  scale_colour_manual(values = c(palette[5], palette[2])) + scale_fill_manual(values = c(palette[5], palette[2])) +
  geom_bracket(xmin = "sgNTA", xmax = "sgPACT-2", y.position = 1250, 
               label = round(t.test(Volume ~ sgRNA, data = tumor_sub)$p.value, 4), colour = "black",
               label.size = 2) 
p2 
ggsave("tumor_endpoint.tiff", height =2, width = 2.9, units = "in", dpi = 300)

prow <- p1 + p2 + plot_layout(widths = c(0.8, 0.35), guides = "collect") & theme_science(base_size = 7)
prow

ggsave("tumor.tiff", height = 2, width = 4.7, units = "in", dpi = 300)
ggsave("tumor.pdf", height = 2, width = 4.7, units = "in", dpi = 300)

