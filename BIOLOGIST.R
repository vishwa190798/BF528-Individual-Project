#calling and loading packages
install.packages('gplots')
install.packages('gridExtra')
install.packages('reshape2')
install.packages("ggpubr")
install.packages('pheatmap')
library(ggpubr)
library(gplots)
library(gridExtra)
library(reshape2)
library(tidyverse)
library(pheatmap)

#reading the tables
P_0_1 <- read.table("P0_1_cufflinks/genes.fpkm_tracking",header = T)
head(P_0_1)
P_0_2 <- read.table("/project/bf528/project_2/data/samples/P0_2/genes.fpkm_tracking",header = T)
AD_1 <- read.table("/project/bf528/project_2/data/samples/Ad_1/genes.fpkm_tracking",header = T)
AD_2 <- read.table("/project/bf528/project_2/data/samples/Ad_2/genes.fpkm_tracking",header = T)
P_4_1 <- read.table("/project/bf528/project_2/data/samples/P4_1/genes.fpkm_tracking",header = T)
P_4_2 <- read.table("/project/bf528/project_2/data/samples/P4_2/genes.fpkm_tracking",header = T)
P_7_1 <- read.table("/project/bf528/project_2/data/samples/P7_1/genes.fpkm_tracking",header = T)
P_7_2 <- read.table("/project/bf528/project_2/data/samples/P7_2/genes.fpkm_tracking",header = T)

#reading genes from paper for sarcomere, mitochondria and cell cycle
sarcomere_gene <- c("Pdlim5","Pygm","Myoz2","Des","Csrp3","Tcap","Cryab")
mito_gene <- c("Brp44l","Prdx3","Acat1","Echs1","Slc25a11","Phyh")
cell_gene <- c("Cdc7","E2f8","Cdk7","Cdc26","Cdc6","E2f1","Cdc27","6720463M24Rik","Cdc45","Rad51","Aurkb","Cdc23")

#subset sarcomere genes between all the samples
sarco_fpkm <- c(P_0_1[P_0_1$gene_short_name %in% sarcomere_gene,"FPKM"],
                P_0_2[P_0_2$gene_short_name %in% sarcomere_gene,"FPKM"],
                AD_1[AD_1$gene_short_name %in% sarcomere_gene,"FPKM"],
                AD_2[AD_2$gene_short_name %in% sarcomere_gene,"FPKM"],
                P_4_1[P_4_1$gene_short_name %in% sarcomere_gene,"FPKM"],
                P_4_2[P_4_2$gene_short_name %in% sarcomere_gene,"FPKM"],
                P_7_1[P_7_1$gene_short_name %in% sarcomere_gene,"FPKM"],
                P_7_2[P_7_2$gene_short_name %in% sarcomere_gene,"FPKM"])

#converting the subset to table
sarcomere_table <- data.frame("FPKM" = sarco_fpkm, "Sample" = rep(c("p0","ad","p4","p7"), each = 14),
                          "Gene" = rep(P_0_1[P_0_1$gene_short_name %in% sarcomere_gene, "gene_short_name"], 8))

#creating new column with levels
sarcomere_table$Sample <- factor(sarcomere_table$Sample, levels = c("p0","p4","p7","ad"))

#creating line graph for sarcomere genes
sarco_plot <- sarcomere_table %>% group_by(Sample, Gene) %>% summarise("mean" = mean(FPKM), "Sd" = sd(FPKM)) %>% ggplot(aes(x = Sample,y = mean, group = Gene)) + geom_line(aes(color = Gene)) + geom_point(aes(color = Gene)) + labs(title="Sarcomere", x ="Sample", y = "FPKM") + theme_bw()

##subset mitochondria genes between all the samples
mito_fpkm <- c(P_0_1[P_0_1$gene_short_name %in% mito_gene, "FPKM"],
               P_0_2[P_0_2$gene_short_name %in% mito_gene,"FPKM"],
               AD_1[AD_1$gene_short_name %in% mito_gene,"FPKM"],
               AD_2[AD_2$gene_short_name %in% mito_gene,"FPKM"],
               P_4_1[P_4_1$gene_short_name %in% mito_gene,"FPKM"],
               P_4_2[P_4_2$gene_short_name %in% mito_gene,"FPKM"],
               P_7_1[P_7_1$gene_short_name %in% mito_gene,"FPKM"],
               P_7_2[P_7_2$gene_short_name %in% mito_gene,"FPKM"])

##converting the subset to table
mito_table <- data.frame("FPKM" = mito_fpkm, "Sample" = rep(c("p0","ad","p4","p7"), each = 12),
                         "Gene" = rep(P_0_1[P_0_1$gene_short_name %in% mito_gene,"gene_short_name"], 8))

##creating new column with levels
mito_table$Sample <- factor(mito_table$Sample, levels = c("p0","p4","p7","ad"))

##creating line graph for mitochondria genes
mito_plot <- mito_table %>% group_by(Sample, Gene) %>% summarise("mean" = mean(FPKM),"Sd" = sd(FPKM)) %>%ggplot(aes(x = Sample,y = mean, group = Gene)) + geom_line(aes(color = Gene)) + geom_point(aes(color = Gene)) + labs(title="Mitochondria", x ="Sample", y = "FPKM") + scale_color_discrete(name = "Gene", labels = c("Acat1","Mpc1","Echs1", "Phyh","Prdx3","Slc25a11")) + theme_bw()

#subset cell cycle genes between all the samples
cell_fpkm <- c(P_0_1[P_0_1$gene_short_name %in% cell_gene,"FPKM"],
               P_0_2[P_0_2$gene_short_name %in% cell_gene,"FPKM"],
               AD_1[AD_1$gene_short_name %in% cell_gene,"FPKM"],
               AD_2[AD_2$gene_short_name %in% cell_gene,"FPKM"],
               P_4_1[P_4_1$gene_short_name %in% cell_gene,"FPKM"],
               P_4_2[P_4_2$gene_short_name %in% cell_gene,"FPKM"],
               P_7_1[P_7_1$gene_short_name %in% cell_gene,"FPKM"],
               P_7_2[P_7_2$gene_short_name %in% cell_gene,"FPKM"])

#converting the subset to table
cell_table <- data.frame("FPKM" = cell_fpkm, "Sample" = rep(c("p0","ad","p4","p7"), each = 24),"Gene" = rep(P_0_1[P_0_1$gene_short_name %in% cell_gene,"gene_short_name"], 8))

#creating new column with levels
cell_table$Sample <- factor(cell_table$Sample, levels = c("p0","p4","p7","ad"))

new_cell_gene <- cell_gene %>% sort()
new_cell_gene[1] <- "Bora"

#creating line graph for cell cycle genes
cell_plot <- cell_table %>% group_by(Sample, Gene) %>% summarise("mean" = mean(FPKM),"Sd" = sd(FPKM)) %>% ggplot(aes(x = Sample,y = mean, group = Gene)) +geom_line(aes(color = Gene)) +geom_point(aes(color = Gene)) +labs(title="Cell Cycle",x ="Sample", y = "FPKM") + scale_color_discrete(name = "Gene", labels = new_cell_gene) + guides(color=guide_legend(ncol=2)) + theme_bw()

#arranging the plots
ggarrange(sarco_plot, mito_plot, cell_plot, labels = c("A", "B","C"), nrow = 1)

ref_table <- read.table("/project/bf528/project_2/data/fpkm_matrix.csv", header = T)
head(ref_table)

cuffdiff <- read.table("cuffdiff_out/gene_exp.diff", header = T)
head(cuffdiff)

#sorting the DEG on q-value and extracting top 1000 genes
sorted <- order(cuffdiff$q_value, decreasing = F)
top_1000_gene <- cuffdiff[sorted,"gene"] %>% unique()
top_1000_gene <- top_1000_gene[1:1000]

top_id <- P_0_1[P_0_1$gene_short_name %in% top_1000_gene,"tracking_id"] %>% unique()

gene_frame <- P_0_1[, c("tracking_id", "FPKM")] %>%inner_join(ref_table, by = "tracking_id")
colnames(gene_frame)[2] <- "P0_1_FPKM"
head(gene_frame)
top_1000 <- gene_frame[gene_frame$tracking_id %in% top_id,c(1:5)]
head(top_1000)
top_1000 <- top_1000[rowSums(top_1000[,-1]) != 0, ]

#creating heatmap
pheatmap(as.matrix(top_1000[,-1]), scale = "row",show_rownames = F, labels_col = c("P0_1","P0_2","Ad_1","Ad_2")
