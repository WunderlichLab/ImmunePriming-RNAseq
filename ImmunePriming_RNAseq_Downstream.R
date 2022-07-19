#### Required Packages ####

library(edgeR)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library("AnnotationDbi")
library("org.Dm.eg.db")
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(corrplot)
library(gplots)
library(circlize)
library(ggpubr)
library(VennDiagram)
library(sva)

#### Colors For Plotting ####
# Changing the colors attributed to each condition for your plots
CondColors <- c("gray","grey35","#648287","#83bd39","#1cd9b6","#4b6e1e","#a7f542","#93c1c9","#355a78",
                "grey70","grey1")
names(CondColors) <- c("NoInj_d0", "NoInj_d7", "DblPBS_d8", 
                       "EfaeLo_d7", "EfaeTrain_d8", 
                       "EfaeHi_d1", "EfaeLo_d1",
                       "PBS_d1", "EfaeMock_d8",
                       "NoInj(1FB)_d0", "NoInj(6FB)_d0")
CondColors <- list(FB_Condition = CondColors)

FB_logFCColors <- c("#83bd39","#1cd9b6","#4b6e1e","#a7f542","#355a78")
names(FB_logFCColors) <- c("FB_EFAELO_d7", "FB_EFAETRAIN_d8", 
                        "FB_EFAEHI_d1", "FB_EFAELO_d1","FB_EFAEMOCK_d8")
FB_logFCColors <- list(Condition = FB_logFCColors)

#### Reading in the Data and Descriptions about the Data -------------------------------------------------------------------
setwd("/Users/cabweebwub/Desktop/FBCounts/")

# Make a dataframe that contains all the metadata about your samples, used for plotting later
FB_sampleInfo <- readr::read_csv("/Users/cabweebwub/Desktop/FB_PlayRoom.csv",
                                   col_names=TRUE)
FB_sampleInfo <- unite(FB_sampleInfo, col= "Condition", c("Condition","Time"), sep = "_")
FB_sampleInfo <- FB_sampleInfo[order(FB_sampleInfo$Sample),]

# Make objects for the specific variables in your metadata
FB_Genotype <- as.factor(FB_sampleInfo$Genotype)
FB_Condition <- as.factor(FB_sampleInfo$Condition)
FB_Batch <- as.factor(FB_sampleInfo$Batch)
FB_Age <- as.factor(FB_sampleInfo$Age)

FB_sample2cond <- as.data.frame(FB_Condition)
row.names(FB_sample2cond) <- FB_sampleInfo$Sample

#### Reading data into EdgeR for differential gene expression analysis ------------------------------------------------------

# Read in your files and group them acording to their condition
FB_DGE <- readDGE(files = list.files(), columns = c(1,7), skip = 1 , group = FB_Condition)
FB_DGE$genes <- read.delim("/Users/cabweebwub/Desktop/Counts/KC20-0100_counts.txt", skip = 1) %>%
  .[,c("Geneid","Length")]

# Make a design matrix that has your batch & condition data
FB_design <- model.matrix(~0+FB_Condition)
rownames(FB_design) <- colnames(FB_DGE)

#### Filtering out non-informative data -------------------------------------------------------------------------------------

# Filtering out non-informative & testes-specific genes
FB_keep <- filterByExpr(FB_DGE, design = FB_design)

FB_DGE <- FB_DGE[FB_keep, , keep.lib.sizes=FALSE] #Setting lib.sizes to F re-calculates lib size after filtering

FB_DGE$genes <- dplyr::filter(FB_DGE$genes, FB_DGE$genes$Geneid %in% row.names(FB_DGE$counts)) #removes genes that have been filtered out from your gene annotations (needed for RPKM calc later)


#### Normalizing your data -------------------------------------------------------------------------------------

# Normalizes your data using the Trimmed Mean of M Values (TMM)
FB_DGE <- calcNormFactors(FB_DGE, method = "TMM")

# Estimate dispersion
FB_DGE <- estimateDisp(FB_DGE, design = FB_design, robust=TRUE)

#### Making a data frame with TPMs & CPMs -------------------------------------------------------------------------------------

# Function for calculating TPMs
calc_tpm <- function(x, gene.length) {
  x <- as.matrix(x)
  len.norm.lib.size <- colSums(x / gene.length)
  return((t(t(x) / len.norm.lib.size) * 1e06) / gene.length)
}

# Manually calculating out the TPMs
FB_tpm <- calc_tpm(FB_DGE, gene.length = FB_DGE$gene$Length)
colnames(FB_tpm) <- sub("_counts", "", colnames(FB_tpm))
FB_tpm <- as.data.frame(FB_tpm)
FB_tpm <- FB_tpm[grep('^FB', rownames(FB_tpm)),]
FB_tpm$name = mapIds(org.Dm.eg.db,
                       keys=row.names(FB_tpm), 
                       column="GENENAME",
                       keytype="ENSEMBL",
                       multiVals="first")
FB_tpm$symbol <- mapIds(org.Dm.eg.db,
                          keys=row.names(FB_tpm),
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
FB_tpm$flybaseID <- row.names(FB_tpm)

# Make a new dataframe with log2(CPM+1) values
FB_cpm <- cpm(FB_DGE, normalized.lib.sizes = TRUE)
FB_cpm <- log2(FB_cpm+1)
colnames(FB_cpm) <- sub("_counts", "", colnames(FB_cpm))
FB_cpm <- as.data.frame(FB_cpm)
FB_cpm <- FB_cpm[grep('^FB', rownames(FB_cpm)),] #Removes non-gene related rows

FB_cpm$name = mapIds(org.Dm.eg.db,
                       keys=row.names(FB_cpm), 
                       column="GENENAME",
                       keytype="ENSEMBL",
                       multiVals="first")
FB_cpm$symbol <- mapIds(org.Dm.eg.db,
                          keys=row.names(FB_cpm),
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
FB_cpm$flybaseID <- row.names(FB_cpm)

#### PCA Plot ####
## PCA Plot ##
FB_data4PCA <- FB_DGE$counts
colnames(FB_data4PCA) <- sub("_counts.txt", "", colnames(FB_data4PCA))

FB_res.pca <- PCA(t(data.matrix(FB_data4PCA)), scale.unit = T, ncp = 5, graph = F)
fviz_pca_ind(FB_res.pca, axes = c(1,2), repel = TRUE, 
             pointsize = 2.5)

## PCA Plot w/ ggplot ##
FB_PCA_info <- as.data.frame(FB_res.pca$ind$coord)
rownames(FB_PCA_info) <- sub("_counts", "", rownames(FB_PCA_info))
FB_PCA_info <- merge(FB_PCA_info, FB_sample2cond, by = "row.names", all = T)
FB_PCA_info$Age <- FB_sampleInfo$Age
FB_PCA_info$Batch <- FB_sampleInfo$Batch
FB_PCA_info <- na.omit(FB_PCA_info)

# Plotting PCA w/ ggplot2 instead (it looks nicer)
#png(filename = "PCA2.png",width = 1900, height=1350, res=300)
ggplot(FB_PCA_info, aes(x=`Dim.1`, y=`Dim.2`, color = FB_PCA_info$FB_Condition)) + 
  theme(panel.background = element_blank(),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16, face = "bold")) +
  geom_point(size = 5) + 
  geom_hline(yintercept = 0, linetype="dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype="dashed", color = "grey") +
  geom_text(label=FB_PCA_info$Row.names, nudge_y = 3.5, alpha = 0.7) +
  xlab(paste("Dim1 (",round(FB_res.pca$eig[1,2], digits=1),"% var)", sep="")) +
  ylab(paste("Dim2 (",round(FB_res.pca$eig[2,2], digits=1),"% var)", sep="")) + 
  labs(shape = "Batch") +
  scale_colour_manual(values=CondColors$FB_Condition,
                      name = "Condition",
                      limits = c("NoInj_d0","NoInj_d7",
                                 "PBS_d1","DblPBS_d8",
                                 "EfaeHi_d1",
                                 "EfaeLo_d1","EfaeLo_d7",
                                 "EfaeMock_d8",
                                 "EfaeTrain_d8"),
                      labels = c("No Injection - d0","No Injection - d7",
                                 "PBS - d1","Double PBS - d8",
                                 "Efae HiDose - d1",
                                 "EfaeLo - d1","EfaeLo - d7",
                                 "Efae Mock-Trained - d8",
                                 "Efae Trained - d8"))
#### Correlation Plot ####
FB_TPM_reordered <- FB_tpm[c("KC20-0100","KC20-0201","KC20-0202",
                             "KC20-0204","KC20-0234","KC20-0404",
                             "KC20-0397","KC20-0401","KC21-0110",
                             "KC20-0393","KC20-0400","KC20-0414",
                             "KC20-0394","KC20-0407","KC20-0408","KC21-0113","KC21-0135",
                             "KC21-0138","KC21-0139","KC21-0140",
                             "KC20-0402","KC21-0111","KC21-0112")]

# How many breaks you want in your color scale (usually change it so it matches the max and min of my data)
SpearmanBreaks <- seq(from = 0.5, to = 1.0, by = c(0.01))
# Check out how many values are in myBreaks and change the paletteLength
SpearmanpaletteLength <- length(SpearmanBreaks)
# Choose what colors you want as your lower, middle, and upper cutoffs
SpearmanColor <- colorRampPalette(c("blue","white","red"))(SpearmanpaletteLength)

png(filename = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/Figures/FB_Spearman.png",width = 1700, height=1500, res=300)
pheatmap(cor(as.matrix(FB_TPM_reordered), method = "spearman"),
         annotation_col = FB_sample2cond,
         annotation_row = FB_sample2cond,
         color=SpearmanColor, breaks=SpearmanBreaks,
         cluster_rows = F,
         cluster_cols = F,
         annotation_colors = CondColors,
         show_rownames = T,
         show_colnames = T,
         annotation_names_row = F,
         annotation_names_col = F,
         display_numbers = F,
         gaps_row = c(3,6,9,12,17,20),
         gaps_col = c(3,6,9,12,17,20),
         cex = 0.75,
         annotation_legend = F,
         fontsize = 15)
dev.off()

#### DGE Testing ####

FB_EFAEHI_d1vsCTRL_d0 <- exactTest(object = FB_DGE, pair = c("NoInj_d0","EfaeHi_d1"))
FB_EFAELO_d1vsCTRL_d0 <- exactTest(object = FB_DGE, pair = c("NoInj_d0","EfaeLo_d1"))
FB_EFAELO_d7vsCTRL_d7 <- exactTest(object = FB_DGE, pair = c("NoInj_d7","EfaeLo_d7"))
FB_EFAETRAIN_d8vsCTRL_d7 <- exactTest(object = FB_DGE, pair = c("NoInj_d7","EfaeTrain_d8"))
FB_EFAEMOCK_d8vsCTRL_d7 <- exactTest(object = FB_DGE, pair = c("NoInj_d7","EfaeMock_d8"))
FB_EFAELO_d1vsEFAELO_d7 <- exactTest(object = y, pair = c("EfaeLo_d7","EfaeLo_d1"))
FB_EFAETRAIN_d8vsEFAELO_d7 <- exactTest(object = FB_DGE, pair = c("EfaeTrain_d8","EfaeLo_d7"))
FB_EFAETRAIN_d8vsEFAELO_d1 <- exactTest(object = FB_DGE, pair = c("EfaeTrain_d8","EfaeLo_d1"))
FB_EFAETRAIN_d8vsEFAEHI_d1 <- exactTest(object = FB_DGE, pair = c("EfaeTrain_d8","EfaeHi_d1"))
FB_EFAETRAIN_d8vsEFAEMOCK_d8 <- exactTest(object = FB_DGE, pair = c("EfaeTrain_d8","EfaeMock_d8"))

## Extracting the top n hits from your comparisons ##
FBRes_EFAEHI_d1vsCTRL_d0 <- topTags(FB_EFAEHI_d1vsCTRL_d0,5000)
FBRes_EFAELO_d1vsCTRL_d0 <- topTags(FB_EFAELO_d1vsCTRL_d0,5000)
FBRes_EFAELO_d7vsCTRL_d7 <- topTags(FB_EFAELO_d7vsCTRL_d7,5000)
FBRes_EFAETRAIN_d8vsCTRL_d7 <- topTags(FB_EFAETRAIN_d8vsCTRL_d7,5000)
FBRes_EFAEMOCK_d8vsCTRL_d7 <- topTags(FB_EFAEMOCK_d8vsCTRL_d7,5000)
FBRes_EFAELO_d1vsEFAELO_d7 <- topTags(FB_EFAELO_d1vsEFAELO_d7,5000)
FBRes_EFAETRAIN_d8vsEFAELO_d7 <- topTags(FB_EFAETRAIN_d8vsEFAELO_d7,5000)
FBRes_EFAETRAIN_d8vsEFAELO_d1 <- topTags(FB_EFAETRAIN_d8vsEFAELO_d1,5000)
FBRes_EFAETRAIN_d8vsEFAEHI_d1 <- topTags(FB_EFAETRAIN_d8vsEFAEHI_d1,5000)
FBRes_EFAETRAIN_d8vsEFAEMOCK_d8 <- topTags(FB_EFAETRAIN_d8vsEFAEMOCK_d8,5000)

#### Volcano Plot ####

volcano_plot <- FBRes_EFAETRAIN_d8vsEFAEMOCK_d8$table
ggplot() +
  geom_point(volcano_plot, 
             mapping = aes(x=logFC, y=-log10(FDR)), color = "grey29") +
  geom_point(subset(volcano_plot, FDR<.05), 
             mapping = aes(x=logFC, y=-log10(FDR)),color = "grey35") +
  geom_point(subset(volcano_plot, abs(logFC)>2), 
             mapping = aes(x=logFC, y=-log10(FDR)),color = "grey35") +
  geom_point(subset(volcano_plot, abs(logFC)>2 & FDR<.05), 
             mapping = aes(x=logFC, y=-log10(FDR)),color = "magenta") +
  geom_point(subset(volcano_plot, 
                    rownames(volcano_plot) %in% FB_overlap$a3), 
             mapping = aes(x=logFC, y=-log10(FDR)),color = "green") +
  ggtitle("FBs MockPrimed-d8 vs EfaePrimed-d8") +
  theme(panel.background = element_blank(),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16, face = "bold"),
        axis.line = element_line(size = 1, colour = "grey25"),
        plot.title = element_text(size=18, face="bold", hjust = 0.5)) + 
  geom_hline(yintercept = 1.30, linetype="twodash", color = "red", size = 1.5) +
  geom_vline(xintercept = 2, linetype="twodash", color = "blue", size = 1.5) +
  geom_vline(xintercept = -2, linetype="twodash", color = "blue", size = 1.5)

sigstats_compare <- FBRes_EFAETRAIN_d8vsEFAEMOCK_d8$table
sigstats_name <- "EfaePrimed-d8_vs_MockPrimed-d8"

### How many genes are significantly DE'd?
sigstats <- subset(sigstats_compare, FDR<0.05)
nrow(sigstats)

## Saving the top n tags as a table on your comp ##
sigstats %>% 
  .[order(-.$logFC),] %>%
  head(. , 200) %>%
  row.names(.) %>%
  as_tibble(.) %>%
  readr::write_tsv(., 
                   path = paste0("/Users/cabweebwub/Desktop/Efae RNA-seq Paper/FB_top200_", sigstats_name, ".tsv", sep=""))


#### Gene Lists ####
immune_list <- readr::read_tsv("/Users/cabweebwub/Documents/Wunderlich Lab Docs/RData/GeneLists/Immunelist_LB_2019_cat.txt")
immune_effetors_list <- subset(immune_list, Involvement=="effector")
immune_signaling_list <- subset(immune_list, Involvement=="signaling")
immune_detection_list <- subset(immune_list, Involvement=="detection")
immune_CrebAPapeCore_list <- subset(immune_list, Core==TRUE)

neg_IMD_list <- readr::read_tsv("/Users/cabweebwub/Documents/Wunderlich Lab Docs/RData/GeneLists/NegativeRegulationofIMD_genelist.txt")
pos_IMD_list <- readr::read_csv("/Users/cabweebwub/Documents/Wunderlich Lab Docs/RData/GeneLists/PositiveRegulatorsofIMD_list_FBgg0001197.csv")
IMD_list <- readr::read_tsv("/Users/cabweebwub/Documents/Wunderlich Lab Docs/RData/GeneLists/imd_list.txt")
sperm_list <- readr::read_tsv("/Users/cabweebwub/Documents/Wunderlich Lab Docs/RData/GeneLists/SpermGenes.txt")
JNKPosReg_list <- readr::read_tsv("/Users/cabweebwub/Documents/Wunderlich Lab Docs/RData/GeneLists/JNK_PosReg.txt")
Bomanins_list <- readr::read_tsv("/Users/cabweebwub/Documents/Wunderlich Lab Docs/RData/GeneLists/Bomanins.txt")
IMD_Masterlist <- read_tsv("/Users/cabweebwub/Documents/Wunderlich Lab Docs/RData/GeneLists/IMD_MasterList.txt") #My own annotated version


#### Making the Dataframe that has logFC values of Conditions over Control ####

FB_tru_logFC <- data.frame(FB_EFAEHI_d1vsCTRL_d0$table$logFC, 
                             FB_EFAELO_d1vsCTRL_d0$table$logFC, 
                             FB_EFAELO_d7vsCTRL_d7$table$logFC,
                             FB_EFAEMOCK_d8vsCTRL_d7$table$logFC,
                             FB_EFAETRAIN_d8vsCTRL_d7$table$logFC)
row.names(FB_tru_logFC) <- make.names(row.names(FB_EFAEHI_d1vsCTRL_d0$table))
colnames(FB_tru_logFC) <- sub("vsCTRL_d0.table.logFC", "", colnames(FB_tru_logFC))
colnames(FB_tru_logFC) <- sub("vsCTRL_d7.table.logFC", "", colnames(FB_tru_logFC))
FB_tru_logFC$name = mapIds(org.Dm.eg.db,
                             keys=row.names(FB_tru_logFC), 
                             column="GENENAME",
                             keytype="ENSEMBL",
                             multiVals="first")
FB_tru_logFC$symbol <- mapIds(org.Dm.eg.db,
                                keys=row.names(FB_tru_logFC),
                                column="SYMBOL",
                                keytype="ENSEMBL",
                                multiVals="first")
FB_tru_logFC$flybaseID <- row.names(FB_tru_logFC)
FB_tru_logFC$name <- make.names(FB_tru_logFC$name, unique = T)
FB_tru_logFC$symbol <- make.names(FB_tru_logFC$symbol, unique = T)

#### TPM Heatmap of a List of Genes ####

## Make the scale as you like
# How many breaks you want in your color scale (usually change it so it matches the max and min of my data)
TPMmyBreaks <- seq(from = 0, to = 10.2, by = c(0.1))
# Check out how many values are in myBreaks and change the paletteLength
TPMpaletteLength <- length(TPMmyBreaks)
# Chose what colores you want as your lower, middle, and upper cutoffs
TPMmyColor <- colorRampPalette(c("white","darkslategray"))(TPMpaletteLength)

FB_TPM_Heatmap <- FB_tpm[c("KC20-0100","KC20-0201","KC20-0202",
                           "KC20-0204","KC20-0234","KC20-0404",
                           "KC20-0397","KC20-0401","KC21-0110",
                           "KC20-0393","KC20-0400","KC20-0414",
                           "KC20-0394","KC20-0407","KC20-0408","KC21-0113","KC21-0135",
                           "KC21-0138","KC21-0139","KC21-0140",
                           "KC20-0402","KC21-0111","KC21-0112",
                               "name","symbol","flybaseID")]

### TPM Heatmap ###
AMP_TPM <- FB_TPM_Heatmap %>%
  dplyr::filter(.,flybaseID %in% AMPs_list$FlyBaseID)
AMP_TPM[11,"name"] <- "Cecropin A2"
AMP_TPM[6,"name"] <- "Attacin-A"
AMP_TPM[8,"name"] <- "Attacin-C"
AMP_TPM <- AMP_TPM %>%
  column_to_rownames(.,"name")%>%
  dplyr::select(.,starts_with("KC"))

Neg_IMD_TPM <- FB_TPM_Heatmap %>%
  dplyr::filter(.,flybaseID %in% neg_IMD_list$FBID_KEY)
Neg_IMD_TPM[12,"symbol"] <- "PGRP-SC1b"
Neg_IMD_TPM <- Neg_IMD_TPM %>%
  column_to_rownames(.,"symbol")%>%
  dplyr::select(.,starts_with("KC"))

Toll_TPM <- FB_TPM_Heatmap %>%
  dplyr::filter(.,flybaseID %in% TollPathway_list$FlyBaseID)
Toll_TPM[29,"symbol"] <- "PGRP-SC1b"
Toll_TPM <- Toll_TPM %>%
  column_to_rownames(.,"symbol")%>%
  dplyr::select(.,starts_with("KC"))

MAPK_TPM <- FB_TPM_Heatmap %>%
  dplyr::filter(.,symbol %in% c("Raf","pirk","puc","Src64B","Src42A","Sac1","CG8184","Traf.like",
                                "PVf1","Tao","Mkp3")) %>%
  column_to_rownames(.,"symbol") %>%
  dplyr::select(.,starts_with("KC"))


IMD_effector_TPM <- FB_TPM_Heatmap %>%
  dplyr::filter(.,symbol %in% c("DptB", "DptA")) %>%
  column_to_rownames(.,"symbol") %>%
  dplyr::select(.,starts_with("KC"))

Sperm_TPM <- FB_TPM_Heatmap %>%
  dplyr::filter(.,flybaseID %in% sperm_list$FlybaseID) %>%
  column_to_rownames(.,"flybaseID") %>%
  dplyr::select(.,starts_with("KC"))

### Making the color bar for your heatmap ###
png(filename = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/Figures/FB_Sperm_TPM.png",width = 1000, height=600, res=150)
pheatmap(log2(as.matrix(Sperm_TPM)+1),
         color=TPMmyColor, 
         breaks=TPMmyBreaks,
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = FB_sample2cond,
         annotation_colors = CondColors,
         show_rownames = F,
         annotation_legend = F,
         annotation_names_col = F,
         gaps_col = c(3,6,9,12,17,20))
dev.off()

#### TPM Plot of a List of Genes ####
list_TPM <- dplyr::filter(FB_TPM_Heatmap, flybaseID %in% 
                            c("FBgn0014865", "FBgn0283461", "FBgn0262881", "FBgn0034407",
                              "FBgn0067905", "FBgn0040653", "FBgn0053470")) %>%
  dplyr::select(.,starts_with("KC")) %>%
  colSums(.)

list_TPM <- dplyr::filter(FB_TPM_Heatmap, flybaseID %in% 
                            c("FBgn0004240","FBgn0034407","FBgn0012042","FBgn0041581",
                              "FBgn0041579","FBgn0038530", "FBgn0010388", "FBgn0010385",
                              "FBgn0000276","FBgn0000277")) %>%
  dplyr::select(.,starts_with("KC")) %>%
  colSums(.)


Gene_Plot <- (list_TPM)
Gene_Plot <- rownames_to_column(as.data.frame(Gene_Plot), var = "Sample")
Gene_Plot <- as.data.frame(Gene_Plot)
Gene_Plot <- rownames_to_column(FB_sample2cond, var = "Sample") %>%
  as.data.frame(.) %>%
  full_join(.,Gene_Plot, by = "Sample")
colnames(Gene_Plot) <- c("Sample","FB_Condition","Gene")
Gene_Plot <- group_by(Gene_Plot, FB_Condition) %>%
  summarise(Mean = exp(mean(log(Gene))), Se = sd(Gene)/sqrt(n()))
Gene_Plot <- Gene_Plot %>%
  mutate(FB_Condition =  factor(FB_Condition, levels = c("NoInj_d0","NoInj_d7",
                                                         "EfaeLo_d1", "EfaeHi_d1",
                                                         "EfaeMock_d8",
                                                         "EfaeTrain_d8"))) %>%
  arrange(FB_Condition)
Gene_Plot <- na.omit(Gene_Plot)

Gene_Plot_Reps <- (list_TPM)
Gene_Plot_Reps <- rownames_to_column(as.data.frame(Gene_Plot_Reps), var = "Sample")
Gene_Plot_Reps <- as.data.frame(Gene_Plot_Reps)
Gene_Plot_Reps <- rownames_to_column(FB_sample2cond, var = "Sample") %>%
  as.data.frame(.) %>%
  full_join(.,Gene_Plot_Reps, by = "Sample")
colnames(Gene_Plot_Reps) <- c("Sample","FB_Condition","Gene")
Gene_Plot_Reps <- Gene_Plot_Reps %>%
  mutate(FB_Condition =  factor(FB_Condition, levels = c("NoInj_d0","NoInj_d7",
                                                         "EfaeLo_d1", "EfaeHi_d1",
                                                         "EfaeMock_d8",
                                                         "EfaeTrain_d8"))) %>%
  arrange(FB_Condition)
Gene_Plot_Reps <- na.omit(Gene_Plot_Reps)


#### TPM Plot of a Specific Gene ####

# Start by making a data frame with the info for the gene you are interested in
SingleGene <- "nub"

SingleGene_TPM <- FB_TPM_Heatmap %>%
  dplyr::filter(.,symbol==SingleGene) %>%
  column_to_rownames(.,"symbol") %>%
  dplyr::select(.,starts_with("KC"))

SingleGene_TPM <- FB_TPM_Heatmap %>%
  dplyr::filter(.,flybaseID== "FBgn0288447") %>%
  remove_rownames(.) %>%
  column_to_rownames(.,"symbol") %>%
  dplyr::select(.,starts_with("KC"))

Gene_Plot <- t(SingleGene_TPM)
Gene_Plot <- rownames_to_column(as.data.frame(Gene_Plot), var = "Sample")
Gene_Plot <- as.data.frame(Gene_Plot)
Gene_Plot <- rownames_to_column(FB_sample2cond, var = "Sample") %>%
  as.data.frame(.) %>%
  full_join(.,Gene_Plot, by = "Sample")
colnames(Gene_Plot) <- c("Sample","FB_Condition","Gene")
Gene_Plot <- group_by(Gene_Plot, FB_Condition) %>%
  summarise(Mean = mean(Gene), Se = sd(Gene)/sqrt(n()))
Gene_Plot <- Gene_Plot %>%
  mutate(FB_Condition =  factor(FB_Condition, levels = c("NoInj_d7",
                                                         "EfaeMock_d8",
                                                         "EfaeTrain_d8"))) %>%
  arrange(FB_Condition)
Gene_Plot <- na.omit(Gene_Plot)

Gene_Plot_Reps <- t(SingleGene_TPM)
Gene_Plot_Reps <- rownames_to_column(as.data.frame(Gene_Plot_Reps), var = "Sample")
Gene_Plot_Reps <- as.data.frame(Gene_Plot_Reps)
Gene_Plot_Reps <- rownames_to_column(FB_sample2cond, var = "Sample") %>%
  as.data.frame(.) %>%
  full_join(.,Gene_Plot_Reps, by = "Sample")
colnames(Gene_Plot_Reps) <- c("Sample","FB_Condition","Gene")
Gene_Plot_Reps <- Gene_Plot_Reps %>%
  mutate(FB_Condition =  factor(FB_Condition, levels = c("NoInj_d7",
                                                         "EfaeMock_d8",
                                                         "EfaeTrain_d8"))) %>%
  arrange(FB_Condition)
Gene_Plot_Reps <- na.omit(Gene_Plot_Reps)

# Plot the TPM values for your gene of interest
png(filename = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/Figures/FB_TPM_FBCommonEfaeBoms.png",width = 600, height=700, res=100)
ggplot(Gene_Plot) +
  geom_col(mapping = aes(x = FB_Condition, y= Mean, fill = FB_Condition)) +
  geom_jitter(data = Gene_Plot_Reps, mapping = aes(x = FB_Condition, y= Gene),
              cex=4, alpha = 0.6, position = position_dodge(0.89)) +
  geom_errorbar(aes(x = FB_Condition, ymin = Mean - Se, ymax = Mean + Se), 
                width = .2, position = position_dodge(.9)) +
  theme(panel.background = element_blank(),
        axis.title=element_text(size=20,face="bold"),
        legend.position = "none",
        legend.text=element_text(size=18),
        legend.title = element_text(size=20, face = "bold"),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=16),
        axis.line = element_line(size = 1, colour = "grey25")) +
  ylab(" ") +
  xlab(" ") +
  ylim(0, 100000) + 
  scale_fill_manual(values=CondColors$FB_Condition,
                    name = "Condition",
                    limits = c("EfaeHi_d1",
                               "EfaeLo_d1","EfaeLo_d7",
                               "EfaeMock_d8",
                               "EfaeTrain_d8"),
                    labels = c("Hi Dose - d1",
                               "Lo Dose- d1","EfaeLo - d7",
                               "Mock-Trained - d8",
                               "E.fae-Trained - d8"))
dev.off()



TPM_Mock <- filter(Gene_Plot_Reps, FB_Condition == "EfaeMock_d8") %>% .$Gene
TPM_Train <- filter(Gene_Plot_Reps, FB_Condition == "EfaeTrain_d8") %>% .$Gene

t.test(TPM_Mock, TPM_Train)

#### log-Fold Change Heatmap ####

## Use this to plot logFC heatmaps for particular gene lists

# Define the scale for your heatmap
LFCmyBreaks <- seq(from = -13.6, to = 13.6, by = c(0.1))
LFCpaletteLength <- length(LFCmyBreaks)
LFCmyColor <- colorRampPalette(c("blue", "white", "firebrick2"))(LFCpaletteLength)

# Make a dumb reference for naming conventions, so your heatmap annotates your conditions properly
FB_logFC_stupid_ref <- t(rbind(c("FB_EFAELO_d1", "FB_EFAEHI_d1", "FB_EFAELO_d7", "FB_EFAEMOCK_d8", "FB_EFAETRAIN_d8"),
                               c("FB_EFAELO_d1", "FB_EFAEHI_d1", "FB_EFAELO_d7", "FB_EFAEMOCK_d8", "FB_EFAETRAIN_d8")))
colnames(FB_logFC_stupid_ref) <- c("name","Condition")
FB_logFC_stupid_ref <-column_to_rownames(as.data.frame(FB_logFC_stupid_ref), "name")

# Filter your logFC matrix based on the lists of genes you're interested in

# IMD genes, as listed in Flybase Gene Group: IMD Signaling
IMD_testo <- FB_tru_logFC %>%
  dplyr::filter(.,flybaseID %in% IMD_list$`#SUBMITTED ID`)
IMD_testo <- IMD_testo %>%
  column_to_rownames(.,"symbol")%>%
  dplyr::select(.,starts_with("FB"))
colnames(IMD_testo) <- c("EfaeHi_d1","EfaeLo_d1","EfaeLo_d7","EfaeMock_d8","EfaeTrain_d8")
IMD_testo <- dplyr::select(IMD_testo, -c("EfaeLo_d7"))

# AMPs as defined by Flybase Gene Group: Antimicrobial Peptides
AMP_testo <- FB_tru_logFC %>%
  dplyr::filter(.,flybaseID %in% AMPs_list$FlyBaseID)
AMP_testo[11,"symbol"] <- "CecA2"
AMP_testo[12,"symbol"] <- "CecA1"
AMP_testo[6,"symbol"] <- "AttA"
AMP_testo[8,"symbol"] <- "AttC"
AMP_testo <- AMP_testo %>%
  remove_rownames(.) %>%
  column_to_rownames(.,"symbol")%>%
  dplyr::select(.,starts_with("FB"))
AMP_testo <- dplyr::select(AMP_testo, -c("EfaeLo_d7"))

# Previously annotated immune genes as used in Trohae, et al. 2018 & Ramirez-Corona, et al. 2021
immunospecific_testo <- FB_tru_logFC %>%
  dplyr::filter(.,flybaseID %in% FB_Immune_DE_TRAIN$.)
immunospecific_testo[19,"symbol"] <- "Acbp6"
immunospecific_testo <- immunospecific_testo %>%
  column_to_rownames(.,"symbol")%>%
  dplyr::select(.,starts_with("FB"))
colnames(immunospecific_testo) <- c("EfaeHi_d1","EfaeLo_d1","EfaeLo_d7","EfaeMock_d8","EfaeTrain_d8")
immunospecific_testo <- dplyr::select(immunospecific_testo, -c("EfaeLo_d7"))

# IMD effectors as defined by FlyBase Gene Group: IMD Signaling
IMD_effector_testo <- FB_tru_logFC %>%
  dplyr::filter(.,flybaseID %in% c("FBgn0004240","FBgn0034407","FBgn0012042","FBgn0041581",
                                   "FBgn0041579","FBgn0038530", "FBgn0010388", "FBgn0010385",
                                   "FBgn0000276","FBgn0000277"))
IMD_effector_testo[5,"symbol"] <- "AttA"
IMD_effector_testo[7,"symbol"] <- "AttC"
IMD_effector_testo[10,"symbol"] <- "CecA2"
IMD_effector_testo <- IMD_effector_testo[order(IMD_effector_testo$symbol),]
rownames(IMD_effector_testo) <- NULL
IMD_effector_testo <- IMD_effector_testo %>%
  column_to_rownames(.,"symbol") %>%
  dplyr::select(.,starts_with("FB"))
IMD_effector_testo <- dplyr::select(IMD_effector_testo, -c("FB_EFAELO_d7"))

# Core response genes from the Troha, et al. 2018 paper
CrebA_testo <- FB_tru_logFC %>%
  dplyr::filter(.,flybaseID %in% immune_CrebAPapeCore_list$FlybaseID)
rownames(CrebA_testo) <- NULL
CrebA_testo <- CrebA_testo %>%
  column_to_rownames(.,"symbol") %>%
  dplyr::select(.,starts_with("FB"))
CrebA_testo <- dplyr::select(CrebA_testo, -c("FB_EFAELO_d7"))

# Toll genes
Toll_testo <-FB_tru_logFC %>%
  dplyr::filter(.,flybaseID %in% TollPathway_list$FlyBaseID) %>%
  column_to_rownames(.,"symbol") %>%
  dplyr::select(.,starts_with("FB"))
colnames(Toll_testo) <- c("EfaeHi_d1","EfaeLo_d1","EfaeLo_d7","EfaeMock_d8","EfaeTrain_d8")
Toll_testo <- dplyr::select(Toll_testo, -c("EfaeLo_d7"))


png(filename = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/Figures/FB_CrebA_Hmap.png",width = 1100, height=1700, res=300)
pheatmap(as.matrix(CrebA_testo),
         color=LFCmyColor, 
         breaks=LFCmyBreaks,
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = FB_logFC_stupid_ref,
         annotation_colors = FB_logFCColors,
         legend = T,
         show_colnames = F,
         show_rownames = F,
         gaps_col = c(3),
         annotation_legend = F,
         scale = "none")
dev.off()




#### LOITERING TEST ####

## Filter based on conditionals that define immune loitering genes

### Finds and removes genes that have are not turned on relative to ctrl in training-d1, training-d7, primed-d8 
# 2,895 genes
FB_loiter_test <- transform(FB_tru_logFC, `FB_EFAELO_d7` = 
                                ifelse( `FB_EFAELO_d7` > 0, # Gene levels at d7 have to be more than d0
                                        `FB_EFAELO_d7`, NA)) %>% 
  transform(., `FB_EFAELO_d1` = 
              ifelse( `FB_EFAELO_d1` > 0, # Gene levels have to go up at d1 compared to d0
                      `FB_EFAELO_d1`, NA)) %>%
  transform(., `FB_EFAETRAIN_d8` = 
              ifelse( `FB_EFAETRAIN_d8` > 0, # Gene levels have to go up at d8 compared to d0
                      `FB_EFAETRAIN_d8`, NA)) %>%
  drop_na(.)

### Creates a dataframe for the conditional that genes have to significantly go up from ctrl to training-d1
# 1,961 genes
FB_loitering_cond_1 <- FBRes_EFAELO_d1vsCTRL_d0$table %>%
  rownames_to_column('gene') %>%
  dplyr::filter(logFC > 0 & FDR < 0.05, ) %>%
  column_to_rownames('gene')

FB_loitering_cond_1$flybaseID <- row.names(FB_loitering_cond_1)

### Creates a dataframe for the conditional that genes have to SIGNIFICANTLY higher at d8 than d0
#5,000 genes
FB_loitering_cond_2 <- FBRes_EFAETRAIN_d8vsCTRL_d7$table %>%
  rownames_to_column('gene') %>%
  dplyr::filter(FDR < 0.05, ) %>%
  column_to_rownames('gene')

FB_loitering_cond_2$flybaseID <- row.names(FB_loitering_cond_2)

### Creates a dataframe for the conditional that genes have to either stay at the same level or go up from training-d7 to primed-d8
#2,906 genes
FB_loitering_cond_3 <- FBRes_EFAETRAIN_d8vsEFAELO_d7$table %>%
  rownames_to_column('gene') %>%
  dplyr::filter(logFC>=0, ) %>%
  column_to_rownames('gene')

FB_loitering_cond_3$flybaseID <- row.names(FB_loitering_cond_3)


### Number of Loitering Genes: 33
FB_loiter_genes <- FB_loiter_test %>%
  dplyr::filter(., .$flybaseID %in% rownames(FB_loitering_cond_1)) %>%
  dplyr::filter(., .$flybaseID %in% rownames(FB_loitering_cond_2)) %>%
  dplyr::filter(., .$flybaseID %in% rownames(FB_loitering_cond_3))

readr::write_tsv(FB_loiter_genes, "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/FB_loiter_genes.tsv")

#### RECALL RESPONSE TEST ####

## Filter based on conditionals that define a potentiated recall response

# logFC has to be higher at Trained-d8 vs Training-d1
# 5480 genes
FB_twofast_cond_1 <- dplyr::filter(FB_tru_logFC, FB_EFAETRAIN_d8 > FB_EFAELO_d1) %>%
  tibble::remove_rownames() %>%
  column_to_rownames("flybaseID")

# Genes need to be significantly turned on upon re-infection
# 862 genes
FB_twofast_cond_2 <- FBRes_EFAETRAIN_d8vsCTRL_d7$table %>%
  dplyr::filter(logFC > 0 & FDR < 0.1)

## Potentially remove FDR threshold
# Genes need to be turned on after the initial infection
# 3803 genes
FB_twofast_cond_3 <- dplyr::filter(FB_tru_logFC, FB_EFAELO_d1 > 0.5) %>%
  tibble::remove_rownames() %>%
  column_to_rownames("flybaseID")

# Genes are returned to age-control levels or lower right before re-infection
# 6784 genes
FB_twofast_cond_4 <- dplyr::filter(FB_tru_logFC, FB_EFAELO_d7 <= 0) %>%
  tibble::remove_rownames() %>%
  column_to_rownames("flybaseID")

## Use this to filter out any genes that are significantly up-regulated in Mock-primed samples
FB_twofast_cond_5 <- FBRes_EFAEMOCK_d8vsCTRL_d7$table %>%
  dplyr::filter(FDR < 0.05)

# Number of Faster Genes: 161
FB_faster_genes <- FB_tru_logFC %>%
  dplyr::filter(., .$flybaseID %in% rownames(FB_twofast_cond_1)) %>%
  dplyr::filter(., .$flybaseID %in% rownames(FB_twofast_cond_2)) %>%
  dplyr::filter(., .$flybaseID %in% rownames(FB_twofast_cond_3)) %>%
  dplyr::filter(., .$flybaseID %in% rownames(FB_twofast_cond_4)) %>%
  dplyr::filter(., !.$flybaseID %in% rownames(FB_twofast_cond_5))

readr::write_tsv(FB_faster_genes, path = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/FB_faster_genes.tsv")

#### SPECIFIC RESPONSE TEST ####

## Pick out what's significantly up-regulated in each 24hr condition 
FB_sig_TRAINvsCTRL <- FBRes_EFAETRAIN_d8vsCTRL_d7$table %>%
  rownames_to_column('gene') %>%
  dplyr::filter(`FDR`<0.05 & `logFC`>1, ) %>%
  column_to_rownames('gene')

FB_sig_LOd1vsCTRL <- FBRes_EFAELO_d1vsCTRL_d0$table %>%
  rownames_to_column('gene') %>%
  dplyr::filter(FDR<0.05 & `logFC`>1, ) %>%
  column_to_rownames('gene')

FB_sig_HId1vsCTRL <- FBRes_EFAEHI_d1vsCTRL_d0$table %>%
  rownames_to_column('gene') %>%
  dplyr::filter(FDR<0.05 & `logFC`>1, ) %>%
  column_to_rownames('gene')

FB_sig_MOCKvsCTRL <- FBRes_EFAEMOCK_d8vsCTRL_d7$table %>%
  rownames_to_column('gene') %>%
  dplyr::filter(FDR<0.05 & `logFC`>1, ) %>%
  column_to_rownames('gene')

FB_sig_MOCKvsEFAEPRIMED <- FBRes_EFAETRAIN_d8vsEFAEMOCK_d8$table %>%
  rownames_to_column('gene') %>%
  dplyr::filter(FDR<0.05 & `logFC`>1, ) %>%
  column_to_rownames('gene')
write_tsv(FB_sig_MOCKvsEFAEPRIMED, file = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/FB_Mock_vs_EfaePrimed.tsv")

FB_sig_EFAEPRIMEDvsMOCK <- FBRes_EFAETRAIN_d8vsEFAEMOCK_d8$table %>%
  rownames_to_column('gene') %>%
  dplyr::filter(FDR<0.05 & `logFC`<1, ) %>%
  column_to_rownames('gene')
write_tsv(FB_sig_EFAEPRIMEDvsMOCK, file = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/FB_EfaePrimed_vs_Mock.tsv")


## Make a venn diagram of up-regulated genes #
FB_vennlist <- list('E.fae-Trained'=unique(rownames(FB_sig_TRAINvsCTRL)), 
                      'Low Dose-d1'=unique(rownames(FB_sig_LOd1vsCTRL)),
                      'Hi Dose-d1'=unique(rownames(FB_sig_HId1vsCTRL)),
                      'Mock-Trained'= unique(rownames(FB_sig_MOCKvsCTRL)))

venn.diagram(FB_vennlist,
             filename = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/Figures/FB_VennUp.png",
             output = TRUE,
             fill = c("#1cd9b6", "#a7f542", "#4b6e1e", "#355a78"),
             cat.cex = 0,
             cex = 0.3,
             height = 3000 , 
             width = 5000 ,
             cat.fontfamily = "sans",
             cat.pos = c(-5, 5, 20, 20),
             fontfamily = "sans")

## Get the list of genes in each part of the venn diagram circle
FB_overlap <- calculate.overlap(FB_vennlist)

write_tsv(as_tibble(FB_overlap$a9), path="/Users/cabweebwub/Desktop/Efae RNA-seq Paper/FB_Down_SpecEfaeTrained.tsv")
write_tsv(as_tibble(FB_overlap$a3), path="/Users/cabweebwub/Desktop/Efae RNA-seq Paper/FB_Down_SpecMock.tsv")
write_tsv(as_tibble(FB_overlap$a1), path="/Users/cabweebwub/Desktop/Efae RNA-seq Paper/FB_Down_SpecHiDose.tsv")
write_tsv(as_tibble(FB_overlap$a14), path="/Users/cabweebwub/Desktop/Efae RNA-seq Paper/FB_Down_SpecLowDose.tsv")
write_tsv(as_tibble(FB_overlap$a6), path="/Users/cabweebwub/Desktop/Efae RNA-seq Paper/FB_Down_Common.tsv")
write_tsv(as_tibble(FB_overlap$a10), path="/Users/cabweebwub/Desktop/Efae RNA-seq Paper/FB_Up_EfaePrimed_MockPrimed_Shared.tsv")


#### HEATMAP PLOTTING ####
# Creating data frames for heatmaps based on differential, loitering, and recall analysis

## Loitering heatmap
FB_loiter_hmap <- dplyr::select(FB_loiter_genes, c('FB_EFAELO_d1',
                                                     'FB_EFAELO_d7',
                                                     'FB_EFAETRAIN_d8',
                                                     'symbol'))

# Renaming gene symbols that weren't properly annotated in the biomart download
FB_loiter_hmap[3,"symbol"] <- "BomBc3"
FB_loiter_hmap[4,"symbol"] <- "GNBP-like 3"
FB_loiter_hmap[5,"symbol"] <- "Bbd"
FB_loiter_hmap[6,"symbol"] <- "BomS1"
FB_loiter_hmap[7,"symbol"] <- "BomS2"
FB_loiter_hmap[8,"symbol"] <- "BomS3"
FB_loiter_hmap[9,"symbol"] <- "BomS5"
FB_loiter_hmap[10,"symbol"] <- "BomBc1"
FB_loiter_hmap[11,"symbol"] <- "Dso1"
FB_loiter_hmap[12,"symbol"] <- "BaraA1"
FB_loiter_hmap[13,"symbol"] <- "lncRNA:CR33942"
FB_loiter_hmap[14,"symbol"] <- "BaraA2"

FB_loiter_hmap <- tibble::remove_rownames(FB_loiter_hmap) %>% 
  tibble::column_to_rownames('symbol')


## Recall response heatmap
FB_faster_hmap <- dplyr::select(FB_faster_genes, c('FB_EFAELO_d1',
                                                   'FB_EFAELO_d7',
                                                   'FB_EFAETRAIN_d8',
                                                   'symbol')) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames('symbol')


## Priming-specific heatmap
FB_priming_hmap <- dplyr::filter(FB_tru_logFC, `flybaseID` %in% unlist(FB_overlap))
FB_priming_hmap <- FB_priming_hmap %>%
  arrange(factor(flybaseID, levels = c(FB_overlap$a9, FB_overlap$a3, FB_overlap$a1, FB_overlap$a14, FB_overlap$a4, FB_overlap$a10,
                                       FB_overlap$a15, FB_overlap$a13, FB_overlap$a8, FB_overlap$a2, FB_overlap$a5, FB_overlap$a7,
                                       FB_overlap$a12, FB_overlap$a11, FB_overlap$a6)))

FB_priming_hmap <- dplyr::select(FB_priming_hmap, c('FB_EFAEHI_d1',
                                                    'FB_EFAELO_d1',
                                                    'FB_EFAEMOCK_d8',
                                                    'FB_EFAETRAIN_d8',
                                                    'symbol')) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames('symbol')

## Making the color bar for the heatmaps #
# How many breaks you want in your color scale (usually change it so it matches the max and min of my data)
FB_loiter_Breaks <- seq(from = 0, to = 11.5, by = c(0.1))
FB_faster_Breaks <- seq(from = -4.2, to = 4.2, by = c(0.1))
FB_priming_Breaks <- seq(from = -9.8, to = 9.8, by = c(0.2))
FB_IMDHmap_Breaks <- seq(from = -6.3, to = 6.3, by = c(0.1))
FB_top50_Breaks <- seq(from = -1.5, to = 1.5, by = c(0.1))

# Check out how many values are in myBreaks and change the paletteLength
FB_loiter_paletteLength <- length(FB_loiter_Breaks)
FB_faster_paletteLength <- length(FB_faster_Breaks)
FB_priming_paletteLength <- length(FB_priming_Breaks)
FB_IMDHmap_paletteLength <- length(FB_IMDHmap_Breaks)
FB_top50_paletteLength <- length(FB_top50_Breaks)

# Choose what colors you want as your lower, middle, and upper cutoffs
FB_loiter_Color <- colorRampPalette(c("white", "firebrick2"))(FB_loiter_paletteLength)
FB_faster_Color <- colorRampPalette(c("navyblue","white", "firebrick2"))(FB_faster_paletteLength)
FB_priming_Color <- colorRampPalette(c("navyblue","white", "firebrick2"))(FB_priming_paletteLength)
FB_IMDHmap_Color <- colorRampPalette(c("navyblue","white", "firebrick2"))(FB_IMDHmap_paletteLength)
FB_top50_Color <- colorRampPalette(c("blue","white", "firebrick2"))(FB_top50_paletteLength)

## Plotting the heatmaps #
png(filename = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/Figures/FB_Loiter_Heatmap.png",width = 1400, height=1200, res=300)
pheatmap(FB_loiter_hmap, cluster_cols = F, 
         cluster_rows = T, 
         color=FB_loiter_Color, breaks=FB_loiter_Breaks,
         #gaps_col = c(3),
         #gaps_row = c(149,557,566,596,779,822),
         #gaps_row = c(388,1119,1121,1122,1299),
         cex=1,
         annotation_col = FB_logFC_stupid_ref,
         annotation_colors = FB_logFCColors,
         show_rownames = T,
         show_colnames = F,
         annotation_legend = F,
         treeheight_row = 0, treeheight_col = 0)
dev.off()  

#### FB vs Hemocyte Testing ####
Both_specificDE_TRAIN <- intersect(unique(FB_specificDE_TRAIN$.), unique(Hemo_specificDE_TRAIN$.))
Both_faster <- intersect(unique(FB_faster_genes$flybaseID), unique(Hemo_faster_genes$flybaseID))
Both_loiter <- intersect(unique(FB_loiter_genes$flybaseID), unique(Hemo_loiter_genes$flybaseID))
