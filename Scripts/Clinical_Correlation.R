
#### Raw RMA non log2 Transformed ####

setwd("~/Desktop/Data array/Arrayer/RIKEN")
packages <- c("corrplot", "readxl", "Hmisc", "seriation", "magrittr", "vroom" )
invisible(lapply(packages, library, character.only = TRUE))
data <- read_xlsx("eset_all_20101119.xlsx")
data_fil <- data[,5:62]
data_fil$ID <- NULL

genes <- c("GLS", "GLUL")
data_GLS_GLUL <- data_fil[data_fil$symbol %in% genes,]
data_GLS_GLUL <- data_GLS_GLUL[-2,]
anno_GLS_GLUL <- data_GLS_GLUL$symbol
data_GLS_GLUL$symbol <- NULL

library(data.table)

data_GLS_GLUL <- transpose(data_GLS_GLUL)
colnames(data_GLS_GLUL) <- anno_GLS_GLUL
data_GLS_GLUL$Qratio <- data_GLS_GLUL$GLS/data_GLS_GLUL$GLUL
data_GLS_GLUL$Qratio <- log2(data_GLS_GLUL$Qratio)

### Log2 transformed data ###

data <- read_xlsx("eset_all_20101119.xlsx", sheet = 4)
data$symbol <- str_replace_all(data$symbol,"---" ," ")
data <- data[which(data$symbol != " "),]
data_fil <- data[,5:62]
data_fil$ID <- NULL
anno <- data_fil$symbol
data_fil$symbol <- NULL

data_transposed <- transpose(data_fil)
colnames(data_transposed) <- anno
data_transposed$Qratio <- data_GLS_GLUL$Qratio
data_transposed$Qratio
rm <- rcorr(as.matrix(data_transposed), type = "spearman")

Corr <- as.data.frame(rm$r)
Pval <- as.data.frame(rm$P)
genes <- rownames(Corr)
genes <- as.data.frame(genes)

Qratio_Corr <- as.data.frame(Corr$Qratio)
Qratio_Pval <- as.data.frame(Pval$Qratio)
Qratio_Pval_Cor <- cbind(Qratio_Corr, Qratio_Pval, genes)
Qratio_Pval_Cor <- Qratio_Pval_Cor[-32374,]

library(dplyr)

Qratio_Pval_Cor_Pval_Cor_Ranked <- Qratio_Pval_Cor %>% mutate(genes = factor(genes, levels = genes[order(Qratio_Pval_Cor$`Pval$Qratio`, Qratio_Pval_Cor$`Corr$Qratio`)]))

colnames(Qratio_Pval_Cor_Pval_Cor_Ranked) <- c("r", "Pval", "genes")

Qratio_Pval_Cor_Pval_Cor_Ranked$adjPval <- p.adjust(Qratio_Pval_Cor_Pval_Cor_Ranked$Pval, method = "BH")

Qratio_fil <- Qratio_Pval_Cor_Pval_Cor_Ranked[which(Qratio_Pval_Cor_Pval_Cor_Ranked$adjPval < 0.05),]
Qratio_fil <- Qratio_fil[Qratio_fil$genes %in% data$symbol[!data$symbol== "---"],]

Qratio_Pos_cor <- Qratio_fil[which(Qratio_fil$r > 0.3),]
Qratio_Neg_cor <- Qratio_fil[which(Qratio_fil$r < -0.3),]

write.table(Qratio_Pos_cor, file = "Qratio_Pos_Corr.txt", row.names = F, col.names = T, sep = "\t")
write.table(Qratio_Neg_cor, file = "Qratio_Neg_Corr.txt", row.names = F, col.names = T, sep = "\t")

Pos <- read.delim("Qratio_Pos_Corr.txt")
Neg <- read.delim("Qratio_Neg_Corr.txt")

list_Pos <- Pos$genes
list_Neg <- Neg$genes

library(enrichR)

databases <- c("KEGG_2021_Human", "GO_Biological_Process_2021")
enriched_Pos <- enrichr(list_Pos, databases)
enriched_Neg <- enrichr(list_Neg, databases)

KEGG_Neg <- enriched_Neg$KEGG_2021_Human
GO_Neg <- enriched_Neg$GO_Biological_Process_2021

KEGG_Pos <- enriched_Pos$KEGG_2021_Human
GO_Pos <- enriched_Pos$GO_Biological_Process_2021

write.table(KEGG_Neg, file = "KEGG_Neg_tab.txt", row.names = F, col.names = T, sep = "\t")
write.table(GO_Neg, file = "GO_Neg_Tab.txt", row.names = F, col.names = T, sep = "\t")
write.table(KEGG_Pos, file = "KEGG_Pos_tab.txt", row.names = F, col.names = T, sep = "\t")
write.table(GO_Pos, file = "GO_Pos_Tab.txt", row.names = F, col.names = T, sep = "\t")

plotEnrich(enriched_Pos$KEGG_2021_Human, showTerms = 30)
plotEnrich(enriched_Pos$GO_Biological_Process_2021, showTerms = 30)
plotEnrich(enriched_Neg$KEGG_2021_Human, showTerms = 30)
plotEnrich(enriched_Neg$GO_Biological_Process_2021, showTerms = 30)

library(rrvgo)
library(org.Hs.eg.db)

GO_Neg <- enriched_Neg$GO_Biological_Process_2021
GO_Neg <- substr(x = GO_Neg$Term,start = nchar(GO_Neg$Term)-10,stop = nchar(GO_Neg$Term)-1)

GO_Pos <- enriched_Pos$GO_Biological_Process_2021
GO_Pos <- substr(x = GO_Pos$Term,start = nchar(GO_Pos$Term)-10,stop = nchar(GO_Pos$Term)-1)

simMatrix_Pos <- calculateSimMatrix(GO_Pos,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

simMatrix_Neg <- calculateSimMatrix(GO_Neg,
                                    orgdb="org.Hs.eg.db",
                                    ont="BP",
                                    method="Rel")

scores_Pos <- setNames(-log10(enriched_Pos$GO_Biological_Process_2021$Adjusted.P.value), GO_Pos)
reducedTerms_Pos <- reduceSimMatrix(simMatrix_Pos,
                                scores_Pos,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

scores_Neg <- setNames(-log10(enriched_Neg$GO_Biological_Process_2021$Adjusted.P.value), GO_Neg)
reducedTerms_Neg <- reduceSimMatrix(simMatrix_Neg,
                                    scores_Neg,
                                    threshold=0.7,
                                    orgdb="org.Hs.eg.db")

rrvgo::shiny_rrvgo()

GO_Neg <- GO_Neg[1:30]
GO_Pos <- GO_Pos[1:30]

scores_Neg <- scores_Neg[1:30]
scores_Pos <- scores_Pos[1:30]

HM_Pos <- heatmapPlot(simMatrix_Pos,
            reducedTerms_Pos,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

HM_Neg <- heatmapPlot(simMatrix_Neg,
                      reducedTerms_Neg,
                      annotateParent=TRUE,
                      annotationLabel="parentTerm",
                      fontsize=6)

Scatter_Pos <- scatterPlot(simMatrix_Pos, reducedTerms_Pos, addLabel = T)
Scatter_Neg <- scatterPlot(simMatrix_Neg, reducedTerms_Neg, addLabel = T)

