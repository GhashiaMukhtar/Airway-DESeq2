# Phase 1 — Load the airway data
library(DESeq2)
library(airway)
data(airway)
airway
colData(airway)
head(assay(airway))
dim(assay(airway))


# Phase 2 — Build the DESeq2 object & quality check
dds <- DESeqDataSet(airway, design = ~dex)     # Build DESeq2 obj]
dds
dds$dex <- relevel(dds$dex, ref = 'untrt')     # Set untreated as reference
keep <- rowSums(counts(dds)) >= 10      # Filter the row with less then 10 reads
dds <- dds[keep, ]
dds
table(colData(dds)$dex)

# PCA Quality Check
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "dex")
library(ggplot2)
ggsave("figures/PCA_check.png", width = 9, height = 10)
# Phase 3 — Differential expression 
dds <- DESeq(dds)
res <- results(dds)
res
summary(res)

resOrdered <- res[order(res$padj), ]
head(resOrdered, 10)

resOrdered$symbol <- rowData(airway)$symbol[match(rownames(resOrdered), rownames(airway))]
head(resOrdered, 10)

write.csv(as.data.frame(resOrdered), "results/deseq2_results.csv")  # Save the result


# Phase 4 — The volcano plot 
library(EnhancedVolcano)
EnhancedVolcano(resOrdered,
                lab = resOrdered$symbol,
                x = "log2FoldChange",
                y = "padj",
                title = 'Dexamethasone treated vs untreated',
                pCutoff = 0.05,
                FCcutoff = 1.00)
ggsave("figures/volcano_plot.png", width = 9, height = 7)  # Save the plot


# Phase 5 — Pathway enrichment
# Upregulated genes
library(clusterProfiler)
library(org.Hs.eg.db)

sig_genes_up <- rownames(resOrdered)[which(resOrdered$padj < 0.05 & resOrdered$log2FoldChange > 1.0)]
length(sig_genes_up)
head(sig_genes_up, 10)

ego_up <- enrichGO(gene = sig_genes_up,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENSEMBL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)
head(ego_up, 10)
library(ggplot2)
dotplot(ego_up, showCategory = 15)
ggsave("figures/pathway_enrichment_upregulated.png", width = 9, height = 8)


library(clusterProfiler)
library(org.Hs.eg.db)

sig_genes_down <- rownames(resOrdered)[which(resOrdered$padj < 0.05 & resOrdered$log2FoldChange < -1.0)]
length(sig_genes_down)
head(sig_genes_down, 10)

ego_down <- enrichGO(gene = sig_genes_down,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENSEMBL",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05)
head(ego_down, 10)
library(ggplot2)
dotplot(ego_down, showCategory = 15)
ggsave("figures/pathway_enrichment_downregulated.png", width = 9, height = 8)

as.data.frame(ego_up)[1:15, c("ID","Description","GeneRatio","p.adjust","Count")]
as.data.frame(ego_down)[1:15, c("ID","Description","GeneRatio","p.adjust","Count")]
sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)
sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE)
