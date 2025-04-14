library(DESeq2)
library(ggplot2)
library(ggrepel)

# Define sample metadata
mysamples = c("SRR31781635",
             "SRR31781636",
             "SRR31781640",
             "SRR31781641")
mygeno = rep(c("spn1", "WT"), each = 2)

sample_info <- data.frame(
  sample = mysamples,
  genotype = mygeno
)
rownames(sample_info) <- sample_info$sample

# Load counts
counts <- read.table("results/counts/count_matrix.txt", header=TRUE, row.names=1, comment.char="#")
counts <- counts[, grepl(paste(sample_info$sample,collapse = "|"),
                               colnames(counts))]
colnames(counts) <- mysamples

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ genotype)
## Remove low-expression genes
remove_flag = FALSE
if(remove_flag){
  dds <- dds[rowSums(counts(dds))>= 10, ]
}

dds <- DESeq(dds)

res <- results(dds, contrast=c("genotype", "spn1", "WT"))
write.csv(as.data.frame(res), file="results/counts/deseq2_results.csv")
summary(res)
res_sig <- res[which(res$padj < 0.05 & res$log2FoldChange< -1),]
as_tibble(res_sig)

# Diagnostics
plotDispEsts(dds)




# Plot volcano
res$padj[is.na(res$padj)] <- 1


# Genes highlighted in the Spt6-Spn1 JBC study
study_genes <- c(
  "SNR47",
  "SNR13",
  "TDH3",
  "PDC1",
  "RPS4A",
  "GCN4",
  "HSP82"
)

genes_of_interest <- c(
  "YGR053G",  # non-transcribed gene, nucleosome positioning control
  "YGL202W",  # actively transcribed gene, altered nucleosome pattern
  "SNR13",    # shows up as significantly regulated in your DESeq2 output
  "HTA1",     # histone protein
  "HTA2",     # histone variant
  "HTB1",     # histone protein
  "HTB2",     # histone variant
  "SPT6",     # key chaperone
  "SPN1",     # partner in complex
  "RPB1",     # RNAPII subunit
  "SPT5",     # elongation factor
  "CDC73",    # transcription complex
  "PAF1",     # elongation complex subunit
  "FACT",     # histone chaperone complex
  "H2B",      # histone modulated by this pathway
  "CKII"      # involved in phosphorylation affecting Spt6-Spn1
)


## My better plot
(p1 <- res %>% 
    as_tibble(rownames = "Gene") %>% 
    mutate(Gene = str_to_upper(Gene)) %>% 
    mutate(Significant = ifelse(padj < 0.05,
                                "Significant",
                                "NS"),
           Label = ifelse(-log10(padj)> 50 | log2FoldChange>10 | Gene %in% genes_of_interest, 
                          Gene,
                          NA)
    ) %>% 
    ggplot(aes(log2FoldChange, -log10(padj),
               fill = Significant,
               label=Label)) +
    geom_point(alpha = .7, pch = 21, col = 'gray10') +
    ggtitle("Volcano plot: spn1-R263D vs WT")+
    geom_label_repel() +
    xlab(expression(Log[2]~Fold~Change)) +
    ylab(expression(-Log[10]~p-value)) +    
    scale_fill_viridis_d()
  # scale_fill_manual(values = c("Significant" = "")) +
)

key_genes <- c(
  "HSP82", # stress response
  "SPT6",  # chromatin remodelling
  "RNR1"  # Cell cycle/replication
)



res %>% 
  as_tibble(rownames = 'Gene') %>% 
  mutate(Gene = str_to_upper(Gene)) %>% 
  filter(Gene %in% study_genes)


# Heatmap -----------------------------------------------------------------
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

# Load or reuse your DESeq2 object
vsd <- vst(dds, blind = TRUE)  # variance-stabilized data

# Choose top DE genes (e.g. top 30 most significant)
res_ordered <- res[order(res$padj), ]
top_genes <- head(rownames(res_ordered), 30)

# Extract expression matrix for top genes
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)  # center each gene's row

# Create sample annotation from colData
ann <- as.data.frame(colData(dds)[, "genotype", drop=FALSE])

# Plot
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = ann,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
         fontsize_row = 8,
         main = "Top 30 Differentially Expressed Genes (spn1 vs WT)")



