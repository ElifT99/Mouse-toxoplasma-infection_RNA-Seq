# Overrepresentation analysis

#Load packages
library(DESeq2)
library(data.table) #for clusterProfiler
library(clusterProfiler)
library(org.Mm.eg.db) #mouse organism annotation package
library(enrichplot)
library(DOSE) #


#DESeq2 results table for one contrast (Disease.WT and Control.WT)
res_WT <- results(dds_interaction,
                  contrast = c("group", "Disease.WT", "Control.WT"))

#Mnsembl IDs without versions for GO
#all tested genes = universe (full background set of genes)
universe_MM <- rownames(res_WT)
universe_MM_go <- sub("\\.\\d+$", "", universe_MM)

# DE genes for this contrast
de_WT <- res_WT[which(res_WT$padj < 0.05 & !is.na(res_WT$padj)), ]
gene_WT <- rownames(de_WT)
gene_WT_go <- sub("\\.\\d+$", "", gene_WT)

#EnrichGO with the stripped Ensembl IDs
ego_WT_BP <- enrichGO(
  gene          = gene_WT_go,
  universe      = universe_MM_go,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",#Biological Process GO terms
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_BP_df <- as.data.frame(ego_WT_BP)
head(ego_BP_df)

#General Dotplot
ego_BP_top <- ego_BP_df |>
  dplyr::arrange(p.adjust) |>
  dplyr::slice_head(n = 20)

ggplot(
  ego_BP_top,
  aes(x = -log10(p.adjust),
      y = reorder(Description, p.adjust))
) +
  geom_point(aes(size = Count)) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "GO Biological Process enrichment of DE genes: WT infected vs WT control",
    subtitle = "Over-representation analysis (ORA) using enrichGO (clusterProfiler), Mus musculus, BP ontology",
    x = "-log10(adj. p-value)",
    y = "GO Biological Process term",
    size = "Gene count",
  ) +
  theme_bw(base_size = 10)
#shows very strong enrichment for immune‑cell movement and muscle‑related processes

#Closer look at immune response related GO terms
ego_BP_df <- as.data.frame(ego_BP)

#defining immune GO terms
immune_go <- c(
  "GO:0045087","GO:0006955","GO:0002682","GO:0045321","GO:0046651",
  "GO:0006954","GO:0001816","GO:0043123","GO:0008009",
  "GO:0002227","GO:0009405","GO:0044409","GO:0044419",
  "GO:0006950","GO:0006457","GO:0008152"
)

ego_imm <- ego_BP_df[ego_BP_df$ID %in% immune_go, ]#filtering

ggplot(
  ego_imm,
  aes(
    x = -log10(p.adjust),
    y = reorder(Description, p.adjust)
  )
) +
  geom_point(aes(size = Count)) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "Immune-related GO Biological Processes enriched in DE genes\nWT infected vs WT control (Mus musculus)",
    x = "-log10(adj. p-value)",
    y = "Immune-related GO term",
    size = "Gene count",
    color = "Significance\n(-log10 adj. p)"
  ) +
  theme_bw(base_size = 10)


#2 DESeq2 contrast: Infection.Mutant vs Infection.WT

# DESeq2 results for Disease.Ifnar-/-x Ifngr-/-, Disease.WT
res_mut_inf <- results(
  dds_interaction,
  contrast = c("group", "Disease.Ifnar-/-x Ifngr-/-", "Disease.WT")
)

#background (all tested genes)
universe_MM     <- rownames(res_mut_inf)
universe_MM_go  <- sub("\\.\\d+$", "", universe_MM)

# DE genes for this contrast
de_mut_inf  <- res_mut_inf[which(res_mut_inf$padj < 0.05 & !is.na(res_mut_inf$padj)), ]
gene_MM     <- rownames(de_mut_inf)
gene_MM_go  <- sub("\\.\\d+$", "", gene_MM)

#GO enrichment for this contrast
ego_mut_BP <- enrichGO(
  gene          = gene_MM_go,
  universe      = universe_MM_go,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_BP_mut_df <- as.data.frame(ego_mut_BP)

ego_BP_top <- ego_BP_mut_df |>
  arrange(p.adjust) |>
  slice_head(n = 20)

ggplot(
  ego_BP_top,
  aes(x = -log10(p.adjust),
      y = reorder(Description, p.adjust))
) +
  geom_point(aes(size = Count)) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title    = "GO Biological Process enrichment of DE genes:\nInfection.Mutant vs Infection.WT",
    subtitle = "Over-representation analysis (ORA) using enrichGO (clusterProfiler), Mus musculus, BP ontology",
    x        = "-log10(adj. p-value)",
    y        = "GO Biological Process term",
    size     = "Gene count",
    color    = "Significance\n(-log10 adj. p)"
  ) +
  theme_bw(base_size = 10)

#Immune response GO terms filtering for the second contrast
ego_imm_mut <- ego_BP_mut_df[ego_BP_mut_df$ID %in% immune_go, ]

ego_imm_auto <- ego_BP_mut_df[
  grepl("immune|leukocyte|lymphocyte|chemokine|cytokine|inflamm|NF-kappaB|defense",
        ego_BP_mut_df$Description,
        ignore.case = TRUE),
]

ggplot(
  ego_imm_mut,
  aes(
    x = -log10(p.adjust),
    y = reorder(Description, p.adjust)
  )
) +
  geom_point(aes(size = Count)) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "Immune-related GO Biological Processes enriched in DE genes\n KO infected vs WT infected (Mus musculus)",
    x = "-log10(adj. p-value)",
    y = "Immune-related GO term",
    size = "Gene count",
    color = "Significance\n(-log10 adj. p)"
  ) +
  theme_bw(base_size = 10)

## Visualization through gene-concept network of Infected.WT and Control.WT
# Extract DE gene list

# Clean results
res_clean_WT$log2FoldChange[is.na(res_clean_WT$log2FoldChange)] <- 0
res_clean_WT$padj[is.na(res_clean_WT$padj)] <- 1

# Named gene list: rows as SYMBOL (from rowData), values as log2FC
# Positive = up in Disease.WT (infected)
geneList_WT <- res_clean_WT$log2FoldChange
names(geneList_WT) <- rowData(dds_interaction)$symbol[rownames(res_clean_WT)]
geneList_WT <- sort(geneList_WT, decreasing = TRUE)  # For GSEA later

# Significant DEGs (|log2FC| > 2 & padj < 0.01)
sig_genes_WT_up <- rownames(res_clean_WT)[res_clean_WT$padj < 0.01 & res_clean_WT$log2FoldChange > 2]
sig_genes_WT_down <- rownames(res_clean_WT)[res_clean_WT$padj < 0.01 & res_clean_WT$log2FoldChange < -2]
sig_genes_WT <- c(sig_genes_WT_up, sig_genes_WT_down)

# Convert to ENTREZ IDs
#For sig_genes_WT (as SYMBOL vector)
sig_symbols_WT <- rowData(dds_interaction)$symbol[sig_genes_WT]
sig_entrez_WT <- bitr(sig_symbols_WT, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
sig_entrez_WT <- sig_entrez_WT[!duplicated(sig_entrez_WT)]

#ORA Enrichment
universe <- rowData(dds_interaction)$symbol[!is.na(rowData(dds_interaction)$symbol)]  # Background

ego_WT <- enrichGO(
  gene = sig_entrez_WT,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # BP/MF/CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  universe = bitr(universe, "SYMBOL", "ENTREZID", OrgDb=org.Mm.eg.db)$ENTREZID,
  readable = TRUE  # Genes as SYMBOL
)
head(ego_WT)
