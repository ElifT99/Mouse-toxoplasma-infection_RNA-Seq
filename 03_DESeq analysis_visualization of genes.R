#DESeq analysis: possible visualizations

#Load packages
library(DESeq2)
library(dplyr)
library(ggplot2)
library(tibble)
library(ggrepel)
library(ggbreak)

#Vulcano plot with top differentially expressed genes
alpha   <- 0.01
lfc_thr <- 2

volc_WT <- as.data.frame(res_infection_WT_vs_control) %>%
  rownames_to_column("ensembl_with_version") %>%
  mutate(
    ensembl_id = sub("\\.\\d+$", "", ensembl_with_version),
    symbol = AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = ensembl_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    ),
    padj = ifelse(is.na(padj), 1, padj),
    log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
    regulation = case_when(
      padj < alpha & log2FoldChange >  lfc_thr ~ "Up in Disease",
      padj < alpha & log2FoldChange < -lfc_thr ~ "Up in Control",
      TRUE ~ "Not significant"
    )
  )
#Selection of top differentially expressed genes (upregulated)

top_up <- volc_WT %>%
  filter(regulation == "Up in Disease", !is.na(symbol)) %>%
  arrange(padj, desc(log2FoldChange)) %>%
  slice_head(n = 10)

top_down <- volc_WT %>%
  filter(regulation == "Up in Control", !is.na(symbol)) %>%
  arrange(padj, log2FoldChange) %>%
  slice_head(n = 10)

top_genes <- bind_rows(top_up, top_down)

#Vulcano plot with annotated genes
volcano_colors <- c(
  "Up in Disease"     = "hotpink4",
  "Up in Control"     = "steelblue4",
  "Not significant"  = "grey70"
)

ggplot(volc_WT,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = regulation)) +
  
  geom_point(alpha = 0.6, size = 1.4) +
  
  geom_vline(xintercept = c(-lfc_thr, lfc_thr),
             linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(alpha),
             linetype = "dashed", color = "darkgrey") +
  
  geom_text_repel(
    data = top_genes,
    aes(label = symbol),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.color = "grey50"
  ) +
  
  scale_color_manual(values = volcano_colors) +
  
  labs(
    title = "Volcano plot: Disease.WT vs Control.WT",
    subtitle = "Top 10 upregulated (Disease) and downregulated (Control) genes",
    x = "log2 fold change (Disease / Control)",
    y = "-log10 adjusted p-value"
  ) +
  
  theme_minimal()
#Volcano plot showing differential gene expression between Disease.WT and Control.WT. Positive log2 fold changes indicate higher expression in Disease, while negative values indicate higher expression in Control.

#Vulcano plot to look at differences between Infected.KO and Infected.WT

alpha   <- 0.01
lfc_thr <- 2

volc_KO <- as.data.frame(res_disease_KO_vs_WT) %>%
  rownames_to_column("ensembl_with_version") %>%
  mutate(
    ensembl_id = sub("\\.\\d+$", "", ensembl_with_version),
    symbol = AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = ensembl_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    ),
    padj = ifelse(is.na(padj), 1, padj),
    log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
    regulation = case_when(
      padj < alpha & log2FoldChange >  lfc_thr ~ "Higher in KO",
      padj < alpha & log2FoldChange < -lfc_thr ~ "Higher in WT",
      TRUE                                     ~ "Not significant"
    ),
    neglog10_padj_raw  = -log10(padj),
    #neglog10_padj_cap  = pmin(neglog10_padj_raw, 50)
  )

ggplot(volc_KO,
       aes(x = log2FoldChange,
           y = neglog10_padj_raw,
           color = regulation)) +
  geom_point(alpha = 0.6, size = 1.4) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr),
             linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(alpha),
             linetype = "dashed", color = "darkgrey") +
  geom_text_repel(
    data = volc_KO %>% semi_join(top_genes_KO, by = "ensembl_with_version"),
    aes(label = symbol, y = neglog10_padj_raw),
    size = 3,
    direction = "y",
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.color = "grey60",
    segment.size = 0.3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  scale_color_manual(values = volcano_colors_KO) +
  labs(
    title = "Volcano plot: Disease KO vs Disease WT",
    subtitle = "Top 10 genes higher in KO and WT (infected), padj= 0.01 and lfc_thr= 2 ",
    x = "log2 fold change (KO / WT)",
    y = "-log10 adjusted p-value"
  ) +
  theme_minimal()


# my genes of interest in the vulcano plot
goi     <- c("Cxcl10", "Stat1", "Gbp2")

volc_KO <- as.data.frame(res_disease_KO_vs_WT) %>%
  rownames_to_column("ensembl_with_version") %>%
  mutate(
    ensembl_id = sub("\\.\\d+$", "", ensembl_with_version),
    symbol = AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = ensembl_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    ),
    padj = ifelse(is.na(padj), 1, padj),
    log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
    regulation = case_when(
      padj < alpha & log2FoldChange >  lfc_thr ~ "Higher in KO",
      padj < alpha & log2FoldChange < -lfc_thr ~ "Higher in WT",
      TRUE                                     ~ "Not significant"
    ),
    neglog10_padj_raw = -log10(padj),
    GOI = symbol %in% goi
  )

volcano_colors_KO <- c(
  "Higher in KO"   = "turquoise4",
  "Higher in WT"   = "hotpink4",
  "Not significant" = "grey70"
)

ggplot(volc_KO,
       aes(x = log2FoldChange,
           y = neglog10_padj_raw)) +
  
  # all points colored by regulation
  geom_point(aes(color = regulation), alpha = 0.5, size = 1.4) +
  
  # genes of interest overplotted in red
  geom_point(data = subset(volc_KO, GOI),
             color = "red", size = 3) +
  
  #labels for top genes (your previous object) plus GOI
  geom_text_repel(
    data = subset(volc_KO, ensembl_with_version %in% top_genes_KO$ensembl_with_version | GOI),
    aes(label = symbol),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.color = "grey60",
    segment.size = 0.3
  ) +
  
  geom_vline(xintercept = c(-lfc_thr, lfc_thr),
             linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(alpha),
             linetype = "dashed", color = "darkgrey") +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  scale_color_manual(values = volcano_colors_KO, drop = TRUE) +
  
  labs(
    title = paste0("Volcano plot: Disease KO vs Disease WT (alpha = ", alpha, ")"),
    subtitle = "Top 10 genes higher in KO and WT (infected); GOI highlighted in red",
    x = "log2 fold change (KO / WT)",
    y = "-log10 adjusted p-value"
  ) +
  theme_minimal()

