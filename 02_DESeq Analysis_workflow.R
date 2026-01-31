# Differential expression analysis

#Load packages
library(DESeq2) 
library(ggplot2)
library(pheatmap)
library(tidyr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)

#1)Infection effect in WT: Disease.WT vs Control.WT, extraction of DE results for pairwise contrast
# dds_interaction already has design ~ group, where group = Condition:Mouse.genotype
resultsNames(dds_interaction)

res_infection_WT_vs_control <- results(
  dds_interaction,
  contrast = c("group", "Disease.WT", "Control.WT")
)
summary(res_infection_WT_vs_control)

# Count DE genes at padj < 0.01
alpha <- 0.01
res_clean_WT <- res_infection_WT_vs_control[!is.na(res_infection_WT_vs_control$padj), ]

de_WT     <- res_clean_WT$padj < alpha
n_de_WT   <- sum(de_WT)
n_up_WT   <- sum(de_WT & res_clean_WT$log2FoldChange > 2)  # higher in Disease.WT
n_down_WT <- sum(de_WT & res_clean_WT$log2FoldChange < -2)  # higher in Control.WT

cat("WT infection effect (Disease.WT vs Control.WT), padj <", alpha, ":\n")
cat("Total DE genes:", n_de_WT, "\n")
cat("Up in Disease.WT:", n_up_WT, "\n")
cat("Down in Control.WT:", n_down_WT, "\n\n")

summary(res_disease_KO_vs_WT)

#2) Infected KO vs infected WT: Disease.Ifnar-/-x Ifngr-/- vs Disease.WT
# group factor was already defined earlier:
# SampleTable$group <- interaction(SampleTable$Condition, SampleTable$Mouse.genotype)

levels(colData(dds_interaction)$group) #checking if group contains the combination needed in the next step

res_disease_KO_vs_WT <- results(
  dds_interaction,
  contrast = c("group", "Disease.Ifnar-/-x Ifngr-/-", "Disease.WT")
)
summary(res_disease_KO_vs_WT)

# Count DE genes at padj < 0.01
res_clean_KO <- res_disease_KO_vs_WT[!is.na(res_disease_KO_vs_WT$padj), ]

de_KO   <- res_clean_KO$padj < alpha
n_de_KO <- sum(de_KO)
n_up_KO <- sum(de_KO & res_clean_KO$log2FoldChange > 2)  # higher in Disease.Ifnar-/-x Ifngr-/-
n_down_KO <- sum(de_KO & res_clean_KO$log2FoldChange < -2)  # higher in Disease.WT

cat("Infected KO vs infected WT (Disease.Ifnar-/-x Ifngr-/- vs Disease.WT), padj <", alpha, ":\n")
cat("Total DE genes:", n_de_KO, "\n")
cat("Up in Disease.Ifnar-/-x Ifngr-/-:", n_up_KO, "\n")
cat("Up in Disease.WT:", n_down_KO, "\n")


#Summary table
df_summary <- data.frame(
  contrast   = c("Disease.WT vs Control.WT",
                 "Disease.Ifnar-/-x Ifngr-/- vs Disease.WT"),
  padj_cutoff = c(0.01, 0.01),
  DE_genes   = c(n_de_WT,   n_de_KO),
  up_genes   = c(n_up_WT,   n_up_KO),
  down_genes = c(n_down_WT, n_down_KO)
)
df_summary

knitr::kable(df_summary,
             caption = "Summary of differentially expressed genes for each contrast.")

#numbers changed because of threshold

#Additional visualization using volcano plot for differentially expressed genes

#1)Volcano plot: DE genes in infected WT vs control WT

alpha  <- 0.01   # padj cutoff
lfc_thr <- 2     # log2FC threshold (≈ 4‑fold change)

volc_WT <- as.data.frame(res_infection_WT_vs_control)
volc_WT$gene <- rownames(volc_WT)

#NAs replaced
volc_WT$padj[is.na(volc_WT$padj)] <- 1
volc_WT$log2FoldChange[is.na(volc_WT$log2FoldChange)] <- 0

#Regulation categories defined
volc_WT$regulation <- "Not significant"
volc_WT$regulation[volc_WT$padj < alpha & volc_WT$log2FoldChange >  lfc_thr] <- "Up in Disease"
volc_WT$regulation[volc_WT$padj < alpha & volc_WT$log2FoldChange < -lfc_thr] <- "Up in Control"

#Volcano plot specification
ggplot(volc_WT, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "darkgrey") +
  labs(title = "Volcano plot: Disease.WT vs Control.WT",
       x = "log2 fold change (Disease / Control)",
       y = "-log10 adjusted p-value") +
  theme_minimal()

#2)Volcano plot: DE genes in infected Mutant vs infected WT
alpha  <- 0.01 #adjusted p-value
lfc_thr <- 2

volc_infection <- as.data.frame(res_disease_KO_vs_WT)
volc_infection$gene <- rownames(volc_infection)

#NAs replaced
volc_infection$padj[is.na(volc_infection$padj)] <- 1
volc_infection$log2FoldChange[is.na(volc_infection$log2FoldChange)] <- 0

#Regulation categories defined
volc_infection$regulation <- "Not significant"
volc_infection$regulation[volc_infection$padj < alpha & volc_infection$log2FoldChange >  lfc_thr] <- "Up in KO"
volc_infection$regulation[volc_infection$padj < alpha & volc_infection$log2FoldChange < -lfc_thr] <- "Up in WT"

#Volcano plot specification
ggplot(volc_infection, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "darkgrey") +
  labs(title = "Volcano plot: Disease.Mutant vs Disease.WT",
       x = "log2 fold change (KO / WT)",
       y = "-log10 adjusted p-value") +
  theme_minimal()

#Vulcano plot with annotated genes

#Genes of interest: Cxcl10, Stat1 and Gbp2
#gene annotation install: BiocManager::install("org.Mm.eg.db") and BiocManager::install("AnnotationDbi")
#1)Getting Ensembl IDs without version as a new column
ensembl_with_version <- rownames(dds_interaction)
ensembl_id <- sub("\\.\\d+$", "", ensembl_with_version) #remove .1, .2 etc

#2) Map Ensembl IDs to gene symbols
symbol_map <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys     = ensembl_id,
  column   = "SYMBOL",
  keytype  = "ENSEMBL",
  multiVals = "first"
)

#3)Saving annotation on the dds object
rowData(dds_interaction)$ensembl_id <- ensembl_id
rowData(dds_interaction)$symbol     <- symbol_map

#4)Adding symbols to normalized counts:
norm_counts <- counts(dds_interaction, normalized = TRUE)
rownames(norm_counts) <- symbol_map   # rownames are gene symbols, with NAs for unmapped

#5)Dropping rows with NA symbol, to have only genes with symbols
keep <- !is.na(rownames(norm_counts))
norm_counts <- norm_counts[keep, ]
dds_interaction <- dds_interaction[keep, ]

#6) Checking if genes can be found
genes_of_interest <- c("Cxcl10", "Stat1", "Gbp2")
genes_of_interest %in% rownames(norm_counts)

#Gene expression plots for Cxcl10, Stat1 and Gbp2

keep <- rownames(norm_counts) %in% genes_of_interest
df <- as.data.frame(t(norm_counts[keep, ]))
df$group <- colData(dds_interaction)$group
df$sample <- rownames(df)

df_long <- pivot_longer(df,
                        cols = all_of(genes_of_interest),
                        names_to = "gene",
                        values_to = "count")

ggplot(df_long, aes(x = group, y = count, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.4) +
  facet_wrap(~ gene, scales = "free_x") +
  scale_y_log10() +
  labs(y = "Normalized read counts (log10)", x = "group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Analysis of statistical significance in changes of gene expression levels between conditions
#Wilcoxon on log counts, test within each genotype: Disease vs Control.
stats_df <- df_long %>%
  filter(gene %in% genes_of_interest,
         group %in% c("Control.WT","Disease.WT",
                      "Control.Mutant","Disease.Mutant")) %>%
  mutate(genotype = ifelse(grepl("WT", group), "WT", "Mutant"),
         condition = ifelse(grepl("Control", group), "Control", "Disease")) %>%
  group_by(gene, genotype) %>%
  summarize(
    p = wilcox.test(log10(count) ~ condition)$p.value,
    .groups = "drop"
  ) %>%
  mutate(star = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  ))

label_pos <- df_long %>%
  filter(gene %in% genes_of_interest) %>%
  group_by(gene) %>%
  summarize(y = max(count)*1.5, .groups = "drop")

plot_df <- df_long %>%
  filter(gene %in% genes_of_interest)

ggplot(plot_df, aes(x = group, y = count, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.4) +
  facet_wrap(~ gene, scales = "free_x") +
  scale_y_log10() +
  labs(y = "Normalized read counts (log10)", x = "group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(
    data = stats_df %>%
      left_join(label_pos, by = "gene") %>%
      mutate(group = ifelse(genotype == "WT",
                            "Disease.WT",
                            "Disease.Mutant")),
    aes(x = group, y = y, label = star),
    inherit.aes = FALSE,
    vjust = 0
  )
