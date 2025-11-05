# =============================================================================
# Pompe Disease Microarray Transcriptomics Analysis Pipeline
# Author: Fatma Tosun - Harun Bayrak
# Date: 202-07-24 (Revised: 2025-10-28)
# Description: Full pipeline for normalization, DEG analysis, annotation,
#              functional enrichment (GO/KEGG), and multiMiR analysis.
# License: MIT
# =============================================================================


# #############################################################################
# SECTION 0: SETUP AND LIBRARIES
# #############################################################################
# Install necessary CRAN packages
cran_packages <- c("pheatmap", "ggplot2", "ggrepel", "gridExtra", "grid", "umap", "dplyr")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install necessary Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c("affy", "limma", "Biobase", "hgu133plus2.db", "biomaRt", 
                   "clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE", 
                   "GOSemSim", "WebGestaltR", "multiMiR", "AnnotationDbi")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

# Load all libraries
suppressPackageStartupMessages({
  library(affy)
  library(limma)
  library(Biobase)
  library(hgu133plus2.db)
  library(AnnotationDbi)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  library(gridExtra)
  library(grid)
  library(umap)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
  library(GOSemSim)
  library(WebGestaltR)
  library(multiMiR)
})


# #############################################################################
# SECTION 1: DATA LOADING AND NORMALIZATION
# #############################################################################

# 1. Set Working Directory
setwd("C:/Users/harun/Desktop/HARUN POMPE")
cat("Working directory set to: C:/Users/harun/Desktop/HARUN POMPE\n")

# 2. Read .CEL files
gset <- ReadAffy()

# 3. Normalize (RMA)
gset_rma <- rma(gset)

# 4. Extract normalized expression matrix
expr_matrix <- exprs(gset_rma)

# 5. Impute NA/Inf values
if (any(is.na(expr_matrix))) {
  cat("NA values found, imputing with row means...\n")
  row_means <- rowMeans(expr_matrix, na.rm = TRUE)
  na_rows <- which(is.na(expr_matrix), arr.ind = TRUE)[, 1]
  expr_matrix[is.na(expr_matrix)] <- row_means[na_rows]
}
if (any(is.infinite(expr_matrix))) {
  cat("Infinite values found, capping at 99th percentile...\n")
  q99 <- quantile(expr_matrix, 0.99, na.rm = TRUE)
  expr_matrix[is.infinite(expr_matrix)] <- q99
}

cat("Normalization complete. Matrix dimensions:", dim(expr_matrix), "\n")


# #############################################################################
# SECTION 2: GROUP DEFINITION AND DEG ANALYSIS (LIMMA)
# #############################################################################

# 1. Get and clean sample names
sample_names <- colnames(expr_matrix)
sample_names <- gsub("\\.CEL$|\\.cel$", "", sample_names, ignore.case = TRUE)
colnames(expr_matrix) <- sample_names

# 2. Define group information (CORRECTED to match file names)
control_samples <- c("GSM947461", "GSM947462", "GSM947463", "GSM947464", "GSM947465", 
                   "GSM947466", "GSM947467", "GSM947468", "GSM947469", "GSM947470")
pompe_samples <- c("GSM947471", "GSM947472", "GSM947473", "GSM947474", "GSM947475", 
                 "GSM947476", "GSM947477", "GSM947478", "GSM947479")

groups <- ifelse(sample_names %in% control_samples, "Control",
                 ifelse(sample_names %in% pompe_samples, "Pompe", NA))

# 3. Check for mismatches
if (any(is.na(groups))) {
  stop("Unmatched samples found: ", paste(sample_names[is.na(groups)], collapse = ", "))
}
groups <- factor(groups, levels = c("Control", "Pompe"))
cat("Group distribution successfully created:\n")
print(table(groups))

# 4. Design and contrast matrix for limma
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
contrast_matrix <- makeContrasts(Pompe - Control, levels = design)

# 5. Run limma model
fit <- lmFit(expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 6. Get results tables
deg_results_all <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P")
deg_results_sig <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH", p.value = 0.05, lfc = 1)

# 7. Statistical summary
cat("Significant gene summary (LFC > 1 & P.adj < 0.05):\n")
de_summary <- summary(decideTests(fit2, adjust.method = "BH", p.value = 0.05, lfc = 1))
print(de_summary)


# #############################################################################
# SECTION 3: ANNOTATION AND RESULTS
# #############################################################################
cat("Annotating DEG results...\n")

# 1. Get probe IDs
probe_ids <- rownames(deg_results_all)

# 2. Annotate using AnnotationDbi::select
anno_df <- AnnotationDbi::select(hgu133plus2.db,
                                 keys = probe_ids,
                                 columns = c("SYMBOL", "GENENAME"),
                                 keytype = "PROBEID")

# 3. Merge annotations with all results
deg_results_all$PROBEID <- rownames(deg_results_all)
deg_results_all$SYMBOL <- anno_df$SYMBOL[match(deg_results_all$PROBEID, anno_df$PROBEID)]
deg_results_all$GENENAME <- anno_df$GENENAME[match(deg_results_all$PROBEID, anno_df$PROBEID)]
deg_results_all <- deg_results_all[, c("PROBEID", "SYMBOL", "GENENAME", "logFC", "P.Value", "adj.P.Val", "B", "t", "AveExpr")]

# 4. Update significant results table
deg_results_sig$PROBEID <- rownames(deg_results_sig)
deg_results_sig <- merge(deg_results_sig, 
                         deg_results_all[, c("PROBEID", "SYMBOL", "GENENAME")], 
                         by = "PROBEID")
deg_results_sig <- deg_results_sig[, c("PROBEID", "SYMBOL", "GENENAME", "logFC", "P.Value", "adj.P.Val", "B", "t", "AveExpr")]

# 5. Save results to CSV
write.csv(deg_results_all, "Pompe_vs_Control_All_Genes.csv", row.names = FALSE)
write.csv(deg_results_sig, "Pompe_vs_Control_Significant_Genes.csv", row.names = FALSE)

cat("Annotation complete and results saved to CSV.\n")


# #############################################################################
# SECTION 4: QC AND VISUALIZATION (UMAP, VOLCANO, HEATMAP)
# #############################################################################

# -----------------------------------
# 4A. UMAP (Quality Control)
# -----------------------------------
cat("Running UMAP analysis...\n")
expr_matrix_clean <- na.omit(expr_matrix)
expr_matrix_clean <- expr_matrix_clean[!duplicated(rownames(expr_matrix_clean)), ]

ump <- umap(t(expr_matrix_clean), n_neighbors = 8, random_state = 123)
ump_df <- data.frame(
  Sample = rownames(ump$layout),
  UMAP1 = ump$layout[, 1],
  UMAP2 = ump$layout[, 2],
  Group = groups
)
ump_df$Sample_clean <- sapply(strsplit(ump_df$Sample, "_"), `[`, 1)

umap_plot <- ggplot(ump_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = Sample_clean), size = 3.5, max.overlaps = 15, color = "black") +
  labs(title = "UMAP Plot of Samples (n_neighbors = 8)",
       x = "UMAP1", y = "UMAP2") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "right")

ggsave("UMAP_Plot.pdf", plot = umap_plot, width = 8, height = 6)
print(umap_plot)

# -----------------------------------
# 4B. VOLCANO PLOT
# -----------------------------------
cat("Drawing Volcano plot...\n")
deg_results_all$category <- "NS"
deg_results_all$category[deg_results_all$adj.P.Val < 0.05 & abs(deg_results_all$logFC) > 1] <- "p-value and log2 FC"
deg_results_all$category[deg_results_all$adj.P.Val < 0.05 & abs(deg_results_all$logFC) <= 1] <- "p-value"
deg_results_all$category[deg_results_all$adj.P.Val >= 0.05 & abs(deg_results_all$logFC) > 1] <- "log2 FC"

colors <- c("NS" = "gray", "log2 FC" = "skyblue", "p-value" = "red", "p-value and log2 FC" = "blue")

top_genes_labels <- deg_results_all %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::arrange(adj.P.Val) %>%
  head(20)

volcano_plot <- ggplot(deg_results_all, aes(x = logFC, y = -log10(adj.P.Val), color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = top_genes_labels,
    aes(label = SYMBOL),
    color = "black", size = 3.5, box.padding = 0.8, max.overlaps = Inf
  ) +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  labs(
    title = "Volcano Plot: Pompe vs Control",
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~Adjusted~P-Value),
    color = "Category"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

print(volcano_plot)

# -----------------------------------
# 4C. HEATMAP
# -----------------------------------
cat("Drawing Heatmap...\n")
top50_genes_probes <- head(deg_results_sig$PROBEID, 50)
top50_expr <- expr_matrix[top50_genes_probes, ]

top50_anno <- deg_results_all[match(top50_genes_probes, deg_results_all$PROBEID), c("PROBEID", "SYMBOL")]
top50_anno$Label <- ifelse(is.na(top50_anno$SYMBOL) | top50_anno$SYMBOL == "", top50_anno$PROBEID, top50_anno$SYMBOL)
top50_anno$Label <- make.unique(top50_anno$Label)
rownames(top50_expr) <- top50_anno$Label

annotation_col <- data.frame(Group = groups)
rownames(annotation_col) <- colnames(top50_expr)

heatmap_plot_obj <- pheatmap(
  top50_expr,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  main = "Heatmap of Top 50 DEGs",
  fontsize_row = 8,
  silent = TRUE
)

print(heatmap_plot_obj)

# -----------------------------------
# 4D. Combined Plot (Volcano + Heatmap)
# -----------------------------------
cat("Saving combined figure as PDF...\n")
pdf("Figure1_Volcano_Heatmap.pdf", width = 16, height = 8)
grid.arrange(
  arrangeGrob(volcano_plot, top = textGrob("A", gp = gpar(fontsize=20, fontface="bold"), x = unit(0, "npc"), hjust = -0.5)),
  arrangeGrob(heatmap_plot_obj$gtable, top = textGrob("B", gp = gpar(fontsize=20, fontface="bold"), x = unit(0, "npc"), hjust = -0.5)),
  ncol = 2
)
dev.off()


# #############################################################################
# SECTION 5: FUNCTIONAL ENRICHMENT (GO / KEGG)
# #############################################################################
cat("Starting GO and KEGG analysis...\n")

# 1. Prepare significant gene list
significant_genes <- deg_results_sig$SYMBOL
significant_genes <- unique(na.omit(significant_genes))
cat(length(significant_genes), "unique significant gene symbols found for enrichment.\n")

# 2. GO Biological Process (BP) Analysis
go_bp <- enrichGO(
  gene = significant_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
go_bp_df <- as.data.frame(go_bp)
write.csv(go_bp_df, "GO_Biological_Process.csv", row.names = FALSE)

# 3. KEGG Pathway Analysis
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = significant_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
entrez_ids <- na.omit(entrez_ids)

kegg <- enrichKEGG(
  gene = entrez_ids,
  organism = 'hsa',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
kegg_df <- as.data.frame(kegg)
write.csv(kegg_df, "KEGG_Enrichment.csv", row.names = FALSE)

# 4. GO/KEGG Visualization
if (!is.null(go_bp) && nrow(go_bp_df) > 0) {
  p_go_dot <- dotplot(go_bp, showCategory = 15, title = "GO Biological Process")
  ggsave("GO_BP_dotplot.pdf", plot = p_go_dot, width = 10, height = 8)
  print(p_go_dot)
} else {
  cat("No significant GO results found.\n")
}

if (!is.null(kegg) && nrow(kegg_df) > 0) {
  p_kegg_dot <- dotplot(kegg, showCategory = 15, title = "KEGG Pathways")
  ggsave("KEGG_dotplot.pdf", plot = p_kegg_dot, width = 10, height = 8)
  print(p_kegg_dot)
} else {
  cat("No significant KEGG results found.\n")
}


# #############################################################################
# SECTION 6: miRNA TARGET ANALYSIS (multiMiR)
# #############################################################################
cat("Starting multiMiR analysis for", length(significant_genes), "DEGs...\n")
cat("WARNING: This section involves heavy queries and may take several hours.\n")

# -----------------------------------
# 6.1. Initial Query (Predicted & Validated Targets for all DEGs)
# -----------------------------------
# Query 1: Predicted
mm_predicted <- tryCatch({
  get_multimir(
    org = "hsa",
    target = significant_genes,
    table = "predicted",
    summary = TRUE
  )
}, error = function(e) {
  message("multiMiR (predicted) query error: ", e$message)
  return(NULL)
})

# Query 2: Validated
mm_validated <- tryCatch({
  get_multimir(
    org = "hsa",
    target = significant_genes,
    table = "validated",
    summary = TRUE
  )
}, error = function(e) {
  message("multiMiR (validated) query error: ", e$message)
  return(NULL)
})

# -----------------------------------
# 6.2. Combine Results and Summarize
# -----------------------------------
if (is.null(mm_predicted) && is.null(mm_validated)) {
  stop("FATAL ERROR: No results returned from multiMiR. Check internet connection.")
} else {
  
  data_list <- list()
  if (!is.null(mm_predicted)) { data_list$pred <- mm_predicted@data }
  if (!is.null(mm_validated)) { data_list$val <- mm_validated@data }
  
  # Combine all results
  mm_all_data <- dplyr::bind_rows(data_list)
  
  if (nrow(mm_all_data) == 0) {
    cat("Query successful, but no miRNA interactions found for these DEGs.\n")
  } else {
    cat("multiMiR query complete.", nrow(mm_all_data), "total interactions found.\n")
    cat("Summarizing miRNAs by unique DEG targets...\n")
    
    # Create the main summary table (mirna_target_counts)
    # This table is sorted by 'unique_DEG_target_count' by default
    mirna_target_counts <- mm_all_data %>%
      dplyr::filter(target_symbol %in% significant_genes) %>%
      dplyr::group_by(mature_mirna_id) %>%
      dplyr::summarise(
        unique_DEG_target_count = n_distinct(target_symbol),
        total_interactions = n()
      ) %>%
      dplyr::arrange(desc(unique_DEG_target_count), desc(total_interactions))
    
    write.csv(mirna_target_counts, "miRNA_DEG_Target_Counts_Summary.csv", row.names = FALSE)
    cat("miRNA summary saved to 'miRNA_DEG_Target_Counts_Summary.csv'\n")
  }
}

# -----------------------------------
# 6.3. Create Final Plot: Top 20 by Total Interactions (Light Colors)
# -----------------------------------
cat("Generating Top 20 miRNA plot (sorted by Total Interactions)...\n")

# 1. Prepare data for stacking
# !!! MODIFICATION: hsa-let-7 filter is REMOVED as requested !!!
data_for_stacking <- mm_all_data %>%
  dplyr::filter(target_symbol %in% significant_genes)

# 2. Get database-specific counts
mirna_db_counts <- data_for_stacking %>%
  dplyr::group_by(mature_mirna_id, database) %>%
  dplyr::summarise(unique_DEG_target_count = n_distinct(target_symbol), .groups = 'drop')

# 3. Create sorting score based on TOTAL interactions
mirna_total_interactions_score <- data_for_stacking %>%
  dplyr::group_by(mature_mirna_id) %>%
  dplyr::summarise(total_interactions = n()) %>%
  dplyr::arrange(desc(total_interactions))

# 4. Get names of the Top 20 miRNAs by this new score
top_20_mirna_names_by_total <- head(mirna_total_interactions_score$mature_mirna_id, 20)

# 5. Filter the database-specific data for just these Top 20
top_20_stacked_data_by_total <- mirna_db_counts %>%
  dplyr::filter(mature_mirna_id %in% top_20_mirna_names_by_total)

# 6. Set factor levels to ensure correct plot order
top_20_stacked_data_by_total$mature_mirna_id <- factor(
  top_20_stacked_data_by_total$mature_mirna_id,
  levels = rev(top_20_mirna_names_by_total) 
)

# 7. Define the light color palette
new_light_palette <- c(
  "diana_microt" = "#E57373", "elmmo" = "#FFB74D", "microcosm" = "#AED581",
  "miranda" = "#81C784", "mirdb" = "#4DB6AC", "mirecords" = "#4DD0E1",
  "mirtarbase" = "#64B5F6", "pictar" = "#9575CD", "pita" = "#7E57C2",
  "tarbase" = "#BA68C8", "targetscan" = "#5C6BC0"
)

# 8. Create the final plot (CORRECTED subtitle)
p_mirna_stacked_final_subtitle <- ggplot(top_20_stacked_data_by_total, 
                                         aes(x = mature_mirna_id, y = unique_DEG_target_count, fill = database)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  labs(
    title = "Top 20 miRNA",
    subtitle = "by unique DEG targets", # Corrected subtitle
    x = "miRNA",
    y = "Number of Unique DEG Targets",
    fill = "Database"
  ) +
  scale_fill_manual(values = new_light_palette, drop = FALSE) +
  theme_minimal(base_size = 12) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) 

print(p_mirna_stacked_final_subtitle)
ggsave("miRNA_Top20_Final_Plot.pdf", plot = p_mirna_stacked_final_subtitle, width = 12, height = 8)
cat("Top 20 stacked miRNA plot (Final Version) saved as 'miRNA_Top20_Final_Plot.pdf'\n")


# #############################################################################
# SECTION 7: SPECIFIC GENE QUERY (GAA)
# #############################################################################

# -----------------------------------
# 7.1. GAA Target Analysis (Corrected Logic: New Query)
# -----------------------------------
cat("Identifying Top 100 miRNAs from *original* DEG list for GAA query...\n")

# 1. Get Top 100 list from the ORIGINAL summary (mirna_target_counts)
# This list is sorted by 'unique_DEG_target_count'
top_100_original_list <- head(mirna_target_counts$mature_mirna_id, 100)

cat("Starting NEW query for 'GAA' targets of these Top 100 miRNAs...\n")
cat("(This may take several minutes...)\n")

# 2. Run NEW query using 'mirna' parameter for the Top 100 list
top_100_targets <- tryCatch({
  get_multimir(
    mirna = top_100_original_list,
    table = "predicted",
    summary = TRUE
  )
}, error = function(e) {
  message("multiMiR (GAA query) error: ", e$message)
  return(NULL)
})

# 3. Process the results
if (is.null(top_100_targets) || nrow(top_100_targets@data) == 0) {
  cat("Error: Could not retrieve targets for the Top 100 miRNAs.\n")
} else {
  cat("GAA query complete.", nrow(top_100_targets@data), "total targets found.\n")
  
  # 4. Filter for TargetScan and 'GAA'
  targetscan_gaa_binders <- top_100_targets@data %>%
    dplyr::filter(database == "targetscan") %>%
    dplyr::filter(target_symbol == "GAA")
    
  # 5. Get unique results
  unique_targetscan_binders <- targetscan_gaa_binders %>%
    dplyr::select(mature_mirna_id, target_symbol, database) %>%
    dplyr::distinct() %>%
    dplyr::arrange(mature_mirna_id)
    
  # 6. Report findings
  count_of_binders <- nrow(unique_targetscan_binders)
  
  cat("-----------------------------------------------------------------\n")
  cat("RESULT (TargetScan):\n")
  cat(count_of_binders, "out of the Top 100 DEG-targeting miRNAs also target 'GAA' (according to TargetScan).\n")
  cat("-----------------------------------------------------------------\n\n")
  
  if (count_of_binders > 0) {
    cat("miRNAs targeting 'GAA' (from TargetScan):\n")
    
    # CORRECTED: Use 'as.data.frame' to prevent print errors
    print(as.data.frame(unique_targetscan_binders))
    
    write.csv(unique_targetscan_binders, "Top100_miRNA_GAA_TargetScan_Hits.csv", row.names = FALSE)
    cat("\nList saved to 'Top100_miRNA_GAA_TargetScan_Hits.csv'\n")
  } else {
    cat("No miRNAs from the Top 100 list were found to target 'GAA' in TargetScan.\n")
  }
}

cat("\n--- Full Analysis Pipeline Finished. ---\n")
