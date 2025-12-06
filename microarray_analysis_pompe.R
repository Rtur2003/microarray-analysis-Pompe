# =============================================================================
# Pompe Disease Microarray Transcriptomics Analysis Pipeline
# Authors: Fatma Tosun - Harun Bayrak
# Date: 2025-12-06
# Description: Normalisation, DEG analysis, annotation, enrichment, QC plots,
#              and optional multiMiR queries for GSE38680.
# License: MIT
# =============================================================================


# #############################################################################
# SECTION 0: SETUP AND CONFIG
# #############################################################################

required_packages <- c(
  "yaml", "affy", "limma", "Biobase", "hgu133plus2.db", "AnnotationDbi",
  "pheatmap", "ggplot2", "ggrepel", "gridExtra", "grid", "umap", "dplyr",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE", "GOSemSim", "multiMiR"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Install them via the Conda environment (environment/environment.yml) before running."
  )
}

suppressPackageStartupMessages({
  library(yaml)
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
  library(multiMiR)
})

default_params <- list(
  deg = list(
    design = "~ 0 + group",
    contrast = "Pompe - Control",
    logFC = 1.0,
    padj = 0.05
  ),
  plots = list(
    volcano_topn = 20
  ),
  enrichment = list(
    ontologies = c("BP", "CC", "MF")
  ),
  multimir = list(
    enabled = TRUE
  )
)

config_path <- file.path("config", "params.yaml")
params <- tryCatch(
  modifyList(default_params, yaml::read_yaml(config_path)),
  error = function(e) {
    message("Could not read config/params.yaml, using defaults: ", e$message)
    default_params
  }
)

results_dir <- Sys.getenv("RESULTS_DIR", unset = "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
}
cat("Results will be written to:", normalizePath(results_dir, winslash = "/"), "\n")

set.seed(123)


# #############################################################################
# SECTION 1: DATA LOADING AND NORMALIZATION
# #############################################################################

samplesheet_path <- file.path("config", "samplesheet.csv")
if (!file.exists(samplesheet_path)) {
  stop("Sample sheet not found at ", samplesheet_path)
}

samplesheet <- read.csv(samplesheet_path, stringsAsFactors = FALSE)
required_cols <- c("sample_id", "group", "filename")
missing_cols <- setdiff(required_cols, names(samplesheet))
if (length(missing_cols) > 0) {
  stop("Missing required columns in sample sheet: ", paste(missing_cols, collapse = ", "))
}

cel_files <- file.path("metadata", samplesheet$filename)
if (!length(cel_files)) {
  stop("No CEL files listed in the sample sheet.")
}
missing_cel <- cel_files[!file.exists(cel_files)]
if (length(missing_cel) > 0) {
  stop("Missing CEL files under metadata/: ", paste(basename(missing_cel), collapse = ", "))
}

cat("Reading CEL files...\n")
gset <- ReadAffy(filenames = cel_files)

cat("Running RMA normalization...\n")
gset_rma <- rma(gset)
expr_matrix <- exprs(gset_rma)
colnames(expr_matrix) <- samplesheet$sample_id

# Impute NA/Inf values defensively
if (anyNA(expr_matrix)) {
  cat("NA values found; imputing with row means.\n")
  row_means <- rowMeans(expr_matrix, na.rm = TRUE)
  idx_na <- which(is.na(expr_matrix), arr.ind = TRUE)
  expr_matrix[idx_na] <- row_means[idx_na[, 1]]
}
if (any(!is.finite(expr_matrix))) {
  cat("Infinite values found; capping at 99th percentile.\n")
  q99 <- quantile(expr_matrix[is.finite(expr_matrix)], 0.99)
  expr_matrix[!is.finite(expr_matrix)] <- q99
}

cat("Normalization complete. Matrix dimensions:", paste(dim(expr_matrix), collapse = " x "), "\n")


# #############################################################################
# SECTION 2: GROUP DEFINITION AND DEG ANALYSIS (LIMMA)
# #############################################################################

groups <- factor(samplesheet$group)
if (any(is.na(groups))) {
  stop("Group labels contain NA. Please fix config/samplesheet.csv.")
}
cat("Group distribution:\n")
print(table(groups))

design <- model.matrix(as.formula(params$deg$design), data = data.frame(group = groups))
colnames(design) <- sub("^group", "", colnames(design))
contrast_matrix <- makeContrasts(contrasts = params$deg$contrast, levels = design)

fit <- lmFit(expr_matrix, design)
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

deg_results_all <- topTable(
  fit2,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)
deg_results_sig <- topTable(
  fit2,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  p.value = params$deg$padj,
  lfc = params$deg$logFC
)

cat("Significant genes (BH <", params$deg$padj, " & |log2FC| >=", params$deg$logFC, "):\n")
print(summary(decideTests(fit2, adjust.method = "BH", p.value = params$deg$padj, lfc = params$deg$logFC)))


# #############################################################################
# SECTION 3: ANNOTATION AND RESULTS
# #############################################################################
cat("Annotating DEG results...\n")

probe_ids <- rownames(deg_results_all)
anno_df <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

deg_results_all$PROBEID <- rownames(deg_results_all)
deg_results_all$SYMBOL <- anno_df$SYMBOL[match(deg_results_all$PROBEID, anno_df$PROBEID)]
deg_results_all$GENENAME <- anno_df$GENENAME[match(deg_results_all$PROBEID, anno_df$PROBEID)]
deg_results_all <- deg_results_all[, c("PROBEID", "SYMBOL", "GENENAME", "logFC", "P.Value", "adj.P.Val", "B", "t", "AveExpr")]

deg_results_sig$PROBEID <- rownames(deg_results_sig)
deg_results_sig <- merge(
  deg_results_sig,
  deg_results_all[, c("PROBEID", "SYMBOL", "GENENAME")],
  by = "PROBEID",
  all.x = TRUE
)
deg_results_sig <- deg_results_sig[, c("PROBEID", "SYMBOL", "GENENAME", "logFC", "P.Value", "adj.P.Val", "B", "t", "AveExpr")]

write.csv(deg_results_all, file.path(results_dir, "Pompe_vs_Control_All_Genes.csv"), row.names = FALSE)
write.csv(deg_results_sig, file.path(results_dir, "Pompe_vs_Control_Significant_Genes.csv"), row.names = FALSE)

cat("Annotation complete and results saved to results/.\n")


# #############################################################################
# SECTION 4: QC AND VISUALIZATION (UMAP, VOLCANO, HEATMAP)
# #############################################################################

# 4A. UMAP
cat("Running UMAP analysis...\n")
expr_matrix_clean <- expr_matrix[!duplicated(rownames(expr_matrix)), ]
expr_matrix_clean <- expr_matrix_clean[complete.cases(expr_matrix_clean), ]
ump <- umap(t(expr_matrix_clean), n_neighbors = 8, random_state = 123)
ump_df <- data.frame(
  Sample = rownames(ump$layout),
  UMAP1 = ump$layout[, 1],
  UMAP2 = ump$layout[, 2],
  Group = groups
)

umap_plot <- ggplot(ump_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 15, color = "black") +
  labs(title = "UMAP of Samples", x = "UMAP1", y = "UMAP2") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "right")

ggsave(file.path(results_dir, "UMAP_Plot.pdf"), plot = umap_plot, width = 8, height = 6)

# 4B. VOLCANO
cat("Drawing Volcano plot...\n")
deg_results_all$category <- "NS"
deg_results_all$category[deg_results_all$adj.P.Val < params$deg$padj & abs(deg_results_all$logFC) > params$deg$logFC] <- "p-value and log2 FC"
deg_results_all$category[deg_results_all$adj.P.Val < params$deg$padj & abs(deg_results_all$logFC) <= params$deg$logFC] <- "p-value"
deg_results_all$category[deg_results_all$adj.P.Val >= params$deg$padj & abs(deg_results_all$logFC) > params$deg$logFC] <- "log2 FC"

colors <- c("NS" = "gray", "log2 FC" = "skyblue", "p-value" = "red", "p-value and log2 FC" = "blue")

top_genes_labels <- deg_results_all %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::arrange(adj.P.Val) %>%
  head(params$plots$volcano_topn)

volcano_plot <- ggplot(deg_results_all, aes(x = logFC, y = -log10(adj.P.Val), color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = top_genes_labels,
    aes(label = SYMBOL),
    color = "black",
    size = 3.5,
    box.padding = 0.8,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = c(-params$deg$logFC, params$deg$logFC), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(params$deg$padj), linetype = "dashed", color = "gray50") +
  labs(
    title = "Volcano Plot: Pompe vs Control",
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~Adjusted~P-Value),
    color = "Category"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(file.path(results_dir, "Volcano_Plot.pdf"), plot = volcano_plot, width = 8, height = 6)

# 4C. HEATMAP
cat("Drawing Heatmap...\n")
if (nrow(deg_results_sig) >= 2) {
  top50_genes_probes <- head(deg_results_sig$PROBEID, 50)
  top50_expr <- expr_matrix[top50_genes_probes, , drop = FALSE]

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

  ggsave(
    filename = file.path(results_dir, "Heatmap_Top50_DEGs.pdf"),
    plot = heatmap_plot_obj$gtable,
    width = 8,
    height = 8
  )
} else {
  cat("Not enough significant genes for heatmap; skipping.\n")
}

# 4D. Combined Plot (Volcano + Heatmap)
if (exists("heatmap_plot_obj")) {
  pdf(file.path(results_dir, "Figure1_Volcano_Heatmap.pdf"), width = 16, height = 8)
  grid.arrange(
    arrangeGrob(volcano_plot, top = textGrob("A", gp = gpar(fontsize = 20, fontface = "bold"), x = unit(0, "npc"), hjust = -0.5)),
    arrangeGrob(heatmap_plot_obj$gtable, top = textGrob("B", gp = gpar(fontsize = 20, fontface = "bold"), x = unit(0, "npc"), hjust = -0.5)),
    ncol = 2
  )
  dev.off()
}


# #############################################################################
# SECTION 5: FUNCTIONAL ENRICHMENT (GO / KEGG)
# #############################################################################
cat("Starting GO and KEGG analysis...\n")

significant_genes <- unique(na.omit(deg_results_sig$SYMBOL))
cat(length(significant_genes), "unique significant gene symbols found for enrichment.\n")

if (length(significant_genes) > 0) {
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
  write.csv(go_bp_df, file.path(results_dir, "GO_Biological_Process.csv"), row.names = FALSE)

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
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  kegg_df <- as.data.frame(kegg)
  write.csv(kegg_df, file.path(results_dir, "KEGG_Enrichment.csv"), row.names = FALSE)

  if (!is.null(go_bp) && nrow(go_bp_df) > 0) {
    p_go_dot <- dotplot(go_bp, showCategory = 15, title = "GO Biological Process")
    ggsave(file.path(results_dir, "GO_BP_dotplot.pdf"), plot = p_go_dot, width = 10, height = 8)
  } else {
    cat("No significant GO results found.\n")
  }

  if (!is.null(kegg) && nrow(kegg_df) > 0) {
    p_kegg_dot <- dotplot(kegg, showCategory = 15, title = "KEGG Pathways")
    ggsave(file.path(results_dir, "KEGG_dotplot.pdf"), plot = p_kegg_dot, width = 10, height = 8)
  } else {
    cat("No significant KEGG results found.\n")
  }
} else {
  cat("No significant genes for enrichment; skipping.\n")
}


# #############################################################################
# SECTION 6: miRNA TARGET ANALYSIS (multiMiR)
# #############################################################################
skip_multimir <- !isTRUE(params$multimir$enabled) || identical(Sys.getenv("SKIP_MULTIMIR"), "1")
if (skip_multimir) {
  cat("Skipping multiMiR analysis (disabled via config or SKIP_MULTIMIR env var).\n")
} else if (length(significant_genes) == 0) {
  cat("No significant genes available; skipping multiMiR analysis.\n")
} else {
  cat("Starting multiMiR analysis for", length(significant_genes), "DEGs...\n")
  cat("Note: This step may be slow and requires internet connectivity.\n")

  mm_predicted <- tryCatch({
    get_multimir(org = "hsa", target = significant_genes, table = "predicted", summary = TRUE)
  }, error = function(e) {
    message("multiMiR (predicted) query error: ", e$message)
    NULL
  })

  mm_validated <- tryCatch({
    get_multimir(org = "hsa", target = significant_genes, table = "validated", summary = TRUE)
  }, error = function(e) {
    message("multiMiR (validated) query error: ", e$message)
    NULL
  })

  if (is.null(mm_predicted) && is.null(mm_validated)) {
    warning("No results returned from multiMiR.")
  } else {
    data_list <- list()
    if (!is.null(mm_predicted)) data_list$pred <- mm_predicted@data
    if (!is.null(mm_validated)) data_list$val <- mm_validated@data

    mm_all_data <- dplyr::bind_rows(data_list)

    if (nrow(mm_all_data) == 0) {
      cat("multiMiR query successful but no interactions found.\n")
    } else {
      cat("multiMiR query complete.", nrow(mm_all_data), "total interactions found.\n")

      mirna_target_counts <- mm_all_data %>%
        dplyr::filter(target_symbol %in% significant_genes) %>%
        dplyr::group_by(mature_mirna_id) %>%
        dplyr::summarise(
          unique_DEG_target_count = n_distinct(target_symbol),
          total_interactions = n(),
          .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(unique_DEG_target_count), dplyr::desc(total_interactions))

      write.csv(mirna_target_counts, file.path(results_dir, "miRNA_DEG_Target_Counts_Summary.csv"), row.names = FALSE)

      data_for_stacking <- mm_all_data %>%
        dplyr::filter(target_symbol %in% significant_genes)

      mirna_db_counts <- data_for_stacking %>%
        dplyr::group_by(mature_mirna_id, database) %>%
        dplyr::summarise(unique_DEG_target_count = n_distinct(target_symbol), .groups = "drop")

      mirna_total_interactions_score <- data_for_stacking %>%
        dplyr::group_by(mature_mirna_id) %>%
        dplyr::summarise(total_interactions = n(), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(total_interactions))

      top_20_mirna_names_by_total <- head(mirna_total_interactions_score$mature_mirna_id, 20)

      top_20_stacked_data_by_total <- mirna_db_counts %>%
        dplyr::filter(mature_mirna_id %in% top_20_mirna_names_by_total)

      top_20_stacked_data_by_total$mature_mirna_id <- factor(
        top_20_stacked_data_by_total$mature_mirna_id,
        levels = rev(top_20_mirna_names_by_total)
      )

      new_light_palette <- c(
        "diana_microt" = "#E57373", "elmmo" = "#FFB74D", "microcosm" = "#AED581",
        "miranda" = "#81C784", "mirdb" = "#4DB6AC", "mirecords" = "#4DD0E1",
        "mirtarbase" = "#64B5F6", "pictar" = "#9575CD", "pita" = "#7E57C2",
        "tarbase" = "#BA68C8", "targetscan" = "#5C6BC0"
      )

      p_mirna_stacked_final_subtitle <- ggplot(
        top_20_stacked_data_by_total,
        aes(x = mature_mirna_id, y = unique_DEG_target_count, fill = database)
      ) +
        geom_bar(stat = "identity", position = "stack") +
        coord_flip() +
        labs(
          title = "Top 20 miRNA",
          subtitle = "by unique DEG targets",
          x = "miRNA",
          y = "Number of Unique DEG Targets",
          fill = "Database"
        ) +
        scale_fill_manual(values = new_light_palette, drop = FALSE) +
        theme_minimal(base_size = 12)

      ggsave(file.path(results_dir, "miRNA_Top20_Final_Plot.pdf"), plot = p_mirna_stacked_final_subtitle, width = 12, height = 8)
    }
  }
}


# #############################################################################
# SECTION 7: SESSION INFO
# #############################################################################

cat("Writing session info...\n")
capture.output(sessionInfo(), file = file.path(results_dir, "sessionInfo.txt"))

cat("\n--- Full Analysis Pipeline Finished. ---\n")
