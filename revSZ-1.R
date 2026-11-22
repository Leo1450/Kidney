
# ============================================================
# DepMap gene expression extraction
# 17 candidate genes in 4 renal cancer cell lines
#
# Expression file:
# OmicsExpressionProteinCodingGenesTPMLogp1.csv
#
# Model metadata:
# Model.csv
#
# Expression values are used as provided by DepMap.
# No additional log transformation is performed.
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 1. Working directory and input files
# ------------------------------------------------------------

setwd("E:/Setareh")

expression_file <- "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
model_file      <- "Model.csv"

# ------------------------------------------------------------
# 2. Target genes and cell lines
# ------------------------------------------------------------

genes <- c(
  "SH2D2A",
  "VEGFA",
  "RAC2",
  "NFATC2",
  "PIK3R2",
  "SPHK1",
  "PLA2G4D",
  "JMJD7-PLA2G4B",
  "PLA2G4B",
  "PRKCG",
  "PLA2G4F",
  "SPHK2",
  "PLA2G4A",
  "PRKCA",
  "PLCG2",
  "PTGS2",
  "RAC3"
)

target_cell_lines <- c(
  "786-O",
  "769-P",
  "Caki-1",
  "OS-RC-2"
)

# ------------------------------------------------------------
# 3. Read Model.csv
# ------------------------------------------------------------

model <- fread(model_file)

required_columns <- c("ModelID", "CellLineName")

if (!all(required_columns %in% colnames(model))) {
  stop(
    "Model.csv must contain ModelID and CellLineName columns."
  )
}

model[, ModelID := as.character(ModelID)]
model[, CellLineName := trimws(as.character(CellLineName))]

# Standardized names for matching
model[, CellLineName_clean := toupper(CellLineName)]

target_clean <- toupper(trimws(target_cell_lines))

# ------------------------------------------------------------
# 4. Read DepMap expression matrix
# ------------------------------------------------------------

expr <- fread(expression_file)

# In your file, the first column is V1 and contains ModelIDs.
if (names(expr)[1] != "V1") {
  stop(
    "The first expression matrix column is not V1. ",
    "Please inspect the file structure."
  )
}

setnames(expr, "V1", "ModelID")

expr[, ModelID := as.character(ModelID)]

# Check ModelID uniqueness
if (anyDuplicated(expr$ModelID)) {
  stop("Duplicated ModelIDs found in expression matrix.")
}

cat("\nExpression matrix dimensions:\n")
print(dim(expr))

cat("\nFirst ModelIDs:\n")
print(head(expr$ModelID))

# ------------------------------------------------------------
# 5. Match the four cell lines to ModelIDs
# ------------------------------------------------------------

candidate_models <- model[
  CellLineName_clean %chin% target_clean
]

# Retain only ModelIDs present in the expression matrix
matched_models <- candidate_models[
  ModelID %chin% expr$ModelID
]

# Remove duplicate metadata rows for the same ModelID
matched_models <- unique(
  matched_models,
  by = c("CellLineName_clean", "ModelID")
)

# Count distinct ModelIDs per cell line
model_counts <- matched_models[
  ,
  .(Number_of_ModelIDs = uniqueN(ModelID)),
  by = "CellLineName_clean"
]

# Include all four requested cell lines in the report
cell_line_check <- data.table(
  CellLineName_clean = target_clean
)

cell_line_check <- merge(
  cell_line_check,
  model_counts,
  by = "CellLineName_clean",
  all.x = TRUE
)

cell_line_check[
  is.na(Number_of_ModelIDs),
  Number_of_ModelIDs := 0L
]

cat("\nCell line matching report:\n")
print(cell_line_check)

# Stop if matching is missing or ambiguous
if (any(cell_line_check$Number_of_ModelIDs != 1)) {
  
  cat("\nCandidate records in Model.csv:\n")
  
  report_columns <- intersect(
    c(
      "CellLineName",
      "ModelID",
      "DepMapModelType",
      "OncotreeSubtype",
      "Expression_Available"
    ),
    colnames(candidate_models)
  )
  
  print(candidate_models[, ..report_columns])
  
  stop(
    "At least one cell line has zero or multiple ModelIDs ",
    "in the expression matrix. Inspect the matching report."
  )
}

# Select the unique model for each cell line
cell_mapping <- matched_models[
  match(target_clean, CellLineName_clean)
]

# Store the requested cell line names in the desired order
cell_mapping[, RequestedCellLine := target_cell_lines]

# Show metadata
metadata_columns <- intersect(
  c(
    "RequestedCellLine",
    "CellLineName",
    "ModelID",
    "DepMapModelType",
    "OncotreeSubtype",
    "Expression_Available"
  ),
  colnames(cell_mapping)
)

cell_mapping_report <- cell_mapping[
  ,
  ..metadata_columns
]

cat("\nFinal cell line mapping:\n")
print(cell_mapping_report)

fwrite(
  cell_mapping_report,
  "DepMap_CellLine_ModelID_Mapping.csv"
)

# ------------------------------------------------------------
# 6. Identify gene symbols in expression matrix
# ------------------------------------------------------------

# Examples of original column names:
# TSPAN6 (7105)
# TNMD (64102)
# DPM1 (8813)

expression_columns <- setdiff(
  colnames(expr),
  "ModelID"
)

gene_column_map <- data.table(
  ExpressionColumn = expression_columns
)

# Remove the trailing Entrez ID in parentheses
gene_column_map[
  ,
  GeneSymbol := sub(
    "\\s*\\([^)]*\\)$",
    "",
    ExpressionColumn
  )
]

# Clean symbols
gene_column_map[
  ,
  GeneSymbol := trimws(GeneSymbol)
]

# Find columns corresponding to the target genes
target_gene_columns <- gene_column_map[
  GeneSymbol %chin% genes
]

cat("\nMatched gene columns:\n")
print(target_gene_columns)

# Check duplicate symbols
duplicate_symbols <- target_gene_columns[
  ,
  .N,
  by = "GeneSymbol"
][N > 1]

if (nrow(duplicate_symbols) > 0) {
  
  cat("\nDuplicated gene symbols:\n")
  print(duplicate_symbols)
  
  stop(
    "Duplicate gene symbols detected. ",
    "Resolve them before extracting expression values."
  )
}

# ------------------------------------------------------------
# 7. Gene availability report
# ------------------------------------------------------------

gene_availability <- data.table(
  Gene = genes
)

gene_availability[
  ,
  ExpressionColumn := gene_column_map[
    match(genes, GeneSymbol),
    ExpressionColumn
  ]
]

gene_availability[
  ,
  Gene_Found := !is.na(ExpressionColumn)
]

cat("\nGene availability:\n")
print(gene_availability)

cat(
  "\nNumber of genes found:",
  sum(gene_availability$Gene_Found),
  "out of",
  length(genes),
  "\n"
)

fwrite(
  gene_availability,
  "DepMap_17Genes_Availability.csv"
)

# ------------------------------------------------------------
# 8. Extract expression values
# ------------------------------------------------------------

# Match selected ModelIDs to expression matrix rows
selected_rows <- match(
  cell_mapping$ModelID,
  expr$ModelID
)

if (anyNA(selected_rows)) {
  stop("Some selected ModelIDs are missing from expression matrix.")
}

# Match genes to their expression columns
selected_columns <- gene_availability$ExpressionColumn

# Initialize matrix:
# rows = genes
# columns = cell lines

expression_matrix <- matrix(
  NA_real_,
  nrow = length(genes),
  ncol = length(target_cell_lines),
  dimnames = list(
    genes,
    target_cell_lines
  )
)

# Extract values only for genes found in the matrix
for (i in seq_along(genes)) {
  
  current_column <- selected_columns[i]
  
  if (!is.na(current_column)) {
    
    expression_matrix[i, ] <- as.numeric(
      expr[[current_column]][selected_rows]
    )
    
  }
}

# Convert to data.table
expression_wide <- as.data.table(
  expression_matrix,
  keep.rownames = "Gene"
)

# Display result
cat("\nFinal expression matrix:\n")
print(expression_wide)

# Save wide-format table
fwrite(
  expression_wide,
  "DepMap_17Genes_Expression_4_ccRCC.csv"
)

# ------------------------------------------------------------
# 9. Create long-format expression table
# ------------------------------------------------------------

expression_long <- melt(
  expression_wide,
  id.vars = "Gene",
  variable.name = "CellLineName",
  value.name = "Expression"
)

fwrite(
  expression_long,
  "DepMap_17Genes_LongFormat.csv"
)

# ------------------------------------------------------------
# 10. Cell line availability report
# ------------------------------------------------------------

cell_line_availability <- data.table(
  CellLineName = target_cell_lines,
  ModelID = cell_mapping$ModelID,
  Expression_data_found = cell_mapping$ModelID %chin% expr$ModelID
)

if ("DepMapModelType" %in% colnames(cell_mapping)) {
  cell_line_availability[
    ,
    DepMapModelType := cell_mapping$DepMapModelType
  ]
}

if ("OncotreeSubtype" %in% colnames(cell_mapping)) {
  cell_line_availability[
    ,
    OncotreeSubtype := cell_mapping$OncotreeSubtype
  ]
}

if ("Expression_Available" %in% colnames(cell_mapping)) {
  cell_line_availability[
    ,
    Expression_Available := cell_mapping$Expression_Available
  ]
}

fwrite(
  cell_line_availability,
  "DepMap_CellLine_Availability.csv"
)

# ------------------------------------------------------------
# 11. Final quality control
# ------------------------------------------------------------

cat("\n============================================\n")
cat("DEPMap extraction completed\n")
cat("============================================\n")

cat(
  "\nGenes detected:",
  sum(gene_availability$Gene_Found),
  "/",
  length(genes)
)

cat(
  "\nCell lines matched:",
  nrow(cell_mapping),
  "/",
  length(target_cell_lines)
)

cat(
  "\nMissing gene symbols:\n"
)

print(
  gene_availability[Gene_Found == FALSE, Gene]
)

cat("\nOutput files:\n")
cat("1. DepMap_17Genes_Expression_4_ccRCC.csv\n")
cat("2. DepMap_17Genes_LongFormat.csv\n")
cat("3. DepMap_17Genes_Availability.csv\n")
cat("4. DepMap_CellLine_ModelID_Mapping.csv\n")
cat("5. DepMap_CellLine_Availability.csv\n")




library(data.table)
library(ggplot2)

# ============================================================
# Heatmap of 17 candidate genes in 4 RCC cell lines
# Row-wise Z-score
# ============================================================

# Convert expression table to matrix
heatmap_matrix <- as.matrix(
  expression_wide[, -1, with = FALSE]
)

rownames(heatmap_matrix) <- expression_wide$Gene

# Remove genes with missing values
heatmap_matrix <- heatmap_matrix[
  complete.cases(heatmap_matrix),
  ,
  drop = FALSE
]

# Calculate row-wise Z-score
heatmap_zscore <- t(
  scale(t(heatmap_matrix))
)

# Remove rows with undefined Z-scores
heatmap_zscore <- heatmap_zscore[
  complete.cases(heatmap_zscore),
  ,
  drop = FALSE
]

# Convert to long format
heatmap_dt <- as.data.table(
  as.table(heatmap_zscore)
)

setnames(
  heatmap_dt,
  c("Gene", "CellLine", "Zscore")
)

heatmap_dt[, Gene := factor(
  Gene,
  levels = rev(rownames(heatmap_zscore))
)]

heatmap_dt[, CellLine := factor(
  CellLine,
  levels = colnames(heatmap_zscore)
)]

# Plot heatmap
p <- ggplot(
  heatmap_dt,
  aes(
    x = CellLine,
    y = Gene,
    fill = Zscore
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  
  geom_text(
    aes(label = sprintf("%.2f", Zscore)),
    size = 3.5
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Row Z-score"
  ) +
  
  labs(
    title = "Expression profile of candidate genes in RCC cell lines",
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      face = "italic"
    ),
    
    panel.grid = element_blank(),
    
    legend.position = "right"
  )

print(p)

# Save figures
ggsave(
  "DepMap_RCC_17Genes_Heatmap.pdf",
  p,
  width = 8,
  height = 7
)

ggsave(
  "DepMap_RCC_17Genes_Heatmap.png",
  p,
  width = 8,
  height = 7,
  dpi = 600
)
