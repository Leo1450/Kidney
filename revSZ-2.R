setwd("E:/Setareh/3/")

library(AnnotationDbi)
library(org.Hs.eg.db)

expr <- read.csv(
  "expression matrix(normalized-KIRC).csv",
  row.names = 1,
  check.names = FALSE
)

genes <- c(
  "SH2D2A", "PDCD1", "LAG3", "CTLA4",
  "CD3D", "CD3E", "PTPRC", "CD8A"
)

ensembl_ids <- sub("\\.[0-9]+$", "", rownames(expr))

mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = genes,
  keytype = "SYMBOL",
  columns = "ENSEMBL"
)

mapping <- unique(mapping[!is.na(mapping$ENSEMBL), ])

hits <- merge(
  data.frame(
    matrix_row = seq_len(nrow(expr)),
    ENSEMBL = ensembl_ids
  ),
  mapping,
  by = "ENSEMBL",
  sort = FALSE
)

hits <- hits[
  order(match(hits$SYMBOL, genes), hits$matrix_row),
]

missing_genes <- setdiff(genes, hits$SYMBOL)

if (length(missing_genes) > 0) {
  stop("Genes not found: ", paste(missing_genes, collapse = ", "))
}

out <- data.frame(
  Gene = hits$SYMBOL,
  Ensembl_ID = hits$ENSEMBL,
  expr[hits$matrix_row, , drop = FALSE],
  check.names = FALSE,
  row.names = NULL
)

write.csv(
  out,
  "KIRC_reviewer_genes_logTMM.csv",
  row.names = FALSE
)

table(out$Gene)
dim(out)
