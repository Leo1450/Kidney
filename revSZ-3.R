setwd("E:/Setareh/3/")

library(readxl)
library(dplyr)
library(survival)

# 1. Read input files
clinical_file <- "KIRC_Complete_Clinical_Data.xlsx"
expression_file <- "KIRC_reviewer_genes_logTMM.xlsx"

clinical <- as.data.frame(
  read_excel(clinical_file, sheet = "Clinical_Main")
)

expr <- as.data.frame(
  read_excel(expression_file)
)

genes <- c(
  "SH2D2A", "PDCD1", "LAG3", "CTLA4",
  "CD3D", "CD3E", "PTPRC", "CD8A"
)

markers <- c("CD3D", "CD3E", "PTPRC", "CD8A")

stopifnot(
  !anyDuplicated(expr$Gene),
  all(genes %in% expr$Gene)
)

# 2. Convert expression data to sample × gene format
sample_columns <- setdiff(
  colnames(expr),
  c("Gene", "Ensembl_ID")
)

expression_matrix <- as.matrix(
  expr[
    match(genes, expr$Gene),
    sample_columns,
    drop = FALSE
  ]
)

storage.mode(expression_matrix) <- "numeric"

expression_data <- as.data.frame(
  t(expression_matrix)
)

colnames(expression_data) <- genes
expression_data$Sample_ID <- rownames(expression_data)
rownames(expression_data) <- NULL

stopifnot(
  !anyDuplicated(expression_data$Sample_ID),
  !anyDuplicated(clinical$Sample_ID)
)

# 3. Match samples and retain primary tumors
stopifnot(
  all(expression_data$Sample_ID %in% clinical$Sample_ID)
)

merged <- inner_join(
  clinical,
  expression_data,
  by = "Sample_ID"
) %>%
  filter(shortLetterCode == "TP")

clinical_variables <- c(
  "ajcc_pathologic_stage",
  "tumor_grade",
  "vital_status",
  "days_to_death",
  "days_to_last_follow_up",
  "age_at_diagnosis",
  "sex_at_birth"
)

# Check clinical consistency across repeated tumor samples
clinical_check <- merged %>%
  group_by(patient) %>%
  summarise(
    across(
      all_of(clinical_variables),
      ~ n_distinct(.x, na.rm = FALSE)
    ),
    .groups = "drop"
  )

stopifnot(
  all(
    as.matrix(
      clinical_check[, clinical_variables, drop = FALSE]
    ) == 1
  )
)

patient_clinical <- merged %>%
  select(patient, all_of(clinical_variables)) %>%
  distinct()

# Average repeated profiles on the supplied log-expression scale
patient_data <- merged %>%
  group_by(patient) %>%
  summarise(
    across(all_of(genes), mean),
    .groups = "drop"
  ) %>%
  left_join(patient_clinical, by = "patient")

# 4. Define clinical variables
stage_map <- c(
  "Stage I" = 1,
  "Stage IA" = 1,
  "Stage IB" = 1,
  "Stage II" = 2,
  "Stage III" = 3,
  "Stage IIIB" = 3,
  "Stage IV" = 4
)

grade_map <- c(
  "G1" = 1,
  "G2" = 2,
  "G3" = 3,
  "G4" = 4
)

patient_data$stage <- unname(
  stage_map[patient_data$ajcc_pathologic_stage]
)

patient_data$grade <- unname(
  grade_map[patient_data$tumor_grade]
)

patient_data$age10 <- (
  as.numeric(patient_data$age_at_diagnosis) / 365.25 / 10
)

patient_data$male <- unname(
  c("male" = 1, "female" = 0)[patient_data$sex_at_birth]
)

patient_data$event <- unname(
  c("Dead" = 1, "Alive" = 0)[patient_data$vital_status]
)

patient_data$OS_days <- ifelse(
  patient_data$event == 1,
  as.numeric(patient_data$days_to_death),
  as.numeric(patient_data$days_to_last_follow_up)
)

# Standardize expression without re-normalization or another log transform
for (gene in genes) {
  patient_data[[paste0(gene, "_z")]] <- as.numeric(
    scale(patient_data[[gene]])
  )
}

# 5. Raw and adjusted checkpoint correlations
checkpoint_results <- lapply(
  c("PDCD1", "LAG3", "CTLA4"),
  function(target) {
    
    analysis_data <- patient_data[
      ,
      c("SH2D2A", target, markers),
      drop = FALSE
    ]
    
    analysis_data <- analysis_data[
      complete.cases(analysis_data),
      ,
      drop = FALSE
    ]
    
    model_SH2D2A <- lm(
      reformulate(markers, response = "SH2D2A"),
      data = analysis_data
    )
    
    model_checkpoint <- lm(
      reformulate(markers, response = target),
      data = analysis_data
    )
    
    partial_r <- cor(
      residuals(model_SH2D2A),
      residuals(model_checkpoint),
      method = "pearson"
    )
    
    n <- nrow(analysis_data)
    degrees_freedom <- n - length(markers) - 2
    
    t_statistic <- partial_r * sqrt(
      degrees_freedom / (1 - partial_r^2)
    )
    
    p_value <- 2 * pt(
      abs(t_statistic),
      df = degrees_freedom,
      lower.tail = FALSE
    )
    
    data.frame(
      Checkpoint = target,
      N = n,
      Raw_r = cor(
        analysis_data$SH2D2A,
        analysis_data[[target]],
        method = "pearson"
      ),
      Adjusted_r = partial_r,
      P = p_value
    )
  }
) %>%
  bind_rows()

checkpoint_results$FDR_BH <- p.adjust(
  checkpoint_results$P,
  method = "BH"
)

# 6. Adjusted logistic models: stage and grade
covariates <- c(
  "SH2D2A_z",
  paste0(markers, "_z"),
  "age10",
  "male"
)

clinical_results <- lapply(
  c("stage", "grade"),
  function(outcome) {
    
    analysis_data <- patient_data[
      ,
      c(outcome, covariates),
      drop = FALSE
    ]
    
    analysis_data <- analysis_data[
      complete.cases(analysis_data),
      ,
      drop = FALSE
    ]
    
    # Stage III-IV vs I-II; grade G3-G4 vs G1-G2
    analysis_data$high <- as.integer(
      analysis_data[[outcome]] > 2
    )
    
    model <- glm(
      reformulate(covariates, response = "high"),
      data = analysis_data,
      family = binomial()
    )
    
    beta <- unname(coef(model)["SH2D2A_z"])
    
    standard_error <- sqrt(
      vcov(model)["SH2D2A_z", "SH2D2A_z"]
    )
    
    saveRDS(
      model,
      paste0("KIRC_", outcome, "_model.rds")
    )
    
    data.frame(
      Outcome = outcome,
      N = nrow(analysis_data),
      High_group_N = sum(analysis_data$high),
      OR = exp(beta),
      Lower_95CI = exp(beta - 1.96 * standard_error),
      Upper_95CI = exp(beta + 1.96 * standard_error),
      P = coef(summary(model))[
        "SH2D2A_z",
        "Pr(>|z|)"
      ]
    )
  }
) %>%
  bind_rows()

# 7. Stage-stratified adjusted Cox model
survival_data <- patient_data[
  ,
  c("OS_days", "event", "stage", covariates),
  drop = FALSE
]

survival_data <- survival_data[
  complete.cases(survival_data),
  ,
  drop = FALSE
]

survival_data <- survival_data[
  survival_data$OS_days > 0,
  ,
  drop = FALSE
]

cox_formula <- as.formula(
  paste(
    "Surv(OS_days, event) ~",
    paste(
      c(covariates, "strata(stage)"),
      collapse = " + "
    )
  )
)

cox_model <- coxph(
  cox_formula,
  data = survival_data,
  ties = "efron",
  x = TRUE,
  y = TRUE
)

cox_beta <- unname(
  coef(cox_model)["SH2D2A_z"]
)

cox_standard_error <- sqrt(
  vcov(cox_model)["SH2D2A_z", "SH2D2A_z"]
)

survival_results <- data.frame(
  Outcome = "Overall survival",
  N = nrow(survival_data),
  Deaths = sum(survival_data$event),
  HR = exp(cox_beta),
  Lower_95CI = exp(
    cox_beta - 1.96 * cox_standard_error
  ),
  Upper_95CI = exp(
    cox_beta + 1.96 * cox_standard_error
  ),
  P = summary(cox_model)$coefficients[
    "SH2D2A_z",
    "Pr(>|z|)"
  ]
)

# 8. Proportional hazards diagnostic
PH_check <- cox.zph(cox_model)

# 9. Save tables, models and reproducibility information
write.csv(
  checkpoint_results,
  "KIRC_checkpoint_adjusted.csv",
  row.names = FALSE
)

write.csv(
  clinical_results,
  "KIRC_stage_grade_adjusted.csv",
  row.names = FALSE
)

write.csv(
  survival_results,
  "KIRC_OS_adjusted.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(PH_check$table),
  "KIRC_Cox_PH_check.csv"
)

write.csv(
  patient_data,
  "KIRC_patient_analysis_data.csv",
  row.names = FALSE
)

saveRDS(
  cox_model,
  "KIRC_OS_model.rds"
)

write.csv(
  data.frame(
    File = c(clinical_file, expression_file),
    MD5 = unname(
      tools::md5sum(
        c(clinical_file, expression_file)
      )
    )
  ),
  "KIRC_input_checksums.csv",
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = "KIRC_sessionInfo.txt"
)

# 10. Display main results
cat("\n--- Patient counts ---\n")
cat("Primary tumor samples:", nrow(merged), "\n")
cat("Unique patients:", nrow(patient_data), "\n")

cat("\n--- Checkpoint correlations ---\n")
print(checkpoint_results, row.names = FALSE)

cat("\n--- Stage and grade models ---\n")
print(clinical_results, row.names = FALSE)

cat("\n--- Survival model ---\n")
print(survival_results, row.names = FALSE)

cat("\n--- Proportional hazards check ---\n")
print(PH_check)
