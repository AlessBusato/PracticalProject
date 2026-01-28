require(umx)
library(readxl)
library(stringr)

# ------------------------------
# Load data
# ------------------------------
twinData <- read_excel("C:/Users/busat/Desktop/PP/prac_proj/data/nsp_balance_with_HCP_info_zscores_wide_all_FC.xlsx")

# Identify all FC variables (e.g., FC_rs_T1, FC_rs_T2, ...)
fc_cols <- grep("^FC_.*_T[12]$", names(twinData), value = TRUE)

# Extract base phenotype names (FC_xxx)
phenotypes <- unique(str_replace(fc_cols, "_T[12]$", ""))

# Ensure numeric
twinData[fc_cols] <- lapply(twinData[fc_cols], function(x) as.numeric(trimws(x)))

# Subset by zygosity
dz <- subset(twinData, zyg == "DZ")
mz <- subset(twinData, zyg == "MZ")

# ------------------------------
# Loop over phenotypes
# ------------------------------
for (pheno in phenotypes) {
  
  cat("\n=====================================\n")
  cat("ACE path coefficients for:", pheno, "\n")
  cat("=====================================\n")
  
  # Run ACE
  mod_ace <- try(umxACE(
    selDVs = pheno,
    sep    = "_T",
    dzData = dz,
    mzData = mz
  ), silent = TRUE)
  
  if (inherits(mod_ace, "try-error")) {
    cat("Model failed for", pheno, "\n")
    next
  }
  
  # Extract standardized A, C, E
  est <- umxSummary(mod_ace, returnStd = TRUE)
  
  a1 <- est$std$A[1,1]
  c1 <- est$std$C[1,1]
  e1 <- est$std$E[1,1]
  
  # Print in your desired format
  cat(sprintf("| %-10s | %6.3f | %6.3f | %6.3f |\n",
              pheno, a1, c1, e1))
}