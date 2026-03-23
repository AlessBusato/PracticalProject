require(umx)
library(readxl)
library(stringr)
library(ggplot2)
library(dplyr)

# Load data
twinData <- read_excel("add_your_path_here/ACE_twinData.xlsx")

# Identify all FC variables (e.g., FC_rs_T1, FC_rs_T2, ...)
fc_cols <- grep("^FC_.*_T[12]$", names(twinData), value = TRUE)

# Extract base phenotype names (FC_xxx)
phenotypes <- unique(str_replace(fc_cols, "_T[12]$", ""))

# Ensure numeric
twinData[fc_cols] <- lapply(twinData[fc_cols], function(x) as.numeric(trimws(x)))

# Subset by zygosity
dz <- subset(twinData, zyg == "DZ")
mz <- subset(twinData, zyg == "MZ")

results <- data.frame(
  phenotype = character(),
  A_std = numeric(),
  C_std = numeric(),
  E_std = numeric(),
  h2 = numeric(),
  stringsAsFactors = FALSE
)

# Loop over phenotypes
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
  
  # Extract standardized components
  A_std <- mod_ace$top$A_std$result[1, 1]
  C_std <- mod_ace$top$C_std$result[1, 1]
  E_std <- mod_ace$top$E_std$result[1, 1]
  
  # Heritability
  h2 <- A_std / (A_std + C_std + E_std)
  
  # Save results
  results <- rbind(results, data.frame(
    phenotype = pheno,
    A_std = A_std,
    C_std = C_std,
    E_std = E_std,
    h2 = h2
  ))
}

print(results)

# ---------------------------------------------------------
# RENAME PHENOTYPES FOR PLOTTING
# ---------------------------------------------------------

pretty_names <- c(
  "FC_motor" = "Balance Motor Task",
  "FC_gamb"  = "Balance Gambling Task",
  "FC_social" = "Balance social Task",
  "FC_lang" = "Balance Language Task",
  "FC_wm" = "Balance Working Memory Task",
  "FC_rs" = "Balance Resting State"
)

results$pretty_label <- pretty_names[results$phenotype]

# Reorder by h2 using the pretty labels
results <- results %>%
  arrange(h2) %>%
  mutate(pretty_label = factor(pretty_label, levels = pretty_label))

# ---------------------------------------------------------
# PLOT
# ---------------------------------------------------------

ggplot(results, aes(x = pretty_label, y = h2)) +
  geom_col(fill = "#4C72B0") +
  geom_text(aes(label = round(h2, 2)),
            vjust = -0.5, size = 3.5) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    x = "Endophenotype",
    y = "Heritability (h²)"
  ) +
  ylim(0, max(results$h2) + 0.1)

ggsave("add_your_path_here/heritability_plot.png", width = 8, height = 5, dpi = 300)