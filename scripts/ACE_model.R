# ==============================
# Twin modeling: ACE and ADE
# ==============================

require(umx)
library(readxl)

# ------------------------------
# Load data
# ------------------------------
twinData <- read_excel("C:/Users/busat/Desktop/PP/prac_proj/data/twin_wide_FC_rs_zscores_corrected.xlsx")

# Ensure numeric
twinData$FC_rs_T1 <- as.numeric(trimws(twinData$FC_rs_T1))
twinData$FC_rs_T2 <- as.numeric(trimws(twinData$FC_rs_T2))

# Subset by zygosity
dz <- subset(twinData, zyg == "DZ")
mz <- subset(twinData, zyg == "MZ")

# ------------------------------
# --- ACE model
# ------------------------------
mod_ace <- umxACE(
  selDVs = "FC_rs",
  sep    = "_T",
  dzData = dz,
  mzData = mz
)

# AE (drop C)
mod_ae <- umxModify(mod_ace, update = "c_r1c1", name = "AE")

# CE (drop A)
mod_ce <- umxModify(mod_ace, update = "a_r1c1", name = "CE")

# E (drop A and C)
mod_e <- umxModify(mod_ae, update = "a_r1c1", name = "E")

# ------------------------------
# --- ADE model
# ------------------------------
mod_ade <- umxACE(
  selDVs = "FC_rs",
  sep    = "_T",
  dzData = dz,
  mzData = mz,
  dzCr   = .25   # dominance correlation for DZ
)

# DE (drop A)
mod_de <- umxModify(mod_ade, update = "a_r1c1", name = "DE")

# AE for ADE (drop D)
mod_ade_ae <- umxModify(mod_ade, update = "c_r1c1", name = "ADE_AE")

# E (drop A and D)
mod_ade_e <- umxModify(mod_ade_ae, update = "a_r1c1", name = "ADE_E")

# ------------------------------
# --- Compare all models
# ------------------------------
all_models <- list(mod_ace, mod_ae, mod_ce, mod_e,
                   mod_ade, mod_ade_ae, mod_de, mod_ade_e)

umxCompare(all_models)

# ------------------------------
# --- Extract standardized variance components from best-fitting model
# Replace 'best_model' with whichever fits best according to AIC/likelihood
# ------------------------------
best_model <- mod_ade  # <- example: ADE model is best

# Extract paths
a <- mxEval(top.a, best_model)
d <- mxEval(top.c, best_model)  # D is stored in top.c
e <- mxEval(top.e, best_model)

# Compute variance
A_var <- a^2
D_var <- d^2
E_var <- e^2

total <- A_var + D_var + E_var

# Standardized proportions
var_props <- round(c(
  A = A_var / total,
  D = D_var / total,
  E = E_var / total
), 3)

print("Proportion of variance explained by best-fitting model:")
print(var_props)

# ------------------------------
# --- Optional: plot paths for best model
# ------------------------------
umxPlotACE(
  best_model,
  file = "BestModel_FC_rs",
  std = TRUE,
  means = FALSE,
  digits = 2,
  strip_zero = FALSE
)
