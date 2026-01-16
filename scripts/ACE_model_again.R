require(umx)
library(readxl)

# Load your wide-format twin dataset
twinData <- read_excel("C:/Users/busat/Desktop/PP/practical_project_2/data/twin_wide_FC_rs_zscores_corrected.xlsx")

# Treat DV as numeric
twinData$FC_rs_T1 <- as.numeric(trimws(twinData$FC_rs_T1))
twinData$FC_rs_T2 <- as.numeric(trimws(twinData$FC_rs_T2))

# Subset by zygosity
dz <- subset(twinData, zyg == "DZ")
mz <- subset(twinData, zyg == "MZ")

# --- Fit full ACE model
mod_ace <- umxACE(
  selDVs = "FC_rs",
  sep    = "_T",
  dzData = dz,
  mzData = mz
)

# --- Fit AE model (drop C)
mod_ae <- umxModify(mod_ace, update = "c_r1c1", name = "AE")

# --- Fit CE model (drop A)
mod_ce <- umxModify(mod_ace, update = "a_r1c1", name = "CE")

# --- Fit E model (drop A and C)
mod_e <- umxModify(mod_ae, update = "a_r1c1", name = "E")

# --- Compare models
umxCompare(mod_ace, c(mod_ae, mod_ce, mod_e))

# --- Summary of each model
summary(mod_ace)
summary(mod_ae)
summary(mod_ce)
summary(mod_e)

# --- Plot ACE model (standardized paths)
umxPlotACE(
  mod_ace,
  file = "ACE_FC_rs",
  std = TRUE,
  means = FALSE,
  digits = 2,
  strip_zero = FALSE
)
