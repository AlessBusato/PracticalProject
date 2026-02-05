
# Twin modeling: ACE and ADE

require(umx)
library(readxl)

# Load data
twinData <- read_excel("C:/Users/busat/Desktop/PP/prac_proj/data/twin_wide_FC_rs_zscores_corrected.xlsx")

# Ensure numeric
twinData$FC_rs_T1 <- as.numeric(trimws(twinData$FC_rs_T1))
twinData$FC_rs_T2 <- as.numeric(trimws(twinData$FC_rs_T2))

# Subset by zygosity
dz <- subset(twinData, zyg == "DZ")
mz <- subset(twinData, zyg == "MZ")

# ACE model
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

# ADE model
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

#models_ACE <- list(mod_ace, mod_ae, mod_ce, mod_e)
#models_ADE <- list(mod_ade, mod_ade_ae, mod_ade_e, mod_de)

models_final <- list(mod_ace, mod_ade)

# Compare models
summary(mod_ace)
summary(mod_ade)
umxCompare(models_final)


