from pathlib import Path
import pandas as pd

basedir = Path(__file__).parent.parent

pathin = basedir / "nsp"
pathout = basedir / "data"

sub = pd.read_csv(pathin / "MMP_HCP_753_subs.txt", header=None) 
sub.columns = ['Subject']  
df = pd.read_csv(pathin / "nsp_balance_zscores.csv")
HCP = pd.read_csv(pathin / "HCP_YA_subjects_2025_11_23_16_00_16.csv")

df.insert(0, 'Subject', sub['Subject'].values)

hcp_cols = ["Subject", "HasGT", "ZygositySR", "ZygosityGT", "Family_ID", "Mother_ID", "Father_ID"]

df = pd.merge(df, HCP[hcp_cols], on="Subject", how="left")
df.to_csv(pathout / "nsp_balance_with_HCP_info_zscores.csv", index=False)

# Derive zyg from SR or GT
df["zyg"] = df["ZygositySR"].where(df["ZygositySR"].isin(["MZ","DZ"]), df["ZygosityGT"])
twins = df[df["zyg"].isin(["MZ","DZ"]) & df["Family_ID"].notna()].copy()

# Order twins within family
twins = twins.sort_values(["Family_ID","Subject"])
twins["twin"] = twins.groupby("Family_ID").cumcount() + 1
twins = twins[twins["twin"].isin([1,2])]

# Choose one phenotype to model (e.g., FC_rs)
wide = twins.pivot_table(index=["Family_ID","zyg"], columns="twin", values="FC_rs").reset_index()
wide = wide.dropna(subset=[1,2]).rename(columns={1:"FC_rs_T1", 2:"FC_rs_T2"})
wide.to_excel(pathout / "nsp_balance_with_HCP_info_zscores_wide_all_FC.xlsx", index=False)
#wide.to_csv(pathout / "twin_wide_FC_rs_zscores.csv", index=False)
"""
THIS NEEDS TO BE CHANGED TO SAVE WIDE DF IF NEEDED LATER. RIGHT NOW COLUMNS 
ARE NOT CREATED CORRECTLY FOR SAVING WIDE DF. 

"""


""""
#CREATE DF FOR MULTILEVEL REGRESSION
df_long = df.melt(
    id_vars=["Subject"], 
    value_vars=["FC_lang", "FC_motor", "FC_rs", "FC_social", "FC_wm"],
    var_name="Condition", 
    value_name="Balance"
)

# Clean up condition labels
df_long["Condition"] = df_long["Condition"].str.replace("FC_", "")

# Save the long-format dataset
df_long.to_csv(pathout / "nsp_balance_long.csv", index=False)
"""