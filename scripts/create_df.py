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
#df.to_csv(pathout / "nsp_balance_with_HCP_info_zscores.csv", index=False)

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
wide.to_excel(pathout / "balance_rs.xlsx", index=False)
