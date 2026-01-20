from pathlib import Path
import scipy.io as scio
from nsp import hierarchichal_clustering
from nsp import segint_component
import pandas as pd
import numpy as np
from scipy import stats
from network import NestedSpectralPartition

basedir = Path(__file__).parent.parent

pathin = basedir / "data"
pathout = basedir / "nsp"

fc_files = [fc for fc in pathin.iterdir() if fc.suffix.lower() == '.mat']

FC_data = {}
#nsp_data = {}
nsp_balance = {}

for fc_file in fc_files:
    var_name = fc_file.stem

    mat_contents = scio.loadmat(fc_file)

    var_keys = [k for k in mat_contents.keys() if not k.startswith("__")]
    if not var_keys:
            raise ValueError(f"No variable found in {fc_file.name}") 

    FC_matrix = mat_contents[var_keys[0]]  # should be 400 x 400 x subjects

    FC_data[var_name] = FC_matrix
    print(f"{var_name} loaded with shape {FC_matrix.shape}")
    
    n_subjects = FC_data[var_name].shape[2]

    clus_size_key = var_name + "_clus_size"
    clus_num_key  = var_name + "_clus_num"

    #nsp_data[clus_size_key] = []
    #nsp_data[clus_num_key]  = []
    nsp_balance[var_name] = []

    for s in range(n_subjects):
        fc = FC_data[var_name][:,:,s] #aggiungi deepcoy()
        #np.fill_diagonal(fc,0)
        clus_size, clus_num = hierarchichal_clustering(fc)
        #nsp_data[clus_size_key].append(clus_size)
        #nsp_data[clus_num_key].append(clus_num)
        hint,hseg = segint_component(fc, clus_size, clus_num)
        hbal = hint - hseg
        nsp_balance[var_name].append(hbal)

# Standardize hbal to z-scores for each FC type across all subjects
for var_name in nsp_balance:
    values = np.array(nsp_balance[var_name])
    hbal_zscores = stats.zscore(values)
    nsp_balance[var_name] = hbal_zscores

df_balance = pd.DataFrame({k: pd.Series(v) for k, v in nsp_balance.items()})
csv_path = pathout / "nsp_balance_zscores.csv"
df_balance.to_csv(csv_path, index=False) 