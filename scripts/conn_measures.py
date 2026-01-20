from pathlib import Path
import scipy.io as scio
import bct

basedir = Path(__file__).parent.parent

pathin = basedir / "data"
pathout = basedir / "nsp"

fc_files = [fc for fc in pathin.iterdir() if fc.suffix.lower() == '.mat']

FC_data = {}

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
    #fc[fc < 0] = 0  # Set negative values to zero

    for s in range(n_subjects):
        fc = FC_data[var_name][:,:,s]
        
        betweeness = bct.efficiency_wei(fc)
        
        #DEVI METTERE TUTTI I VALORI NEGATIVI NELLA FC A ZERO
    break  # Remove this break to process all files