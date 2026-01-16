from pathlib import Path
import scipy.io as scio
import pickle
from nsp import hierarchichal_clustering, segint_component

def save_nsp(nsp, filename="nsp_data.pkl", path=None):

    if path is None:
        path = Path(".")  # default: current folder
    else:
        path = Path(path)

    path.mkdir(parents=True, exist_ok=True)  # create folder if it doesn't exist
    full_path = path / filename

    with open(full_path, "wb") as f:
        pickle.dump(nsp, f)

    print(f"nsp_data saved to {full_path}")

def main():
    basedir = Path(__file__).parent.parent

    pathin = basedir / "data"
    pathout = basedir / "nsp"

    fc_files = [fc for fc in pathin.iterdir() if fc.suffix.lower() == '.mat']

    FC_data = {}
    nsp_data = {}

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

        nsp_data[clus_size_key] = []
        nsp_data[clus_num_key]  = []

        for s in range(n_subjects):
            fc = FC_data[var_name][:,:,s]
            clus_size, clus_num = hierarchichal_clustering(fc)
            hseg,hint = segint_component(fc, clus_size, clus_num)
            hbal
            nsp_data[clus_size_key].append(clus_size)
            nsp_data[clus_num_key].append(clus_num)
        
    save_nsp(nsp_data, filename="nsp_data.pkl", path=pathout)

if __name__ == "__main__":
    main()



        