from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    PandasTools,
    Draw,
    Descriptors,
    MACCSkeys,
    rdFingerprintGenerator,)
from tqdm import tqdm

if __name__=="__main__":
    data = pd.read_csv("./data/final_data_added_lipinski.csv", low_memory=False)
    print("data load done")

    query_mol = Chem.MolFromSmiles("COc1ccc([C@@H]2CCCN2C(=O)Nc2cc(C(F)(F)F)ccc2N2CCOCC2)cc1")
    print("query_mol done")

    maccs_fp_qeury = MACCSkeys.GenMACCSKeys(query_mol)
    print("maccs_fp_qeury done")

    tanimoto_maccs_list = []
    for smile in tqdm(data["smiles"], desc="Tanimoto Similarity Calculation..."):
        mol = Chem.MolFromSmiles(smile)
        maccs_fp = MACCSkeys.GenMACCSKeys(mol)
        tanimoto_maccs = DataStructs.TanimotoSimilarity(maccs_fp_qeury, maccs_fp)
        tanimoto_maccs_list.append(tanimoto_maccs)

    data.loc[:, "tanimoto_maccs"] = tanimoto_maccs_list
    data.to_csv("added_tanimoto_maccs.csv")

    # I used MACCS fingerprint type