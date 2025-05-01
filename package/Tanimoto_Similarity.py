import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import (PandasTools, Draw, Descriptors, MACCSkeys, rdFingerprintGenerator)
from tqdm import tqdm

def TanimotoSimilarity(data, query):
    """
    Tanimoto Similarity Calculation using MACCS fingerprint
    
    Parameters
    ----------
    data
        DataFrame containing the compound and its properties.
    query : str
        SMILES string of query compound.

    Returns
    -------
    DataFrame
        DataFrame containing the compound and its properties with tanimoto similarity.
    """

    query_mol = Chem.MolFromSmiles(query)
    maccs_fp_query = MACCSkeys.GenMACCSKeys(query_mol)

    tanimoto_maccs_list = []
    for smile in tqdm(data["smiles"], desc="Tanimoto Similarity Calculation..."):
        mol = Chem.MolFromSmiles(smile)
        maccs_fp = MACCSkeys.GenMACCSKeys(mol)
        tanimoto_maccs = DataStructs.TanimotoSimilarity(maccs_fp_query, maccs_fp)
        tanimoto_maccs_list.append(tanimoto_maccs)

    data.loc[:, "tanimoto_maccs"] = tanimoto_maccs_list
    data.to_csv("./data/Tanimoto_Similiarty_results.csv")
    return data