from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    PandasTools,
    Draw,
    Descriptors,
    MACCSkeys,
    rdFingerprintGenerator,
)
from tqdm import tqdm

if __name__=="__main__":
    maccs = pd.read_csv("./data/added_tanimoto_maccs.csv", usecols=["smiles", "zinc_id", "tanimoto_maccs"], low_memory=False)
    maccs.head()
    print("data loading done")

    list_1 = []
    list_2 = []
    list_3 = []
    list_4 = []
    list_5 = []
    extra_data = []

    for i in tqdm(range(len(maccs)), desc="Processing..."):
        if 0 <= maccs["tanimoto_maccs"][i] < 0.2:
            list_1.append(maccs.iloc[i])
        elif 0.2 <= maccs["tanimoto_maccs"][i] < 0.4:
            list_2.append(maccs.iloc[i])
        elif 0.4 <= maccs["tanimoto_maccs"][i] < 0.6:
            list_3.append(maccs.iloc[i])
        elif 0.6 <= maccs["tanimoto_maccs"][i] < 0.8:
            list_4.append(maccs.iloc[i])
        elif 0.8 <= maccs["tanimoto_maccs"][i] <= 1:
            list_5.append(maccs.iloc[i])
        else:
            extra_data.append(maccs.iloc[i])

    print(len(list_1))
    print(len(list_2))
    print(len(list_3))
    print(len(list_4))
    print(len(list_5))
    print(len(extra_data))

    list_5 = pd.DataFrame(list_5)
    list_5 = list_5.reset_index(drop=True)

    list_5.to_csv("./data/similarity_with_zn27_by_tanimoto_maccs.csv", index=False)

