import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, PandasTools
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from tqdm import tqdm 
import pandas as pd


class DataPreparation:
    """
    Data Preparation using Lipinski's rule of five and remove compounds including PAINS

    Parameters
    ----------
    data
        DataFrame containing the compound and its properties.
    
    Returns
    -------
    DataFrame
        Final processed DataFrame after filtering based on Lipinski's rule and removing PAINS.

    """

    def __init__(self, data):
        self.data = data

    def preprocess_and_filter(self):
        # Step 1: Data Preprocessing
        print("Step 1: Data Preprocessing")
        data = self.data
        invalid_rows = data[data["smiles"] == "smiles"].index
        data.drop(invalid_rows, axis=0, inplace=True)
        data = data.reset_index(drop=True)

        # Step 2: Lipinski Rule of Five
        print("Step 2: Applying Lipinski's Rule of Five")
        molecular_weight = []
        LogP = []
        HDonors = []
        HAcceptors = []

        for i in tqdm(range(len(data)), desc="Lipinski Processing..."):
            mol = Chem.MolFromSmiles(data["smiles"][i])

            desc_MolWt = Descriptors.MolWt(mol)
            molecular_weight.append(desc_MolWt)

            desc_MolLogP = Descriptors.MolLogP(mol)
            LogP.append(desc_MolLogP)

            desc_NumHDonors = Lipinski.NumHDonors(mol)
            HDonors.append(desc_NumHDonors)

            desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
            HAcceptors.append(desc_NumHAcceptors)

        molecular_weight = pd.DataFrame(molecular_weight, columns=["MolecularWeight"])
        LogP = pd.DataFrame(LogP, columns=["LogP"])
        HDonors = pd.DataFrame(HDonors, columns=["HDonors"])
        HAcceptors = pd.DataFrame(HAcceptors, columns=["HAcceptors"])

        data = pd.concat([data, molecular_weight, LogP, HDonors, HAcceptors], axis=1)
        data.to_csv("./data/Lipinski_descriptors.csv")

        for i in tqdm(range(len(data)), desc="Dropping..."):
            if data["MolecularWeight"][i] > 500:
                data = data.drop(i, axis=0)

            elif data["LogP"][i] > 5:
                data = data.drop(i, axis=0)

            elif data["HDonors"][i] > 5:
                data = data.drop(i, axis=0)

            elif data["HAcceptors"][i] > 10:
                data = data.drop(i, axis=0)

        data = data.reset_index(drop=True)
        data.to_csv("./data/Lipinski_results.csv")

        # Step 3: PAINS Filtering
        print("Step 3: Removing PAINS")
        # initialize filter
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog(params)

        # search for PAINS
        clean = []
        for index, row in tqdm(data.iterrows(), total=data.shape[0]):
            molecule = Chem.MolFromSmiles(row.smiles)
            entry = catalog.GetFirstMatch(molecule) # get the first matching PAINS
            if entry is None:
                # collect indices of molecules without PAINS
                clean.append(index)

        data = data.loc[clean]

        # Final Output
        data = data.reset_index(drop=True)
        data.to_csv("./data/final_filtered_data.csv", index=False)
        print("Processing completed! Filtered data saved.")
        return data
