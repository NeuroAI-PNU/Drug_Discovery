from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

if __name__=="__main__":
    lipinski_data = pd.read_csv("./data/lipinski_data.csv", encoding="UTF-8", low_memory=False, index_col=0)
    print("DataFrame shape:", lipinski_data.shape)
    #lipinski_data.drop(columns=["MolWt", "LogP", "HDonors", "HAcceptors"], inplace=True)
    #print("DataFrame shape after dropping columns:", lipinski_data.shape)

    #### Filter for PAINS
    '''
    PAINS are compounds that often occur as hits in HTS even though they actually are false positives. 
    PAINS show activity at numerous targets rather than one specific target. 
    Such behavior results from unspecific binding or interaction with assay components.
    '''
    # initialize filter
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

    # search for PAINS
    matches = []
    clean = []
    for index, row in tqdm(lipinski_data.iterrows(), total=lipinski_data.shape[0]):
        molecule = Chem.MolFromSmiles(row.smiles)
        entry = catalog.GetFirstMatch(molecule) # Get the first matching PAINS
        if entry is None:
            # collet indices of molecules without PAINS
            clean.append(index)

    #matches = pd.DataFrame(matches)
    lipinski_data = lipinski_data.loc[clean]

    # save data without PAINS
    lipinski_data.to_csv("final_data_added_lipinski.csv", index=False)
    # save data with PAINS
    #matches.to_csv("matched_data_with_pains.csv", index=False)