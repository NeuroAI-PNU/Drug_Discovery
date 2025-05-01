import rdkit
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)

import sys
import os
sys.path.append(os.path.join(os.environ["CONDA_PREFIX"],"share","RDKit","Contrib"))

from SA_Score import sascorer
from NP_Score import npscorer

import pandas as pd

def SAScorer(data):
    """
    Calculate synthetic accessibility score.
    
    Parameters
    ----------
    data : pandas.DataFrame
        Dataframe containing SMILES strings. 

    Returns
    -------
    float
        Synthetic accessibility score.
    """
    sascores = []
    for smiles in data["smiles"]:
        mol = Chem.MolFromSmiles(smiles)
        sascores.append(sascorer.calculateScore(mol))

    data = pd.concat([data, pd.DataFrame(sascores, columns=["sascore"])], axis=1)
    data.to_csv("./data/df_sascore.csv", index=False)
    return data