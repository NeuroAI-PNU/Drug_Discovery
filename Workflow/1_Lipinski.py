import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from tqdm import tqdm 
import pandas as pd

def data_preprocessing(data_df):
    invalid_rows = data_df[data_df["smiles"]=="smiles"].index
    data_df.drop(invalid_rows, axis=0, inplace=True)
    data_df = data_df.reset_index(drop=True)
    return data_df

def Molecular_Weight(self, verbose=False):

    Mt_list = []
    for i in tqdm(range(len(self.data_df)), desc="Molecular Weight Processing"):
        mol=Chem.MolFromSmiles(self.data_df["smiles"][i])

        desc_MolWt = Descriptors.MolWt(mol)
        Mt_list.append(desc_MolWt)

    Mt_list = pd.DataFrame(Mt_list, columns=["MolWt"])
    self.data_df = pd.concat([self.data_df, Mt_list], axis=1)

    self.data_df.to_csv("Molecular_Weight.csv")

def LogP(self, verbose=False):

    LogP_list = []
    for i in tqdm(range(len(self.data_df)), desc="LogP Processing..."):
        mol=Chem.MolFromSmiles(self.data_df["smiles"][i])

        desc_MolLogP = Descriptors.MolLogP(mol)
        LogP_list.append(desc_MolLogP)

    LogP_list = pd.DataFrame(LogP_list, columns=["LogP"])

    self.data_df = pd.concat([self.data_df, LogP_list], axis=1)

    self.data_df.to_csv("LogP.csv")

def HDonor(self, verbose=False):

    HDonors_list = []
    for i in tqdm(range(len(self.data_df)), desc="HDonors Processing..."):
        mol=Chem.MolFromSmiles(self.data_df["smiles"][i])
            
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        HDonors_list.append(desc_NumHDonors)

    HDonors_list = pd.DataFrame(HDonors_list, columns=["HDonors"])

    self.data_df = pd.concat([self.data_df, HDonors_list], axis=1)

    self.data_df.to_csv("HDonors.csv")

def HAcceptor(self, verbose=False):

    HAcceptors_list = []
    for i in tqdm(range(len(self.data_df)), desc="HAcceptors Processing..."):
        mol=Chem.MolFromSmiles(self.data_df["smiles"][i])

        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
        HAcceptors_list.append(desc_NumHAcceptors)

    HAcceptors_list = pd.DataFrame(HAcceptors_list, columns=["HAcceptors"])

    self.data_df = pd.concat([self.data_df, HAcceptors_list], axis=1)

    self.data_df.to_csv("HAcceptors.csv")

if __name__=="__main__":
    data = pd.read_csv("merged_zinc_dataset.txt", sep="\t", encoding="UTF-8", usecols=["zinc_id", "smiles"], low_memory=False)
    data_df = data_preprocessing(data)
    Molecular_Weight(data_df)
    LogP(data_df)
    HDonor(data_df)
    HAcceptor(data_df)
