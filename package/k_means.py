from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, scale
from sklearn.cluster import KMeans
from sklearn import metrics

from scipy.spatial.distance import cdist

import time
import random
from pathlib import Path
from tqdm import tqdm

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina 
from rdkit.Chem import Draw, rdFingerprintGenerator, PandasTools, Descriptors, Lipinski

class KMeansClustering():
    """
    KMeans Clustering for the compounds

    Parameters
    ----------
    data : DataFrame
        DataFrame containing the compound and its properties.
    query : str
        Query compound's smiles string
    
    Returns
    -------
    DataFrame
        DataFrame containing the compound and its properties.

    """

    def __init__(self, data, query):
        self.data = data
        self.query = query

    def normalize(self):
        data = self.data
        query = self.query

        # query compound's properties
        query_mol = Chem.MolFromSmiles(query)
        query_MolWt = Descriptors.MolWt(query_mol)
        query_MolLogP = Descriptors.MolLogP(query_mol)
        query_NumHDonors = Lipinski.NumHDonors(query_mol)
        query_NumHAcceptors = Lipinski.NumHAcceptors(query_mol)
        row = np.array([query_MolWt, query_MolLogP, query_NumHDonors, query_NumHAcceptors])
        row = pd.DataFrame(row.reshape(1, -1), columns=["MolecularWeight", "LogP", "HDonors", "HAcceptors"])
        data = pd.concat([data, row], axis=0)

        data = data.drop(["smiles", "zinc_id"], axis=1)
        normalized_data = data.values
        normalized_data = StandardScaler().fit_transform(normalized_data)

        feat_cols = [f"normalized feature {i}" for i in range(normalized_data.shape[1])]
        normalized_data = pd.DataFrame(normalized_data, columns=feat_cols)
        return normalized_data
    
    @staticmethod
    def ElbowMethod(data):
        sse = []

        kmeans_kawrgs = {
            "init": "random",
            "n_init": 10,
            "random_state": 42,
        }

        for k in tqdm((50, 100, 150, 200, 250, 300, 350, 400, 450, 500), desc="Processing..."):
            kmeans = KMeans(n_clusters=k, **kmeans_kawrgs)
            kmeans.fit(data)
            sse.append(kmeans.inertia_)

        plt.plot((50, 100, 150, 200, 250, 300, 350, 400, 450, 500), sse)
        plt.xticks((50, 100, 150, 200, 250, 300, 350, 400, 450, 500))
        plt.xlabel("Number of Clusters")
        plt.ylabel("SSE")
        plt.show()
        plt.savefig("./data/elbow_method.png")

    @classmethod
    def KMeans(cls, self, n_clusters):
        data = self.data

        kmeans_kawrgs = {
            "init": "random",
            "n_init": 10,
            "random_state": 42,
        }

        kmeans = KMeans(n_clusters=n_clusters, **kmeans_kawrgs)
        kmeans.fit(data)
        kmeans.labels_ = pd.DaraFrame(kmeans.labels_, columns=["Cluster"])
        final_data = pd.concat([data, kmeans.labels_], axis=1)
        return cls(final_data)