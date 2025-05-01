import warnings
import subprocess
from pathlib import Path
import nglview as nv
from openbabel import pybel
from opencadd.structure.core import Structure
import time
import random
import pandas as pd
import numpy
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Draw
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import PandasTools
from tqdm import tqdm
import multiprocessing

# pdb -> pdbqt file translation start
#========================================================================================================================================================
#========================================================================================================================================================
#========================================================================================================================================================
def smiles_to_pdbqt(smiles, pdbqt_path, pH=7.4):
    """
    Convert a SMILES string to a PDBQT file needed by docking programs of the AutoDock family.

    Parameters
    ----------
    smiles: str
        SMILES string.
    pdbqt_path: str or pathlib.path
        Path to output PDBQT file.
    pH: float
        Protonation at given pH.
    """
    molecule = pybel.readstring("smi", smiles)
    # add hydrogens at given pH
    molecule.OBMol.CorrectForPH(pH)
    molecule.addh()
    # generate 3D coordinates
    molecule.make3D(forcefield="mmff94s", steps=10000)
    # add partial charges to each atom
    for atom in molecule.atoms:
        atom.OBAtom.GetPartialCharge()
    molecule.write("pdbqt", str(pdbqt_path), overwrite=True)
    return


# pdb -> pdbqt file translation done
#========================================================================================================================================================
#========================================================================================================================================================
#========================================================================================================================================================
# docking start
        
#### 4.Docking calculation ####
#### conda install -c conda-forge smina

def run_smina(
    ligand_path, protein_path, out_path, pocket_center, pocket_size, num_poses=10, exhaustiveness=8
):
    """
    Perform docking with Smina.

    Parameters
    ----------
    ligand_path: str or pathlib.Path
        Path to ligand PDBQT file that should be docked.
    protein_path: str or pathlib.Path
        Path to protein PDBQT file that should be docked to.
    out_path: str or pathlib.Path
        Path to which docking poses should be saved, SDF or PDB format.
    pocket_center: iterable of float or int
        Coordinates defining the center of the binding site.
    pocket_size: iterable of float or int
        Lengths of edges defining the binding site.
    num_poses: int
        Maximum number of poses to generate.
    exhaustiveness: int
        Accuracy of docking calculations.

    Returns
    -------
    output_text: str
        The output of the Smina calculation.
    """
    output_text = subprocess.check_output(
        [
            # In Max   #########
            #"./smina.osx.12",

            # In Linux #########
            "./smina.linux",
            
            "--ligand",
            str(ligand_path),
            "--receptor",
            str(protein_path),
            "--out",
            str(out_path),
            "--center_x",
            str(pocket_center[0]),
            "--center_y",
            str(pocket_center[1]),
            "--center_z",
            str(pocket_center[2]),
            "--size_x",
            str(pocket_size[0]),
            "--size_y",
            str(pocket_size[1]),
            "--size_z",
            str(pocket_size[2]),
            "--num_modes",
            str(num_poses),
            "--exhaustiveness",
            str(exhaustiveness),
        ],
        universal_newlines=True,  # needed to capture output text
    )
    return output_text


# docking done
#========================================================================================================================================================
#========================================================================================================================================================
#========================================================================================================================================================

if __name__ == "__main__":

    warnings.filterwarnings("ignore")
    ob_log_handler = pybel.ob.OBMessageHandler()
    pybel.ob.obErrorLog.SetOutputLevel(0)
    print("warining setting is done")

    pdb_id = "2jkq"
    #structure = Structure.from_pdbid(pdb_id) # url error로 작동하지않아 로컬 파일 불러와서 진행

    pdb_file = "./data/2jkq.pdb"
    structure = Structure(pdb_file)
    # element information maybe missing, but important for subsequent pdbqt conversion
    if not hasattr(structure.atoms, "elements"):
        structure.add_TopologyAttr("elements", structure.atoms.types)
    structure
    print("structure load done")

    #data = pd.read_csv("./data/rank_top_10_molecules(only_smiles).csv", low_memory=False)
    #print("data load done")

    # Step 1, pdb -> pdbqt file translation using multiprocessing.Process
    #====================================================================================================================================================
    # 한번 코드 진행했고 pdbqt 파일 생성 완료되어 주석처리함
    #for i in range(len(data)):
    #    smiles_to_pdbqt(data["smiles"][i], f"./data/rank10_pdbqt_data/rank_{i}.pdbqt")

    #print("pdbqt files translation done")
    #====================================================================================================================================================
    #====================================================================================================================================================
    #====================================================================================================================================================

    ligand_resname ="VG8"
    ligand = structure.select_atoms(f"resname {ligand_resname}")
    pocket_center = (ligand.positions.max(axis=0) + ligand.positions.min(axis=0)) / 2
    pocket_size = ligand.positions.max(axis=0) - ligand.positions.min(axis=0) + 5
    print(pocket_center, pocket_size)

    # Step 2, docking calculation using multiprocessing.Process
    
    docking = open("docking_score_only_ZN27.txt", mode="w")
    results = run_smina(
        "./data/pdbqt_file/ligand_Zn27.pdbqt",
        "./data/rank10_pdbqt_data/protein_4.pdbqt",
        "docking_poses_Zn27.pdb",
        pocket_center,
        pocket_size,
    )
    docking.write(results)
    docking.close()

    print("docking done")
    #====================================================================================================================================================
    #====================================================================================================================================================
    #====================================================================================================================================================