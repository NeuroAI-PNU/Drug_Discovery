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

def multiprocessing_smiles_to_pdbqt_0(data):
    for i in tqdm(range(130), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_1(data):
    for i in tqdm(range(130, 260), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_2(data):
    for i in tqdm(range(260, 390), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_3(data):
    for i in tqdm(range(390, 520), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_4(data):
    for i in tqdm(range(520, 650), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_5(data):
    for i in tqdm(range(650, 780), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_6(data):
    for i in tqdm(range(780, 910), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_7(data):
    for i in tqdm(range(910, 1040), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_8(data):
    for i in tqdm(range(1040, 1170), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_9(data):
    for i in tqdm(range(1170, 1300), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smils_to_pdbqt_10(data):
    for i in tqdm(range(1300, 1430), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smils_to_pdbqt_11(data):
    for i in tqdm(range(1430, 1560), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smils_to_pdbqt_12(data):
    for i in tqdm(range(1560, 1690), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smils_to_pdbqt_13(data):
    for i in tqdm(range(1690, 1820), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smils_to_pdbqt_14(data):
    for i in tqdm(range(1820, 1950), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smils_to_pdbqt_15(data):
    for i in tqdm(range(1950, 2080), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smils_to_pdbqt_16(data):
    for i in tqdm(range(2080, 2210), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smils_to_pdbqt_17(data):
    for i in tqdm(range(2210, 2340), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smils_to_pdbqt_18(data):
    for i in tqdm(range(2340, 2470), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smils_to_pdbqt_19(data):
    for i in tqdm(range(2470, 2600), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smils_to_pdbqt_20(data):
    for i in tqdm(range(2600, 2730), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_21(data):
    for i in tqdm(range(2730, 2860), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_22(data):
    for i in tqdm(range(2860, 2990), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_23(data):
    for i in tqdm(range(2990, 3120), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_24(data):
    for i in tqdm(range(3120, 3250), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_25(data):
    for i in tqdm(range(3250, 3380), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_26(data):
    for i in tqdm(range(3380, 3510), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_27(data):
    for i in tqdm(range(3510, 3640), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_28(data):
    for i in tqdm(range(3640, 3770), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_29(data):
    for i in tqdm(range(3770, 3900), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_30(data):
    for i in tqdm(range(3900, 4030), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_31(data):
    for i in tqdm(range(4030, 4160), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_32(data):
    for i in tqdm(range(4160, 4290), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_33(data):
    for i in tqdm(range(4290, 4420), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_34(data):
    for i in tqdm(range(4420, 4550), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_35(data):
    for i in tqdm(range(4550, 4680), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_36(data):
    for i in tqdm(range(4680, 4810), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_37(data):
    for i in tqdm(range(4810, 4940), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_38(data):
    for i in tqdm(range(4940, 5070), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

def multiprocessing_smiles_to_pdbqt_39(data):
    for i in tqdm(range(5070, 5159), desc="pdbqt files translation..."):
        smiles_to_pdbqt(data["smiles"][i], f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt")

# pdb -> pdbqt file translation done
#========================================================================================================================================================
#========================================================================================================================================================
#========================================================================================================================================================
# docking start
        
#### 4.Docking calculation ####
#### conda install -c conda-forge smina

def run_smina(
    ligand_path, protein_path, pocket_center, pocket_size, num_poses=10, exhaustiveness=8
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
            #"--out",
            #str(out_path),
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

def multiprocessing_run_smina_0():
    docking_0 = open("docking_score_0.txt", mode="w")
    for i in tqdm(range(130), desc="docking__0..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_0.write(results)
    docking_0.close()

def multiprocessing_run_smina_1():
    docking_1 = open("docking_score_1.txt", mode="w")
    for i in tqdm(range(130, 260), desc="docking__1..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_1.write(results)
    docking_1.close()

def multiprocessing_run_smina_2():
    docking_2 = open("docking_score_2.txt", mode="w")
    for i in tqdm(range(260, 390), desc="docking__2..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_2.write(results)
    docking_2.close()

def multiprocessing_run_smina_3():
    docking_3 = open("docking_score_3.txt", mode="w")
    for i in tqdm(range(390, 520), desc="docking__3..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_3.write(results)
    docking_3.close()

def multiprocessing_run_smina_4():
    docking_4 = open("docking_score_4.txt", mode="w")
    for i in tqdm(range(520, 650), desc="docking__4..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_4.write(results)
    docking_4.close()

def multiprocessing_run_smina_5():
    docking_5 = open("docking_score_5.txt", mode="w")
    for i in tqdm(range(650, 780), desc="docking__5..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_5.write(results)
    docking_5.close()

def multiprocessing_run_smina_6():
    docking_6 = open("docking_score_6.txt", mode="w")
    for i in tqdm(range(780, 910), desc="docking__6..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_6.write(results)
    docking_6.close()

def multiprocessing_run_smina_7():
    docking_7 = open("docking_score_7.txt", mode="w")
    for i in tqdm(range(910, 1040), desc="docking__7..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_7.write(results)
    docking_7.close()

def multiprocessing_run_smina_8():
    docking_8 = open("docking_score_8.txt", mode="w")
    for i in tqdm(range(1040, 1170), desc="docking__8..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_8.write(results)
    docking_8.close()

def multiprocessing_run_smina_9():
    docking_9 = open("docking_score_9.txt", mode="w")
    for i in tqdm(range(1170, 1300), desc="docking__9..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_9.write(results)
    docking_9.close()

def multiprocessing_run_smina_10():
    docking_10 = open("docking_score_10.txt", mode="w")
    for i in tqdm(range(1300, 1430), desc="docking__10..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_10.write(results)
    docking_10.close()

def multiprocessing_run_smina_11():
    docking_11 = open("docking_score_11.txt", mode="w")
    for i in tqdm(range(1430, 1560), desc="docking__11..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_11.write(results)
    docking_11.close()

def multiprocessing_run_smina_12():
    docking_12 = open("docking_score_12.txt", mode="w")
    for i in tqdm(range(1560, 1690), desc="docking__12..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_12.write(results)
    docking_12.close()

def multiprocessing_run_smina_13():
    docking_13 = open("docking_score_13.txt", mode="w")
    for i in tqdm(range(1690, 1820), desc="docking__13..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_13.write(results)
    docking_13.close()

def multiprocessing_run_smina_14():
    docking_14 = open("docking_score_14.txt", mode="w")
    for i in tqdm(range(1820, 1950), desc="docking__14..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_14.write(results)
    docking_14.close()

def multiprocessing_run_smina_15():
    docking_15 = open("docking_score_15.txt", mode="w")
    for i in tqdm(range(1950, 2080), desc="docking__15..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_15.write(results)
    docking_15.close()

def multiprocessing_run_smina_16():
    docking_16 = open("docking_score_16.txt", mode="w")
    for i in tqdm(range(2080, 2210), desc="docking__16..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_16.write(results)
    docking_16.close()

def multiprocessing_run_smina_17():
    docking_17 = open("docking_score_17.txt", mode="w")
    for i in tqdm(range(2210, 2340), desc="docking__17..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_17.write(results)
    docking_17.close()

def multiprocessing_run_smina_18():
    docking_18 = open("docking_score_18.txt", mode="w")
    for i in tqdm(range(2340, 2470), desc="docking__18..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_18.write(results)
    docking_18.close()

def multiprocessing_run_smina_19():
    docking_19 = open("docking_score_19.txt", mode="w")
    for i in tqdm(range(2470, 2600), desc="docking__19..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_19.write(results)
    docking_19.close()

def multiprocessing_run_smina_20():
    docking_20 = open("docking_score_20.txt", mode="w")
    for i in tqdm(range(2600, 2730), desc="docking__20..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_20.write(results)
    docking_20.close()

def multiprocessing_run_smina_21():
    docking_21 = open("docking_score_21.txt", mode="w")
    for i in tqdm(range(2730, 2860), desc="docking__21..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_21.write(results)
    docking_21.close()

def multiprocessing_run_smina_22():
    docking_22 = open("docking_score_22.txt", mode="w")
    for i in tqdm(range(2860, 2990), desc="docking__22..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_22.write(results)
    docking_22.close()

def multiprocessing_run_smina_23():
    docking_23 = open("docking_score_23.txt", mode="w")
    for i in tqdm(range(2990, 3120), desc="docking__23..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_23.write(results)
    docking_23.close()

def multiprocessing_run_smina_24():
    docking_24 = open("docking_score_24.txt", mode="w")
    for i in tqdm(range(3120, 3250), desc="docking__24..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_24.write(results)
    docking_24.close()

def multiprocessing_run_smina_25():
    docking_25 = open("docking_score_25.txt", mode="w")
    for i in tqdm(range(3250, 3380), desc="docking__25..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_25.write(results)
    docking_25.close()

def multiprocessing_run_smina_26():
    docking_26 = open("docking_score_26.txt", mode="w")
    for i in tqdm(range(3380, 3510), desc="docking__26..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_26.write(results)
    docking_26.close()

def multiprocessing_run_smina_27():
    docking_27 = open("docking_score_27.txt", mode="w")
    for i in tqdm(range(3510, 3640), desc="docking__27..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_27.write(results)
    docking_27.close()

def multiprocessing_run_smina_28():
    docking_28 = open("docking_score_28.txt", mode="w")
    for i in tqdm(range(3640, 3770), desc="docking__28..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_28.write(results)
    docking_28.close()

def multiprocessing_run_smina_29():
    docking_29 = open("docking_score_29.txt", mode="w")
    for i in tqdm(range(3770, 3900), desc="docking__29..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_29.write(results)
    docking_29.close()

def multiprocessing_run_smina_30():
    docking_30 = open("docking_score_30.txt", mode="w")
    for i in tqdm(range(3900, 4030), desc="docking__30..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_30.write(results)
    docking_30.close()

def multiprocessing_run_smina_31():
    docking_31 = open("docking_score_31.txt", mode="w")
    for i in tqdm(range(4030, 4160), desc="docking__31..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_31.write(results)
    docking_31.close()

def multiprocessing_run_smina_32():
    docking_32 = open("docking_score_32.txt", mode="w")
    for i in tqdm(range(4160, 4290), desc="docking__32..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_32.write(results)
    docking_32.close()

def multiprocessing_run_smina_33():
    docking_33 = open("docking_score_33.txt", mode="w")
    for i in tqdm(range(4290, 4420), desc="docking__33..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_33.write(results)
    docking_33.close()

def multiprocessing_run_smina_34():
    docking_34 = open("docking_score_34.txt", mode="w")
    for i in tqdm(range(4420, 4550), desc="docking__34..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_34.write(results)
    docking_34.close()

def multiprocessing_run_smina_35():
    docking_35 = open("docking_score_35.txt", mode="w")
    for i in tqdm(range(4550, 4680), desc="docking__35..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_35.write(results)
    docking_35.close()

def multiprocessing_run_smina_36():
    docking_36 = open("docking_score_36.txt", mode="w")
    for i in tqdm(range(4680, 4810), desc="docking__36..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_36.write(results)
    docking_36.close()

def multiprocessing_run_smina_37():
    docking_37 = open("docking_score_37.txt", mode="w")
    for i in tqdm(range(4810, 4940), desc="docking__37..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_37.write(results)
    docking_37.close()

def multiprocessing_run_smina_38():
    docking_38 = open("docking_score_38.txt", mode="w")
    for i in tqdm(range(4940, 5070), desc="docking__38..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_38.write(results)
    docking_38.close()

def multiprocessing_run_smina_39():
    docking_39 = open("docking_score_39.txt", mode="w")
    for i in tqdm(range(5070, 5159), desc="docking__39..."):
        results = run_smina(
            f"./data/pdbqt_data/ligand_based_tanimoto_similarity/ligand_{i}.pdbqt",
            "./data/pdbqt_data/protein_4.pdbqt",
            #f"docking_poses_{i}",
            pocket_center,
            pocket_size,
        )
        docking_39.write(results)
    docking_39.close()

# docking done
#========================================================================================================================================================
#========================================================================================================================================================
#========================================================================================================================================================

if __name__ == "__main__":

    warnings.filterwarnings("ignore")
    ob_log_handler = pybel.ob.OBMessageHandler()
    pybel.ob.obErrorLog.SetOutputLevel(0)
    print("warining setting is done")

    # Structure.from_pdbid()가 작동안하면 local에서 pdb 파일 불러와서 진행하는 방법도 있음
    pdb_id = "2jkq"
    structure = Structure.from_pdbid(pdb_id)
    # element information maybe missing, but important for subsequent pdbqt conversion
    if not hasattr(structure.atoms, "elements"):
        structure.add_TopologyAttr("elements", structure.atoms.types)
    structure
    print("structure load done")

    data = pd.read_csv("./data/added_tanimoto_maccs.csv", low_memory=False)
    print("data load done")

    # Step 1, pdb -> pdbqt file translation using multiprocessing.Process
    #====================================================================================================================================================
    procs = []

    p0 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_0, args=(data,))
    p0.start()
    procs.append(p0)

    p1 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_1, args=(data,))
    p1.start()
    procs.append(p1)

    p2 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_2, args=(data,))
    p2.start()
    procs.append(p2)

    p3 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_3, args=(data,))
    p3.start()
    procs.append(p3)

    p4 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_4, args=(data,))
    p4.start()
    procs.append(p4)

    p5 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_5, args=(data,))
    p5.start()
    procs.append(p5)

    p6 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_6, args=(data,))
    p6.start()
    procs.append(p6)

    p7 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_7, args=(data,))
    p7.start()
    procs.append(p7)

    p8 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_8, args=(data,))
    p8.start()
    procs.append(p8)

    p9 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_9, args=(data,))
    p9.start()
    procs.append(p9)

    p10 = multiprocessing.Process(target = multiprocessing_smils_to_pdbqt_10, args=(data,))
    p10.start()
    procs.append(p10)

    p11 = multiprocessing.Process(target = multiprocessing_smils_to_pdbqt_11, args=(data,))
    p11.start()
    procs.append(p11)

    p12 = multiprocessing.Process(target = multiprocessing_smils_to_pdbqt_12, args=(data,))
    p12.start()
    procs.append(p12)

    p13 = multiprocessing.Process(target = multiprocessing_smils_to_pdbqt_13, args=(data,))
    p13.start()
    procs.append(p13)

    p14 = multiprocessing.Process(target = multiprocessing_smils_to_pdbqt_14, args=(data,))
    p14.start()
    procs.append(p14)

    p15 = multiprocessing.Process(target = multiprocessing_smils_to_pdbqt_15, args=(data,))
    p15.start()
    procs.append(p15)

    p16 = multiprocessing.Process(target = multiprocessing_smils_to_pdbqt_16, args=(data,))
    p16.start()
    procs.append(p16)

    p17 = multiprocessing.Process(target = multiprocessing_smils_to_pdbqt_17, args=(data,))
    p17.start()
    procs.append(p17)

    p18 = multiprocessing.Process(target = multiprocessing_smils_to_pdbqt_18, args=(data,))
    p18.start()
    procs.append(p18)

    p19 = multiprocessing.Process(target = multiprocessing_smils_to_pdbqt_19, args=(data,))
    p19.start()
    procs.append(p19)

    p20 = multiprocessing.Process(target = multiprocessing_smils_to_pdbqt_20, args=(data,))
    p20.start()
    procs.append(p20)

    p21 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_21, args=(data,))
    p21.start()
    procs.append(p21)

    p22 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_22, args=(data,))
    p22.start()
    procs.append(p22)

    p23 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_23, args=(data,))
    p23.start()
    procs.append(p23)

    p24 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_24, args=(data,))
    p24.start()
    procs.append(p24)

    p25 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_25, args=(data,))
    p25.start()
    procs.append(p25)

    p26 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_26, args=(data,))
    p26.start()
    procs.append(p26)

    p27 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_27, args=(data,))
    p27.start()
    procs.append(p27)

    p28 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_28, args=(data,))
    p28.start()
    procs.append(p28)

    p29 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_29, args=(data,))
    p29.start()
    procs.append(p29)

    p30 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_30, args=(data,))
    p30.start()
    procs.append(p30)

    p31 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_31, args=(data,))
    p31.start()
    procs.append(p31)

    p32 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_32, args=(data,))
    p32.start()
    procs.append(p32)

    p33 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_33, args=(data,))
    p33.start()
    procs.append(p33)

    p34 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_34, args=(data,))
    p34.start()
    procs.append(p34)

    p35 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_35, args=(data,))
    p35.start()
    procs.append(p35)

    p36 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_36, args=(data,))
    p36.start()
    procs.append(p36)

    p37 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_37, args=(data,))
    p37.start()
    procs.append(p37)

    p38 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_38, args=(data,))
    p38.start()
    procs.append(p38)

    p39 = multiprocessing.Process(target = multiprocessing_smiles_to_pdbqt_39, args=(data,))
    p39.start()
    procs.append(p39)

    for p in procs:
        p.join()

    print("pdbqt files translation done")
    #====================================================================================================================================================
    #====================================================================================================================================================
    #====================================================================================================================================================

    ligand_resname ="VG8"
    ligand = structure.select_atoms(f"resname {ligand_resname}")
    pocket_center = (ligand.positions.max(axis=0) + ligand.positions.min(axis=0)) / 2
    pocket_size = ligand.positions.max(axis=0) - ligand.positions.min(axis=0) + 5
    print(pocket_center, pocket_size)

    # Step 2, docking calculation using multiprocessing.Process
    procs=[]

    p0 = multiprocessing.Process(target = multiprocessing_run_smina_0)
    p0.start()
    procs.append(p0)

    p1 = multiprocessing.Process(target = multiprocessing_run_smina_1)
    p1.start()
    procs.append(p1)

    p2 = multiprocessing.Process(target = multiprocessing_run_smina_2)
    p2.start()
    procs.append(p2)

    p3 = multiprocessing.Process(target = multiprocessing_run_smina_3)
    p3.start()
    procs.append(p3)

    p4 = multiprocessing.Process(target = multiprocessing_run_smina_4)
    p4.start()
    procs.append(p4)

    p5 = multiprocessing.Process(target = multiprocessing_run_smina_5)
    p5.start()
    procs.append(p5)

    p6 = multiprocessing.Process(target = multiprocessing_run_smina_6)
    p6.start()
    procs.append(p6)

    p7 = multiprocessing.Process(target = multiprocessing_run_smina_7)
    p7.start()
    procs.append(p7)

    p8 = multiprocessing.Process(target = multiprocessing_run_smina_8)
    p8.start()
    procs.append(p8)

    p9 = multiprocessing.Process(target = multiprocessing_run_smina_9)
    p9.start()
    procs.append(p9)

    p10 = multiprocessing.Process(target = multiprocessing_run_smina_10)
    p10.start()
    procs.append(p10)

    p11 = multiprocessing.Process(target = multiprocessing_run_smina_11)
    p11.start()
    procs.append(p11)

    p12 = multiprocessing.Process(target = multiprocessing_run_smina_12)
    p12.start()
    procs.append(p12)

    p13 = multiprocessing.Process(target = multiprocessing_run_smina_13,)
    p13.start()
    procs.append(p13)

    p14 = multiprocessing.Process(target = multiprocessing_run_smina_14)
    p14.start()
    procs.append(p14)

    p15 = multiprocessing.Process(target = multiprocessing_run_smina_15)
    p15.start()
    procs.append(p15)

    p16 = multiprocessing.Process(target = multiprocessing_run_smina_16)
    p16.start()
    procs.append(p16)

    p17 = multiprocessing.Process(target = multiprocessing_run_smina_17)
    p17.start()
    procs.append(p17)

    p18 = multiprocessing.Process(target = multiprocessing_run_smina_18)
    p18.start()
    procs.append(p18)

    p19 = multiprocessing.Process(target = multiprocessing_run_smina_19)
    p19.start()
    procs.append(p19)

    p20 = multiprocessing.Process(target = multiprocessing_run_smina_20)
    p20.start()
    procs.append(p20)

    p21 = multiprocessing.Process(target = multiprocessing_run_smina_21)
    p21.start()
    procs.append(p21)

    p22 = multiprocessing.Process(target = multiprocessing_run_smina_22)
    p22.start()
    procs.append(p22)

    p23 = multiprocessing.Process(target = multiprocessing_run_smina_23)
    p23.start()
    procs.append(p23)

    p24 = multiprocessing.Process(target = multiprocessing_run_smina_24)
    p24.start()
    procs.append(p24)

    p25 = multiprocessing.Process(target = multiprocessing_run_smina_25)
    p25.start()
    procs.append(p25)

    p26 = multiprocessing.Process(target = multiprocessing_run_smina_26)
    p26.start()
    procs.append(p26)

    p27 = multiprocessing.Process(target = multiprocessing_run_smina_27)
    p27.start()
    procs.append(p27)

    p28 = multiprocessing.Process(target = multiprocessing_run_smina_28)
    p28.start()
    procs.append(p28)

    p29 = multiprocessing.Process(target = multiprocessing_run_smina_29)
    p29.start()
    procs.append(p29)

    p30 = multiprocessing.Process(target = multiprocessing_run_smina_30)
    p30.start()
    procs.append(p30)

    p31 = multiprocessing.Process(target = multiprocessing_run_smina_31)
    p31.start()
    procs.append(p31)

    p32 = multiprocessing.Process(target = multiprocessing_run_smina_32)
    p32.start()
    procs.append(p32)

    p33 = multiprocessing.Process(target = multiprocessing_run_smina_33)
    p33.start()
    procs.append(p33)

    p34 = multiprocessing.Process(target = multiprocessing_run_smina_34)
    p34.start()
    procs.append(p34)

    p35 = multiprocessing.Process(target = multiprocessing_run_smina_35)
    p35.start()
    procs.append(p35)

    p36 = multiprocessing.Process(target = multiprocessing_run_smina_36)
    p36.start()
    procs.append(p36)

    p37 = multiprocessing.Process(target = multiprocessing_run_smina_37)
    p37.start()
    procs.append(p37)

    p38 = multiprocessing.Process(target = multiprocessing_run_smina_38)
    p38.start()
    procs.append(p38)

    p39 = multiprocessing.Process(target = multiprocessing_run_smina_39)
    p39.start()
    procs.append(p39)
    
    for p in procs:
        p.join()

    print("docking done")
    #====================================================================================================================================================
    #====================================================================================================================================================
    #====================================================================================================================================================