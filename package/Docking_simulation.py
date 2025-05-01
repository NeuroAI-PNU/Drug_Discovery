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

class Docking():
    """
    Docking class for docking simulations.
    """

    def __init__(self, data):
        self.data = data

    def smiles_to_pdbqt(self, smiles, pdbqt_path, pH=7.4):
        """
        Convert a SMILES string to a PDBQT file needed by docking programs of the AutoDock family.
        
        Parameters
        ----------
        pdb_id : str
            PDB ID of the protein structure.
        smiles : str
            SMILES string.
        pdbqt_path : str or pathlib.path
            Path to output PDBQT file.
        pH : float
            Protonation at given pH.
        """
        molecule = pybel.readstring("smi", smiles)
        # add hydrogens at giben pH
        molecule.OBMol.CorrectForPH(pH)
        molecule.addh()
        # generate 3D coordinates
        molecule.make3D(forecefield='mmff94s', steps=100000)
        # add partial charges to each atom
        for atom in molecule.atoms:
            atom.OBAtom.GetPartialCharge()
        molecule.write("pdbqt", str(pdbqt_path), overwrite=True)
        return
    
    def preparation_for_docking(self, pdb_id, ligand_name, smiles, pdbqt_path):
        """
        generate pdbqt files for ligand and protein.
        and make pocket information.

        Parameters
        ----------
        pdb_id : str
            PDB ID of the protein structure.
        smiles : str
            SMILES string of ligands.
        ligand_name : str
            Residue name of the ligand in the protein structure.
        pdbqt_path : str or pathlib.Path
            Path to the output PDBQT file.

        Returns
        -------
        structure : opencadd.structure.core.Structure
            Protein structure.
        pocket_center : iterable of float
            Coordinates defining the center of the binding site.
        pocket_size : iterable of float
            Lengths of edges defining the binding site.
        """
        warnings.filterwarnings("ignore")
        ob_log_handler = pybel.ob.OBMessageHandler()
        pybel.ob.obErrorLog.SetOutputLevel(0)
        print("warining setting is done")

        structure = Structure.from_pdbid(pdb_id)

        #pdb_file = "./data/2jkq.pdb"
        #structure = Structure(pdb_file)
        # element information maybe missing, but important for subsequent pdbqt conversion
        if not hasattr(structure.atoms, "elements"):
            structure.add_TopologyAttr("elements", structure.atoms.types)
        print("structure load done")

        smiles_to_pdbqt(smiles, pdbqt_path)
    
        ligand = structure.select_atoms(f"resname {ligand_name}")
        pocket_center = (ligand.positions.max(axis=0) + ligand.positions.min(axis=0)) / 2
        pocket_size = ligand.positions.max(axis=0) - ligand.positions.min(axis=0) + 5
        print(pocket_center, pocket_size)
        return  structure, pocket_center, pocket_size

    
    def run_smina(
            ligand_path, protein_path, out_path, pocket_center, pocket_size, num_poses=10, exhaustiveness=8):
        """
        Perform dociking with Smina.
        
        Parameters
        ----------
        ligand_path : str or pathlib.Path
            Path to ligand PDBQT file that should be docked.
        protein_path : str or pathlib.Path
            Path to protein PDBQT file.
        out_path : str or pathlib.Path
            Path to whick docking poses should be saved, SDF or PDB format.
        pocket_center : iterable of float or int
            Coordinates defining the center of the binding site.
        pocket_size : iterable of float or int
            Lengths of edges defining the binding site.
        num_poses : int
            Maximum number of poses to generate.
        exhaustiveness : int
            Accuracy of docking calculations.
        
        Returns
        -------
        output_text : str
            The output of the Smina calculation.
        """
        output_text = subprocess.check_output(
            [
                # In Max
                #'./smina.osx.12",

                # In Lunux
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
            universal_newlines=True, # needed to capture output text
        )
        return output_text
    
    def extract_docking_results(self, n, ligand_path, protein_path, out_path, pocket_center, pocket_size, num_poses=10, exhaustiveness=8):
        """
        Extract docking results from SDF file.
        
        Parameters
        ----------
        n : int
            Number of ligands that were docked.
        ligand_path : str or pathlib.Path
            Path to ligand PDBQT file that should be docked.
        protein_path : str or pathlib.Path
            Path to protein PDBQT file that should be docked to.
        out_path : str or pathlib.Path
            Path to the SDF file containing docking results.
        pocket_center : iterable of float or int
            Coordinates defining the center of the binding site.
        pocket_size : iterable of float or int
            Lengths of edges defining the binding site.
        num_poses : int
            Maximum number of poses to generate.
        exhaustiveness : int
            Accuracy of docking calculations

        Returns
        -------
        docking_results : list of dict
            List of dictionaries containing docking results.
        """
        for i in range(n):
            docking = open(f"docking_score_{i}.txt", mode="w")
            results = run_smina(
                ligand_path, protein_path, out_path, pocket_center, pocket_size
            )
            docking.write(results)
            docking.close()
