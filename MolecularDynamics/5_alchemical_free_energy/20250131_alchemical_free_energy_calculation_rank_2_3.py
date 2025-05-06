import openmm as mm
from openmm import *
import openmm.app as app
from openmm.app import *
import openmm.unit as unit
from openmm.unit import *
from openmmtools import alchemy, integrators, states, multistate, mcmc
from openmmtools.mcmc import LangevinSplittingDynamicsMove
from openmmtools.states import *
from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalState, AlchemicalRegion
from openmmtools.integrators import LangevinIntegrator
from openmmtools.multistate import ReplicaExchangeSampler, MultiStateReporter, MultiStateSamplerAnalyzer
from openmmtools.testsystems import LennardJonesFluid
import pdbfixer
import numpy as np
import warnings
from tqdm import tqdm

from pymbar import MBAR, timeseries

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from IPython.display import display
from openff.toolkit.topology import Molecule, Topology
import mdtraj as md
from openmmforcefields.generators import GAFFTemplateGenerator

def prepare_protein(
    pdb_file, ignore_missing_residues=True, ignore_terminal_missing_residues=True, ph=7.0
):
    """
    Use pdbfixer to prepare the protein from a PDB file. Hetero atoms such as ligands are
    removed and non-standard residues replaced. Missing atoms to existing residues are added.
    Missing residues are ignored by default, but can be included.

    Parameters
    ----------
    pdb_file: pathlib.Path or str
        PDB file containing the system to simulate.
    ignore_missing_residues: bool, optional
        If missing residues should be ignored or built.
    ignore_terminal_missing_residues: bool, optional
        If missing residues at the beginning and the end of a chain should be ignored or built.
    ph: float, optional
        pH value used to determine protonation state of residues

    Returns
    -------
    fixer: pdbfixer.pdbfixer.PDBFixer
        Prepared protein system.
    """
    fixer = pdbfixer.PDBFixer(str(pdb_file))
    fixer.removeHeterogens()  # co-crystallized ligands are unknown to PDBFixer 이 과정에서 ligand 삭제됨(hetro atom 이니까)
    fixer.findMissingResidues()  # identify missing residues, needed for identification of missing atoms

    # if missing terminal residues shall be ignored, remove them from the dictionary
    if ignore_terminal_missing_residues:
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

    # if all missing residues shall be ignored ignored, clear the dictionary
    if ignore_missing_residues:
        fixer.missingResidues = {}

    fixer.findNonstandardResidues()  # find non-standard residue
    fixer.replaceNonstandardResidues()  # replace non-standard residues with standard one
    fixer.findMissingAtoms()  # find missing heavy atoms
    fixer.addMissingAtoms()  # add missing atoms and residues
    fixer.addMissingHydrogens(ph)  # add missing hydrogens
    return fixer

def prepare_ligand(pdb_file, resname, smiles, depict=True):
    """
    Prepare a ligand from a PDB file via adding hydrogens and assigning bond orders. A depiction
    of the ligand before and after preparation is rendered in 2D to allow an inspection of the
    results. Huge thanks to @j-wags for the suggestion.

    Parameters
    ----------
    pdb_file: pathlib.PosixPath
       PDB file containing the ligand of interest.
    resname: str
        Three character residue name of the ligand.
    smiles : str
        SMILES string of the ligand informing about correct protonation and bond orders.
    depict: bool, optional
        show a 2D representation of the ligand

    Returns
    -------
    prepared_ligand: rdkit.Chem.rdchem.Mol
        Prepared ligand.
    """
    # split molecule
    rdkit_mol = Chem.MolFromPDBFile(str(pdb_file), proximityBonding=True, sanitize=False)
    print(rdkit_mol)
    rdkit_mol_split = Chem.rdmolops.SplitMolByPDBResidues(rdkit_mol)

    # extract the ligand and remove any already present hydrogens
    ligand = rdkit_mol_split[resname]
    ligand = Chem.RemoveHs(ligand)

    # assign bond orders from template
    reference_mol = Chem.MolFromSmiles(smiles)
    prepared_ligand = AllChem.AssignBondOrdersFromTemplate(reference_mol, ligand)
    prepared_ligand.AddConformer(ligand.GetConformer(0))

    # protonate ligand
    prepared_ligand = Chem.rdmolops.AddHs(prepared_ligand, addCoords=True)
    prepared_ligand = Chem.MolFromMolBlock(Chem.MolToMolBlock(prepared_ligand))

    # 2D depiction
    if depict:
        ligand_2d = copy.deepcopy(ligand)
        prepared_ligand_2d = copy.deepcopy(prepared_ligand)
        AllChem.Compute2DCoords(ligand_2d)
        AllChem.Compute2DCoords(prepared_ligand_2d)
        display(
            Draw.MolsToGridImage(
                [ligand_2d, prepared_ligand_2d], molsPerRow=2, legends=["original", "prepared"]
            )
        )

    # return ligand
    return prepared_ligand

def rdkit_to_openmm(rdkit_mol, name="LIG"):
    """
    Convert an RDKit molecule to an OpenMM molecule.
    Inspired by @hannahbrucemcdonald and @glass-w.

    Parameters
    ----------
    rdkit_mol: rdkit.Chem.rdchem.Mol
        RDKit molecule to convert.
    name: str
        Molecule name.

    Returns
    -------
    omm_molecule: openmm.app.Modeller
        OpenMM modeller object holding the molecule of interest.
    """
    # convert RDKit to OpenFF
    off_mol = Molecule.from_rdkit(rdkit_mol)
    #print(off_mol)

    # add name for molecule
    off_mol.name = name
    #print(off_mol.name)

    # add names for atoms
    element_counter_dict = {}
    for off_atom, rdkit_atom in zip(off_mol.atoms, rdkit_mol.GetAtoms()):
        element = rdkit_atom.GetSymbol()
        if element in element_counter_dict.keys():
            element_counter_dict[element] += 1
        else:
            element_counter_dict[element] = 1
        off_atom.name = element + str(element_counter_dict[element])
    #print(element_counter_dict)

    # convert from OpenFF to OpenMM
    off_mol_topology = off_mol.to_topology()
    mol_topology = off_mol_topology.to_openmm()
    #print(mol_topology)
    mol_positions = off_mol.conformers[0]
    #print(mol_positions)

    # convert units from Ångström to nanometers
    # since OpenMM works in nm
    mol_positions = mol_positions.to("nanometers")
    #print(mol_positions)

    # combine topology and positions in modeller object
    omm_mol = app.Modeller(mol_topology, mol_positions)

    return omm_mol

def merge_protein_and_ligand(protein, ligand):
    """
    Merge two OpenMM objects.

    Parameters
    ----------
    protein: pdbfixer.pdbfixer.PDBFixer
        Protein to merge.
    ligand: openmm.app.Modeller
        Ligand to merge.

    Returns
    -------
    complex_topology: openmm.app.topology.Topology
        The merged topology.
    complex_positions: openmm.unit.quantity.Quantity
        The merged positions.
    """
    # combine topologies
    md_protein_topology = md.Topology.from_openmm(protein.topology)  # using mdtraj for protein top
    #print(md_protein_topology)
    md_ligand_topology = md.Topology.from_openmm(ligand.topology)  # using mdtraj for ligand top
    #print(md_ligand_topology)
    md_complex_topology = md_protein_topology.join(md_ligand_topology)  # add them together
    #print(md_complex_topology)
    complex_topology = md_complex_topology.to_openmm()
    #print(complex_topology)

    # combine positions
    total_atoms = len(protein.positions) + len(ligand.positions)

    # create an array for storing all atom positions as tupels containing a value and a unit
    # called OpenMM Quantities
    complex_positions = unit.Quantity(np.zeros([total_atoms, 3]), unit=unit.nanometers)
    complex_positions[: len(protein.positions)] = protein.positions  # add protein positions
    complex_positions[len(protein.positions) :] = ligand.positions  # add ligand positions

    return complex_topology, complex_positions

def generate_forcefield(
    rdkit_mol=None, protein_ff="amber14-all.xml", solvent_ff="amber14/tip3pfb.xml"
):
    """
    Generate an OpenMM Forcefield object and register a small molecule.

    Parameters
    ----------
    rdkit_mol: rdkit.Chem.rdchem.Mol
        Small molecule to register in the force field.
    protein_ff: string
        Name of the force field.
    solvent_ff: string
        Name of the solvent force field.

    Returns
    -------
    forcefield: openmm.app.Forcefield
        Forcefield with registered small molecule.
    """
    forcefield = app.ForceField(protein_ff, solvent_ff)

    if rdkit_mol is not None:
        gaff = GAFFTemplateGenerator(
            molecules=Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
        )
        forcefield.registerTemplateGenerator(gaff.generator)

    return forcefield

if __name__ == "__main__":
    warnings.filterwarnings("ignore")

    pdb_path = "./rank2.pdb"
    ligand_name = "UNL"
    smiles = "COc1ccccc1N1CCN(C(=O)Nc2ccccc2F)CC1"

    # Load the protein and ligand system
    protein = prepare_protein(pdb_path, ignore_missing_residues=False)
    # prepare ligand
    rdkit_ligand = prepare_ligand(pdb_path, ligand_name, smiles)
    # rdkit object to openmm object
    omm_ligand = rdkit_to_openmm(rdkit_ligand, ligand_name)
    # merge protein and ligand
    complex_topology, complex_positions = merge_protein_and_ligand(protein, omm_ligand)
    print("Complex topology has", complex_topology.getNumAtoms(), "atoms.")

    # Step 2: Define the force field and create the system
    forcefield = generate_forcefield(rdkit_ligand)
    modeller = app.Modeller(complex_topology, complex_positions)
    # Add solvent (optional if already solvated)
    modeller.addSolvent(forcefield, padding=1.5*unit.nanometers, ionicStrength=0.15*unit.molar)
    modeller.topology.getPeriodicBoxVectors()

    # save the modeller to pdb file
    with open("./rank2_processed_structure.pdb", "w") as output_pdb:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, output_pdb)

    # Create the OpenMM system
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds)

    # Simulation settings
    pressure = 1*unit.atmospheres
    temperature = 309.65*unit.kelvin
    collision_rate = 1.0/unit.picoseconds
    timestep = 2.0*unit.femtoseconds

    sigma = 3.4*unit.angstrom
    epsilon = 0.238 * unit.kilocalories_per_mole

    barostat = MonteCarloBarostat(pressure, temperature)
    system.addForce(barostat)

    # Retrieve the NonbondedForce
    forces = { force.__class__.__name__ : force for force in system.getForces() }
    nbforce = forces['NonbondedForce']

    # Define the alchemical region (ligand)
    ligand_atoms = [atom.index for atom in modeller.topology.atoms() if atom.residue.name == "UNK"]
    # ligand atoms이 잘 정의되었는지 확인하기 위함
    print("ligand_atoms = ", ligand_atoms)
    protein_atoms = set(range(system.getNumParticles())) - set(ligand_atoms)
    alchemical_region = AlchemicalRegion(alchemical_atoms=ligand_atoms)
    alchemical_factory = AbsoluteAlchemicalFactory(consistent_exceptions=False)
    alchemical_system = alchemical_factory.create_alchemical_system(system, alchemical_region)

    # Define the energy function for the CustomNonbondedForce
    # when lambda is 1.0 it is a normal LJ potential, when lambda is 0.0 the interaction vanishes
    energy_function = 'lambda*4*epsilon*x*(x-1.0); x = (sigma/reff_sterics)^6;'
    energy_function += 'reff_sterics = sigma*(0.5*(1.0-lambda) + (r/sigma)^6)^(1/6);'
    energy_function += 'sigma = 0.5*(sigma1+sigma2); epsilon = sqrt(epsilon1*epsilon2);'
    custom_force = CustomNonbondedForce(energy_function)

    # Add lambda as a parameter we can change during the simulation
    custom_force.addGlobalParameter('lambda', 1.0)

    # set the values of sigma and epsilon by copying them from the existing NonBondedForce
    custom_force.addPerParticleParameter('sigma')
    custom_force.addPerParticleParameter('epsilon')
    for index in range(alchemical_system.getNumParticles()):
        [charge, sigma, epsilon] = nbforce.getParticleParameters(index)
        custom_force.addParticle([sigma, epsilon])
        if index in ligand_atoms:
            # remove the alchemical particle from the existing NonBondedForce
            nbforce.setParticleParameters(index, charge*0, sigma, epsilon*0)

    # Set the custom force to occur between just the alchemical particle and the other particles
    custom_force.addInteractionGroup(ligand_atoms, protein_atoms)
    alchemical_system.addForce(custom_force)


    # Create a Lennard Jones test fluid
    sigma = 3.4*unit.angstrom
    epsilon = 0.238 * unit.kilocalories_per_mole

    # Create an integrator
    integrator = LangevinIntegrator(temperature, collision_rate, timestep)

    # set GPU
    platform = mm.Platform.getPlatformByName("CUDA")
    properties = {"DeviceIndex": "0", "Precision": "single"}

    # Create a simulation
    simulation = app.Simulation(modeller.topology, alchemical_system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # Minimize energy
    print('Minimizing energy...')
    #LocalEnergyMinimizer.minimize(context)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(temperature)

    # Collect data

    # number of steps per sample
    nsteps = 5000
    # number of samples to collect per alchemical state
    niterations = 100

    import numpy as np
    lambdas = np.linspace(1.0, 0.0, 10) # alchemical lambda schedule
    nstates = len(lambdas)
    u_kln = np.zeros([nstates,nstates,niterations], np.float64)
    kT = AVOGADRO_CONSTANT_NA * BOLTZMANN_CONSTANT_kB * integrator.getTemperature()
    for k in tqdm(range(nstates), desc="Processing..."):
        for iteration in range(niterations):
            print('state %5d iteration %5d / %5d' % (k, iteration, niterations))
            # Set alchemical state
            simulation.context.setParameter('lambda', lambdas[k])
            # Run some dynamics
            simulation.step(nsteps)
            # Compute energies at all alchemical states
            for l in range(nstates):
                simulation.context.setParameter('lambda', lambdas[l])
                u_kln[k,l,iteration] = simulation.context.getState(getEnergy=True).getPotentialEnergy() / kT

    ## Subsample data to extract uncorrelated equilibrium timeseries
    N_k = np.zeros([nstates], np.int32) # number of uncorrelated samples
    for k in range(nstates):
        [nequil, g, Neff_max] = timeseries.detectEquilibration(u_kln[k,k,:])
        indices = timeseries.subsampleCorrelatedData(u_kln[k,k,:], g=g)
        N_k[k] = len(indices)
        u_kln[k,:,0:N_k[k]] = u_kln[k,:,indices].T

    # Compute free energy differences
    mbar = MBAR(u_kln, N_k)

    # dont compute uncertainties here, if you do it may fail with an error for
    # pymbar versions > 3.0.3. See this issue: https://github.com/choderalab/pymbar/issues/419
    [DeltaF_ij] = mbar.getFreeEnergyDifferences(compute_uncertainty=False)
    results = DeltaF_ij[nstates-1][0] * kT

    print("Free energy change to insert a particle = ", results)

    with open("./20250208_free_energy_results_rank2_3.txt", "w") as f:
        f.write("Free energy change to insert a particle = " + str(DeltaF_ij[nstates-1][0]) + "\n")
        f.write(str(DeltaF_ij) + "\n")
        f.write("steps = " +str(nsteps) + "niterations = " + str(niterations) + "\n")


    # 결과값의 unit : J/mol