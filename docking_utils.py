import os
import shutil
import tempfile
import time
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from tqdm import tqdm
import mdtraj
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina

# Conditionally imported libraries with error handling
try:
    from openbabel import pybel
except ImportError:
    print("Warning: OpenBabel and Pybel are required for full functionality.")

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
except ImportError:
    print("Warning: PDBFixer and OpenMM are required to fix the protein structure.")

class MoleculePreparationConfig:
    """Configuration for molecule preparation."""
    def __init__(self, add_hydrogens=True, make_3d=True, charge_model="gasteiger"):
        self.add_hydrogens = add_hydrogens
        self.make_3d = make_3d
        self.charge_model = charge_model

class MoleculeReader:
    """Reads molecule data from files."""
    @staticmethod
    def read_molecule(file: Union[str, os.PathLike], file_format: str = "pdb") -> List[Any]:
        """Reads molecules from a file."""
        return list(pybel.readfile(format=file_format, filename=str(file)))

class MoleculePreparator:
    """Prepares molecules for docking."""
    def __init__(self, config: MoleculePreparationConfig = MoleculePreparationConfig()):
        self.config = config

    def prepare_molecule(self, molecule: Any) -> Any:
        """Prepares a single molecule for docking."""
        if self.config.add_hydrogens:
            molecule.addh()
        if self.config.make_3d and not molecule.OBMol.HasNonZeroCoords():
            molecule.make3D(forcefield="mmff94s", steps=10000)
        molecule.calccharges(model=self.config.charge_model)
        return molecule

    def save_molecule(self, molecule: Any, outpath: Union[str, os.PathLike], 
                      out_format: str = "pdbqt", overwrite: bool = False, rigid: bool = False) -> None:
        """Saves a molecule to a file."""
        if rigid:
            opt = {"r": None, "c": None, "h": None}
        else:
            opt = {"b": None, "p": None, "h": None}
            
        with pybel.Outputfile(format=out_format, filename=str(outpath), overwrite=overwrite, opt=opt) as out:
            out.write(molecule)

class Preprocessor:
    """Handles preparation of receptor and ligand."""
    def __init__(self, config: MoleculePreparationConfig = MoleculePreparationConfig(), fix_protein: bool = False):
        self.preparator = MoleculePreparator(config)
        self.fix_protein = fix_protein

    def _fix_receptor_structure(self, receptor_path: str, output_path: str, keep_water: bool = False) -> str:
        """Fixes protein structure using PDBFixer."""
        print("Fixing protein structure with PDBFixer...")
        fixer = PDBFixer(filename=receptor_path)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=keep_water)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)
        PDBFile.writeFile(fixer.topology, fixer.positions, open(output_path, 'w'))
        return output_path
    
    def prepare_receptor(self, receptor_path: str, output_path: str, keep_water: bool = False) -> None:
        fixed_receptor_path = receptor_path
        if self.fix_protein:
            fixed_receptor_path = os.path.splitext(receptor_path)[0] + "_fixed.pdb"
            fixed_receptor_path = self._fix_receptor_structure(receptor_path, fixed_receptor_path, keep_water=keep_water)

        molecules = MoleculeReader.read_molecule(fixed_receptor_path)
        prepared_receptor = self.preparator.prepare_molecule(molecules[0])
        self.preparator.save_molecule(prepared_receptor, output_path, rigid=True)
        print(f"Prepared receptor saved to: {output_path}")

    def prepare_ligand(self, ligand_path: str, output_path: str, in_format: str = "pdb") -> None:
        molecules = MoleculeReader.read_molecule(ligand_path, in_format)
        prepared_ligand = self.preparator.prepare_molecule(molecules[0])
        self.preparator.save_molecule(prepared_ligand, output_path, out_format="sdf")
        print(f"Prepared ligand saved to: {output_path}")

class VinaDocking:
    """Performs docking using AutoDock Vina."""
    def __init__(self, receptor_path: str, center: List[float], size: List[float], num_poses: int = 5, exhaustiveness: int = 8):
        self.v = Vina(sf_name='vina', cpu=0, verbosity=0)
        self.v.set_receptor(receptor_path)
        self.center = center
        self.size = size
        self.num_poses = num_poses
        self.exhaustiveness = exhaustiveness

    def dock(self, ligand_path: str, out_path: str = None) -> pd.DataFrame:
        self.v.set_ligand_from_file(ligand_path)
        self.v.compute_vina_maps(center=self.center, box_size=self.size)
        self.v.dock(exhaustiveness=self.exhaustiveness, n_poses=self.num_poses)
        if out_path:
            self.v.write_poses(out_path, n_poses=self.num_poses, overwrite=True)
        energies = self.v.energies(n_poses=self.num_poses)
        return pd.DataFrame(energies, columns=["affinity", "inter", "intra", "torsions", "intra_best_post"])

def get_box_from_ligand(ligand_file: str, padding: float = 5.0) -> Tuple[List[float], List[float]]:
    """Calculates docking box center and size from a ligand file."""
    traj = mdtraj.load(ligand_file)
    xyz = traj.xyz[0] * 10 # to Angstroms
    center = ((xyz.max(axis=0) + xyz.min(axis=0)) / 2).tolist()
    size = (xyz.max(axis=0) - xyz.min(axis=0) + padding).tolist()
    return center, size


class MMGBSACalculator:
    """Calculates binding energy using MM-GBSA (Molecular Mechanics/Generalized Born Surface Area)."""
    
    def __init__(self, implicit_solvent='OBC2'):
        self.implicit_solvent = implicit_solvent
        self._check_dependencies()
    
    def _check_dependencies(self):
        """Verify required dependencies are available."""
        try:
            import openmm
            import openmm.app
            import openmm.unit as unit
            from openforcefield.topology import Molecule as OFFMolecule
            from openforcefield.topology import Topology as OFFTopology
            from openforcefield.typing.engines.smirnoff import ForceField
        except ImportError as e:
            raise ImportError(f"MM-GBSA requires OpenMM and OpenForceField: {e}")
    
    def _load_protein(self, protein_pdb: str):
        """Load protein structure using OpenMM."""
        import openmm.app as app
        pdb = app.PDBFile(protein_pdb)
        return pdb.topology, pdb.positions
    
    def _load_ligand(self, ligand_sdf: str):
        """Load ligand from SDF using RDKit and convert to OpenFF."""
        from rdkit import Chem
        from openforcefield.topology import Molecule as OFFMolecule
        from openforcefield.topology import Topology as OFFTopology
        
        mol = Chem.MolFromMolFile(ligand_sdf, sanitize=False)
        mol.UpdatePropertyCache()
        Chem.rdmolops.AssignStereochemistryFrom3D(mol)
        
        off_mol = OFFMolecule.from_rdkit(mol)
        return off_mol
    
    def _create_system(self, topology, *molecules):
        """Create OpenMM System with implicit solvent."""
        import openmm
        import openmm.app as app
        import openmm.unit as unit
        from openforcefield.topology import Topology as OFFTopology
        
        if molecules:
            off_topology = OFFTopology.from_openmm(topology)
            for mol in molecules:
                off_topology.add_molecule(mol)
            
            ff = ForceField('gaff-2.11.xml')
            system = ff.create_system(off_topology, nonbondedMethod=app.NoCutoff)
        else:
            off_topology = OFFTopology.from_openmm(topology)
            ff = ForceField('gaff-2.11.xml')
            system = ff.create_system(off_topology, nonbondedMethod=app.NoCutoff)
        
        if self.implicit_solvent:
            system.addForce(openmm.GBSACalculationForce(self.implicit_solvent))
        
        return system
    
    def _calculate_energy(self, topology, positions, system):
        """Calculate potential energy of a system."""
        import openmm
        import openmm.app as app
        import openmm.unit as unit
        
        integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
        platform = openmm.Platform.getPlatformByReference()
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions)
        
        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        
        del context, integrator
        return energy
    
    def calculate_score(self, protein_pdb: str, ligand_sdf: str) -> float:
        """
        Calculate MM-GBSA binding energy.
        
        Parameters:
        -----------
        protein_pdb : str
            Path to protein PDB file
        ligand_sdf : str
            Path to ligand SDF file
            
        Returns:
        --------
        float
            Binding energy in kcal/mol (ΔG = E_complex - E_protein - E_ligand)
        """
        import tempfile
        import os
        
        protein_topology, protein_positions = self._load_protein(protein_pdb)
        ligand_off_mol = self._load_ligand(ligand_sdf)
        
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            tmp_protein_pdb = tmp.name
        
        try:
            import openmm.app as app
            app.PDBFile.writeFile(protein_topology, protein_positions, open(tmp_protein_pdb, 'w'))
            
            complex_topology, complex_positions = self._load_protein(tmp_protein_pdb)
            
            complex_system = self._create_system(complex_topology, ligand_off_mol)
            protein_system = self._create_system(protein_topology)
            ligand_topology = ligand_off_mol.to_topology()
            ligand_system = self._create_system(ligand_topology, ligand_off_mol)
            
            e_complex = self._calculate_energy(complex_topology, complex_positions, complex_system)
            e_protein = self._calculate_energy(protein_topology, protein_positions, protein_system)
            e_ligand = self._calculate_energy(ligand_topology, ligand_off_mol.conformers[0], ligand_system)
            
            binding_energy = e_complex - (e_protein + e_ligand)
            return binding_energy
            
        finally:
            if os.path.exists(tmp_protein_pdb):
                os.remove(tmp_protein_pdb)
