# Docking System Overview & Architecture

This document provides a comprehensive summary of the molecular docking system we have built, the underlying technologies, and server performance considerations.

## 1. System Architecture

The docking system is organized as a standalone suite in the `docking_suite/` directory. It follows a modular design:

- **Docking Engine**: [AutoDock Vina 1.2.x](https://vina.scripps.edu/) (Python API).
- **Structure Preparation**:
    - **Protein**: [PDBFixer](https://github.com/openmm/pdbfixer) (OpenMM) for fixing residues and adding atoms.
    - **Ligand**: [Open Babel (Pybel)](https://openbabel.org/) and [Meeko](https://github.com/forlilab/Meeko) for PDBQT conversion.
- **Analysis**: [ProLif](https://prolif.readthedocs.io/) for protein-ligand interaction fingerprints.
- **Workflow Strategy**: Bayesian Optimization is used to accelerate screening by focusing on high-probability candidates.

---

## 2. Server Performance & Concurrency

### Hardware Specs: 2 vCPUs, 8GB RAM

| Factor | Status | Details |
| :--- | :--- | :--- |
| **Concurrency** | **1-2 Users** | Best for 1 user (uses both vCPUs). 2 users can run simultaneously at half speed. |
| **RAM Usage** | **Very Low** | Each docking run uses <500MB. 8GB can easily handle multiple instances. |
| **CPU Load** | **High** | Vina is CPU-bound. 2 vCPUs will be at 100% capacity during a run. |
| **GPU Need** | **None** | The official Vina API runs exclusively on CPU using multithreading. |

### Concurrency Recommendation
- **Ideal**: 1 active docking process at a time.
- **Maximum**: 2 concurrent processes (if using `--cpu 1` each).
- **Bottleneck**: CPU cores are the limiting factor, not RAM.

---

## 3. Advanced Workflows (MM-GBSA)

### Theoretical Background

**MM-GBSA** (Molecular Mechanics/Generalized Born Surface Area) is a post-docking rescoring method that calculates binding free energy using:

$$\Delta G_{binding} = E_{complex} - (E_{protein} + E_{ligand}) + G_{solv}$$

Where:
- $E_{complex}$, $E_{protein}$, $E_{ligand}$ are the molecular mechanics potential energies
- $G_{solv}$ is the solvation free energy estimated using implicit solvent (Generalized Born model)

The method uses:
- **GAFF (Generalized Amber Force Field)** for ligand parameterization via OpenForceField
- **OBC2 (Onufriev-Bashford-Colvo)** implicit solvent model for solvation energy
- Single-point energy evaluation on docked poses (no MD simulation)

### Implementation Details

The `MMGBSACalculator` class in `docking_utils.py`:
1. Loads protein (PDB) and ligand (SDF) structures
2. Parameterizes the ligand using OpenForceField with GAFF
3. Creates an OpenMM System with implicit solvent (OBC2)
4. Calculates potential energy for Complex, Protein, and Ligand separately
5. Returns: $\Delta G = E_{complex} - (E_{protein} + E_{ligand})$

### Performance Considerations

- **Feasibility**: High intensity on 2 vCPUs. Single-point evaluation is much faster than MD-based methods.
- **Recommendation**: Use "single-point" MM-GBSA (on static poses) to stay within the 8GB RAM and 2 vCPU limits.
- **Water Molecules**: By default, protein preparation removes waters for docking efficiency. You can now toggle this using the `keep_water=True` parameter in the `prepare_receptor` method if structural waters are critical.
- **Use Case**: Re-scoring docked poses to improve binding affinity predictions beyond Vina scores.

---

## 4. Component Breakdown

1. **`docking_utils.py`**:
    - `Preprocessor`: Handles structural fixing and PDBQT conversion.
    - `MoleculePreparator`: Adds hydrogens and generates 3D coordinates.
    - `VinaDocking`: Encapsulates the Vina engine and docking logic.
2. **`docking_pipeline.ipynb`**:
    - A plug-and-play Jupyter Notebook that links all components into a user-friendly workflow.
3. **`environment.yml`**:
    - Defines the Conda environment `ml4dd2025` required to run the suite.

---

## 5. Summary of Discussion
We have moved from exploring general notebooks to isolating a specific, optimized pipeline for docking. We identified **ML4DD chapter 9** as the primary technical reference and extracted its best practices into this standalone suite. The setup is designed to be free, open-source, and capable of running on standard server hardware without requiring specialized GPUs.
