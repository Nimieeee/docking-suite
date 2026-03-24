# Docking Suite

Comprehensive utilities and an example pipeline for running molecular docking workflows locally.

## Table of contents

- [Features](#features)
- [Repository structure](#repository-structure)
- [Quick start](#quick-start)
- [Installation & dependencies](#installation--dependencies)
- [Usage examples](#usage-examples)
- [Example API (docking_utils)](#example-api-docking_utils)
- [Jupyter notebook walkthrough](#jupyter-notebook-walkthrough)
- [Testing & development](#testing--development)
- [Contributing](#contributing)
- [License & contact](#license--contact)

## Features

- Helpers to prepare receptors and ligands for docking
- Parsers to extract and summarize docking results
- Example Jupyter notebook demonstrating an end-to-end workflow

## Repository structure

- `docking_pipeline.ipynb` — end-to-end example notebook
- `docking_utils.py` — core helper functions used by the notebook and scripts
- `README.md` — this file

Additions you might expect: `requirements.txt`, `LICENSE`, and CI configs at the project root.

## Quick start

1. Create and activate a Python environment (recommended using `conda`):

   conda create -n docking python=3.9 -y
   conda activate docking

2. Install project dependencies (if a `requirements.txt` is present at repo root):

   pip install -r ../requirements.txt

3. Run the notebook locally:

   jupyter lab docking_pipeline.ipynb

## Installation & dependencies

Core Python packages commonly used in docking workflows (install the ones you need):

- numpy, pandas — data handling
- rdkit — chemistry utilities (if used)
- openbabel / pybel — format conversion (optional)
- biopython — PDB handling (optional)

If you use a docking engine such as AutoDock Vina, DOCK, or Glide, install and configure that engine separately and ensure its binaries are on your `PATH`.

## Usage examples

1) Prepare inputs and run docking (high-level example):

   - Prepare receptor and ligands using the helper functions in `docking_utils.py`.
   - Call your docking engine (example shown for AutoDock Vina):

   vina --receptor receptor.pdbqt --ligand ligand.pdbqt --out out.pdbqt --log out.log

2) Parse results and summarize:

   - Use the parsing utilities to extract poses and scores into a CSV for downstream filtering and analysis.

## Example API (docking_utils)

Below is a suggested usage pattern; check `docking_utils.py` for exact function names and signatures.

```py
from docking_utils import prepare_receptor, prepare_ligand, run_docking, parse_results

# prepare files
prepare_receptor('receptor.pdb', out='receptor.pdbqt')
prepare_ligand('ligand.sdf', out='ligand.pdbqt')

# run docking (this function may wrap a call to a docking engine)
run_docking(receptor='receptor.pdbqt', ligand='ligand.pdbqt', out='results.pdbqt')

# parse results
df = parse_results('results.pdbqt')
print(df.head())
```

If the exact function names differ, refer to `docking_utils.py` for the canonical API.

## Jupyter notebook walkthrough

- Open `docking_pipeline.ipynb` to see a step-by-step example covering:
  - Data preparation (ligands, receptors)
  - Running docking jobs (locally or via a cluster)
  - Parsing and visualizing results

## Testing & development

- Add unit tests for core functions in `docking_utils.py` (example: tests/test_utils.py)
- Use `pytest` to run tests:

  pytest -q

## Contributing

- Fork the repo, create a feature branch, and open a PR with a clear description and tests/examples.
- Suggested improvements: expand `docking_utils` docstrings, add a `requirements.txt`, include CI, and add small example datasets.

## License & contact

Add a `LICENSE` file at the repository root. If none exists, contact the maintainer before reuse.

Maintainer / Contact: Nimieeee — https://github.com/Nimieeee/docking-suite

---

If you'd like, I can:

- add a `requirements.txt` with commonly used packages
- generate a `LICENSE` (MIT / Apache-2.0) and commit it
- add simple unit tests for `docking_utils.py`

Tell me which of the above you'd like next.
