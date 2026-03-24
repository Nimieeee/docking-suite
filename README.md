# Docking Suite

A small docking utilities and pipeline collection used for molecular docking workflows.

## Contents

- `docking_pipeline.ipynb` — Jupyter notebook demonstrating the pipeline.
- `docking_utils.py` — Helper functions for preparing and analyzing docking runs.

## Quick start

1. Create and activate a Python environment (recommended using `conda` or `venv`).

   conda create -n docking python=3.9 -y
   conda activate docking

2. Install required packages (adjust to your environment):

   pip install -r ../requirements.txt

3. Open the notebook:

   jupyter lab docking_pipeline.ipynb

## Usage

- Use `docking_utils.py` to prepare input files, run docking with your preferred engine, and parse results.
- The notebook contains an end-to-end example showing data preparation, running, and result analysis.

## Contributing

If you'd like to improve this suite, please open a PR with a clear description of changes and tests or examples.

## License

Please check the repository `LICENSE` at the project root (if present). If none exists, contact the maintainer.

## Contact

Maintainer: Nimieeee (original target repo: https://github.com/Nimieeee/docking-suite)
