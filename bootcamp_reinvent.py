"""
Bootcamp REINVENT 4 Integration
Wrapper for REINVENT 4 generative molecular design via CLI
Based on: https://github.com/MolecularAI/REINVENT4

INSTALLATION (Colab/Local):
    git clone https://github.com/MolecularAI/REINVENT4.git
    cd REINVENT4 && pip install --no-deps .
"""

import os
import subprocess
import tempfile
import pandas as pd
from typing import Optional, Dict, List
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.QED import qed as calculate_qed

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Default prior model URL (Reinvent de novo prior from Zenodo)
PRIOR_URL = "https://zenodo.org/records/15641297/files/reinvent.prior?download=1"
DEFAULT_PRIOR_PATH = "priors/reinvent.prior"

# ==============================================================================
# PRIOR MODEL HANDLING
# ==============================================================================

def download_prior(url: str = PRIOR_URL, save_path: str = DEFAULT_PRIOR_PATH) -> str:
    """
    Download REINVENT prior model from Zenodo if not present.
    
    Args:
        url: URL to download from
        save_path: Local path to save the model
    
    Returns:
        Path to the downloaded model
    """
    save_path = Path(save_path)
    
    if save_path.exists():
        print(f"✓ Prior model already exists: {save_path}")
        return str(save_path)
    
    print(f"📥 Downloading REINVENT prior model (~400MB)...")
    print(f"   From: {url}")
    
    save_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Use wget or curl for download
    try:
        subprocess.run(
            ['wget', '-q', '--show-progress', '-O', str(save_path), url],
            check=True
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Try curl if wget not available
        subprocess.run(
            ['curl', '-L', '-o', str(save_path), url],
            check=True
        )
    
    print(f"✓ Downloaded to: {save_path}")
    return str(save_path)

def check_reinvent_installed() -> bool:
    """Check if REINVENT 4 is installed and accessible."""
    try:
        result = subprocess.run(
            ['reinvent', '--help'],
            capture_output=True,
            text=True,
            timeout=10
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False

# ==============================================================================
# TOML CONFIG GENERATION
# ==============================================================================

def generate_sampling_config(
    model_file: str,
    output_file: str,
    num_smiles: int = 100,
    unique_molecules: bool = True,
    randomize_smiles: bool = True,
    device: str = "cuda:0"
) -> str:
    """
    Generate a TOML config for REINVENT sampling.
    
    Returns:
        String containing the TOML config
    """
    config = f'''# REINVENT4 Sampling Config (Auto-generated)

run_type = "sampling"
device = "{device}"

[parameters]
model_file = "{model_file}"
output_file = "{output_file}"
num_smiles = {num_smiles}
unique_molecules = {str(unique_molecules).lower()}
randomize_smiles = {str(randomize_smiles).lower()}
'''
    return config

def generate_rl_config(
    model_file: str,
    output_dir: str,
    agent_file: str,
    num_steps: int = 50,
    batch_size: int = 64,
    device: str = "cuda:0"
) -> str:
    """
    Generate a TOML config for REINVENT RL optimization.
    
    Returns:
        String containing the TOML config
    """
    config = f'''# REINVENT4 RL Config (Auto-generated)

run_type = "reinforcement_learning"
device = "{device}"

[parameters]
prior_file = "{model_file}"
agent_file = "{agent_file}"
summary_csv_prefix = "{output_dir}/rl_summary"

[learning_strategy]
type = "dap"
sigma = 128

[stage.1]
max_steps = {num_steps}
batch_size = {batch_size}
chkpt_file = "{output_dir}/agent_checkpoint.chkpt"

    [stage.1.scoring]
    type = "standard"

        [[stage.1.scoring.component]]
        component_type = "custom_product"
        name = "QED_score"

            [[stage.1.scoring.component.custom_alerts]]
            endpoint = "qed"
            weight = 1.0
'''
    return config

# ==============================================================================
# REINVENT EXECUTION
# ==============================================================================

def run_reinvent_sampling(
    n_samples: int = 100,
    model_file: Optional[str] = None,
    device: str = "cuda:0",
    auto_download: bool = True
) -> pd.DataFrame:
    """
    Run REINVENT 4 sampling to generate molecules.
    
    Args:
        n_samples: Number of molecules to generate
        model_file: Path to prior model (downloads default if None)
        device: Torch device (cuda:0, cpu)
        auto_download: Whether to auto-download prior if missing
    
    Returns:
        DataFrame with generated SMILES and NLL scores
    """
    if not check_reinvent_installed():
        raise RuntimeError(
            "REINVENT 4 not installed. Install with:\n"
            "  git clone https://github.com/MolecularAI/REINVENT4.git\n"
            "  cd REINVENT4 && pip install --no-deps ."
        )
    
    # Handle model file
    if model_file is None:
        if auto_download:
            model_file = download_prior()
        else:
            raise ValueError("model_file must be provided or auto_download=True")
    
    if not os.path.exists(model_file):
        raise FileNotFoundError(f"Model file not found: {model_file}")
    
    # Create temp files for config and output
    with tempfile.TemporaryDirectory() as tmpdir:
        config_file = os.path.join(tmpdir, "sampling.toml")
        output_file = os.path.join(tmpdir, "output.csv")
        log_file = os.path.join(tmpdir, "reinvent.log")
        
        # Generate config
        config = generate_sampling_config(
            model_file=model_file,
            output_file=output_file,
            num_smiles=n_samples,
            device=device
        )
        
        with open(config_file, 'w') as f:
            f.write(config)
        
        print(f"🚀 Running REINVENT sampling ({n_samples} molecules)...")
        
        # Run REINVENT
        result = subprocess.run(
            ['reinvent', '-l', log_file, config_file],
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            print(f"⚠️ REINVENT error: {result.stderr}")
            # Try to read log for more info
            if os.path.exists(log_file):
                with open(log_file) as f:
                    print(f"Log: {f.read()[-500:]}")  # Last 500 chars
            raise RuntimeError("REINVENT sampling failed")
        
        # Read output
        if not os.path.exists(output_file):
            raise RuntimeError(f"REINVENT did not produce output file: {output_file}")
        
        df = pd.read_csv(output_file)
        print(f"✓ Generated {len(df)} molecules")
        
        return df

# ==============================================================================
# SCORING FUNCTIONS
# ==============================================================================

def qed_score(smiles: str) -> float:
    """Quantitative Estimate of Drug-likeness."""
    mol = Chem.MolFromSmiles(smiles)
    return calculate_qed(mol) if mol else 0.0

def mw_score(smiles: str, target: float = 350, tolerance: float = 100) -> float:
    """Molecular weight penalty (centered on target)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0
    mw = Descriptors.MolWt(mol)
    deviation = abs(mw - target)
    return max(0, 1 - deviation / tolerance)

def logp_score(smiles: str, min_val: float = 1.0, max_val: float = 4.0) -> float:
    """LogP within desired range."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0
    logp = Descriptors.MolLogP(mol)
    if logp < min_val:
        return max(0, 1 - (min_val - logp))
    elif logp > max_val:
        return max(0, 1 - (logp - max_val))
    return 1.0

def combined_score(smiles: str, weights: Dict[str, float] = None) -> float:
    """Multi-objective scoring function."""
    if weights is None:
        weights = {'qed': 0.5, 'mw': 0.25, 'logp': 0.25}
    
    scores = {
        'qed': qed_score(smiles),
        'mw': mw_score(smiles),
        'logp': logp_score(smiles),
    }
    
    return sum(weights.get(k, 0) * v for k, v in scores.items())

# ==============================================================================
# HIGH-LEVEL API
# ==============================================================================

def run_reinvent_pipeline(
    n_generate: int = 100,
    model_path: Optional[str] = None,
    reference_smiles: Optional[str] = None,
    device: str = "cuda:0"
) -> Dict:
    """
    Complete REINVENT generation pipeline.
    
    Args:
        n_generate: Number of molecules to generate
        model_path: Path to REINVENT prior model
        reference_smiles: Optional reference for similarity (future use)
        device: Torch device
    
    Returns:
        Dictionary with generated molecules and scores
    """
    print("=" * 60)
    print("REINVENT 4 GENERATION PIPELINE")
    print("=" * 60)
    
    # Check REINVENT installation
    if not check_reinvent_installed():
        print("⚠️ REINVENT 4 not found. Installing...")
        # Attempt installation
        try:
            subprocess.run(
                ['pip', 'install', '-q', 'git+https://github.com/MolecularAI/REINVENT4.git'],
                check=True
            )
        except subprocess.CalledProcessError:
            raise RuntimeError(
                "Failed to install REINVENT 4. Please install manually:\n"
                "  git clone https://github.com/MolecularAI/REINVENT4.git\n"
                "  cd REINVENT4 && pip install ."
            )
    
    # Run sampling
    df = run_reinvent_sampling(
        n_samples=n_generate,
        model_file=model_path,
        device=device
    )
    
    # Add scores
    print("\n📊 Scoring generated molecules...")
    df['qed'] = df['SMILES'].apply(qed_score)
    df['mw'] = df['SMILES'].apply(lambda s: Descriptors.MolWt(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else 0)
    df['logp'] = df['SMILES'].apply(lambda s: Descriptors.MolLogP(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else 0)
    df['combined_score'] = df['SMILES'].apply(combined_score)
    
    # Sort by score
    df = df.sort_values('combined_score', ascending=False)
    
    print(f"\n✓ Top molecule: {df.iloc[0]['SMILES']}")
    print(f"  QED: {df.iloc[0]['qed']:.3f}, Score: {df.iloc[0]['combined_score']:.3f}")
    print("=" * 60)
    
    return {
        'all_results': df,
        'top_10': df.head(10),
        'optimized': None,  # Placeholder for RL results
        'reinvent_available': True
    }
