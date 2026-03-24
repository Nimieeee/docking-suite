import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import PandasTools, FilterCatalog, Descriptors, MolFromSmarts
import os

# --- Visual Setup ---
def setup_visuals():
    np.random.seed(42)
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['axes.labelsize'] = 16
    return ["#A20025", "#6C8EBF"]

# --- Core Functions ---

def load_compound_library(file_path, smiles_name='smiles', mol_col_name='mol'):
    """
    Loads a compound library from SDF or CSV.
    """
    print(f"Loading file: {file_path}")
    
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    _, ext = os.path.splitext(file_path)
    
    if ext.lower() == '.sdf':
        df = PandasTools.LoadSDF(file_path, smilesName=smiles_name, molColName=None)
        print("Adding RDKit molecule objects...")
        PandasTools.AddMoleculeColumnToFrame(df, smiles_name, mol_col_name)
    elif ext.lower() == '.csv':
        df = pd.read_csv(file_path)
        print("Adding RDKit molecule objects from SMILES...")
        PandasTools.AddMoleculeColumnToFrame(df, smiles_name, mol_col_name)
    else:
        raise ValueError(f"Unsupported file format: {ext}")

    print(f"\nLoaded {len(df)} compounds")
    print(f"Valid molecules: {df[mol_col_name].notnull().sum()}")
    print(f"Invalid molecules: {df[mol_col_name].isnull().sum()}")
    
    # Remove invalid molecules
    initial_len = len(df)
    df = df.dropna(subset=[mol_col_name])
    if len(df) < initial_len:
        print(f"Removed {initial_len - len(df)} invalid molecules")
        
    return df

def calculate_ro5_descriptors(df, mol_col='mol'):
    """Calculate properties for Lipinski's Rule of 5."""
    RO5_PROPS = ['ExactMolWt', 'NumHAcceptors', 'NumHDonors', 'MolLogP']
    
    def compute_descriptor(mol, desc_name):
        if mol is None: return None
        try: return getattr(Descriptors, desc_name)(mol)
        except: return None
        
    print("Calculating molecular descriptors...")
    for desc in RO5_PROPS:
        df[desc] = df[mol_col].apply(lambda x: compute_descriptor(x, desc))
        
    df = df.dropna(subset=RO5_PROPS)
    return df

def apply_lipinski_filter(df):
    """Filter molecules based on Lipinski's Rule of 5 (Strict: 0 violations)."""
    def count_ro5_violations(row):
        violations = 0
        if row['ExactMolWt'] > 500: violations += 1
        if row['MolLogP'] > 5: violations += 1
        if row['NumHDonors'] > 5: violations += 1
        if row['NumHAcceptors'] > 10: violations += 1
        return violations
        
    df['ro5_violations'] = df.apply(count_ro5_violations, axis=1)
    df['ro5_compliant'] = df['ro5_violations'] == 0
    
    compliant = df[df['ro5_compliant']]
    print(f"RO5 Compliant: {len(compliant)} / {len(df)} ({len(compliant)/len(df)*100:.1f}%)")
    return df, compliant

def apply_pains_filter(df, mol_col='mol'):
    """Apply RDKit's built-in PAINS filter."""
    print("Applying PAINS filter...")
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(params.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)
    
    df['PAINS_match'] = df[mol_col].apply(catalog.HasMatch)
    compliant = df[~df['PAINS_match']]
    
    print(f"PAINS Compliant: {len(compliant)} / {len(df)} ({len(compliant)/len(df)*100:.1f}%)")
    return df, compliant

def apply_brenk_filter(df, mol_col='mol'):
    """Apply RDKit's built-in BRENK filter."""
    print("Applying Brenk filter...")
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(params.FilterCatalogs.BRENK)
    catalog = FilterCatalog.FilterCatalog(params)
    
    df['BRENK_match'] = df[mol_col].apply(catalog.HasMatch)
    compliant = df[~df['BRENK_match']]
    
    print(f"Brenk Compliant: {len(compliant)} / {len(df)} ({len(compliant)/len(df)*100:.1f}%)")
    return df, compliant

def load_glaxo_alerts(file_path):
    """Load Glaxo structural alerts from CSV."""
    if not os.path.exists(file_path):
        print(f"Warning: Glaxo alerts file not found at {file_path}")
        return None
        
    try:
        alerts_df = pd.read_csv(file_path)
        alerts_df['ROMol'] = alerts_df.smarts.apply(MolFromSmarts)
        alerts_df = alerts_df.dropna(subset=['ROMol'])
        print(f"Loaded {len(alerts_df)} Glaxo alerts")
        return alerts_df
    except Exception as e:
        print(f"Error loading alerts: {e}")
        return None

def apply_glaxo_filters(df, alerts_df, mol_col='mol'):
    """Apply custom structural alerts."""
    if alerts_df is None:
        return df, df
        
    print("Applying Glaxo filters...")
    
    def check_match(mol):
        for _, alert in alerts_df.iterrows():
            if mol.HasSubstructMatch(alert.ROMol):
                return True
        return False
        
    df['GLAXO_match'] = df[mol_col].apply(check_match)
    compliant = df[~df['GLAXO_match']]
    
    print(f"Glaxo Compliant: {len(compliant)} / {len(df)} ({len(compliant)/len(df)*100:.1f}%)")
    return df, compliant

def save_compounds(df, output_path, mol_col='mol'):
    """Save dataframe to SDF or CSV."""
    print(f"Exporting to {output_path}...")
    
    # Clean up columns unrelated to output if needed, or keep all
    output_cols = [c for c in df.columns if c != mol_col]
    
    if output_path.endswith('.sdf'):
        PandasTools.WriteSDF(df, output_path, molColName=mol_col, properties=output_cols)
    elif output_path.endswith('.csv'):
        df[output_cols].to_csv(output_path, index=False)
    
    print("✓ Export complete")

def run_filtering_pipeline(input_path, output_path, glaxo_alerts_path='data/ch02/glaxo_structural_alerts.csv'):
    """
    Orchestrates the full filtering pipeline.
    """
    setup_visuals()
    
    # 1. Load
    df = load_compound_library(input_path)
    
    # 2. RO5
    df = calculate_ro5_descriptors(df)
    _, ro5_pass = apply_lipinski_filter(df)
    
    # 3. PAINS
    _, pains_pass = apply_pains_filter(ro5_pass)

    # 4. Brenk
    _, brenk_pass = apply_brenk_filter(pains_pass)
    
    # 5. Glaxo
    alerts = load_glaxo_alerts(glaxo_alerts_path)
    _, final_pass = apply_glaxo_filters(brenk_pass, alerts)
    
    # 5. Save
    print(f"\nFinal Survivor Count: {len(final_pass)}")
    save_compounds(final_pass, output_path)
    
    return final_pass
