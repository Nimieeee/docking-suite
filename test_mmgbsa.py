import pytest
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestMMGBSACalculator:
    """Tests for the MMGBSACalculator class."""
    
    @pytest.fixture
    def sample_protein_pdb(self, tmp_path):
        """Create a minimal test protein PDB file."""
        pdb_content = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C  
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C  
ATOM      4  O   ALA A   1       1.246   2.390   0.000  1.00  0.00           O  
ATOM      5  CB  ALA A   1       1.986  -0.760  -1.225  1.00  0.00           C  
END
"""
        pdb_file = tmp_path / "test_protein.pdb"
        pdb_file.write_text(pdb_content)
        return str(pdb_file)
    
    @pytest.fixture
    def sample_ligand_sdf(self, tmp_path):
        """Create a minimal test ligand SDF file."""
        sdf_content = """     RDKit          3D

  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
"""
        sdf_file = tmp_path / "test_ligand.sdf"
        sdf_file.write_text(sdf_content)
        return str(sdf_file)
    
    def test_import_mmgbsa_calculator(self):
        """Test that MMGBSACalculator can be imported."""
        from docking_utils import MMGBSACalculator
        assert MMGBSACalculator is not None
    
    def test_mmgbsa_calculator_init(self):
        """Test MMGBSACalculator initialization."""
        from docking_utils import MMGBSACalculator
        calc = MMGBSACalculator(implicit_solvent='OBC2')
        assert calc.implicit_solvent == 'OBC2'
    
    def test_calculate_score_returns_float(self, sample_protein_pdb, sample_ligand_sdf):
        """Test that calculate_score returns a numerical value."""
        try:
            from docking_utils import MMGBSACalculator
            calc = MMGBSACalculator()
            score = calc.calculate_score(sample_protein_pdb, sample_ligand_sdf)
            assert isinstance(score, float)
        except ImportError:
            pytest.skip("OpenMM or OpenForceField not available")


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
