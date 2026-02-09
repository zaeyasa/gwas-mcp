"""
Tests for validators module.
"""

import pytest
from pathlib import Path
import sys

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from gwas_mcp.utils.validators import (
    validate_rsid,
    validate_genomic_region,
    validate_pvalue,
    validate_chromosome,
    validate_threshold
)


class TestValidateRsid:
    """Tests for rsID validation."""
    
    def test_valid_rsid(self):
        assert validate_rsid("rs12345") == "rs12345"
        assert validate_rsid("rs1") == "rs1"
        assert validate_rsid("RS12345") == "rs12345"
    
    def test_invalid_rsid(self):
        with pytest.raises(ValueError):
            validate_rsid("12345")
        with pytest.raises(ValueError):
            validate_rsid("snp123")
        with pytest.raises(ValueError):
            validate_rsid("")


class TestValidateGenomicRegion:
    """Tests for genomic region validation."""
    
    def test_valid_regions(self):
        assert validate_genomic_region("chr1:1000-2000") == ("1", 1000, 2000)
        assert validate_genomic_region("1:1000-2000") == ("1", 1000, 2000)
        assert validate_genomic_region("chrX:1000-2000") == ("X", 1000, 2000)
    
    def test_invalid_regions(self):
        with pytest.raises(ValueError):
            validate_genomic_region("chr1:2000-1000")  # start > end
        with pytest.raises(ValueError):
            validate_genomic_region("chr1:1000")  # missing end
        with pytest.raises(ValueError):
            validate_genomic_region("invalid")


class TestValidatePvalue:
    """Tests for p-value validation."""
    
    def test_valid_pvalues(self):
        assert validate_pvalue(0.05) == 0.05
        assert validate_pvalue(1e-8) == 1e-8
        assert validate_pvalue(1.0) == 1.0
    
    def test_invalid_pvalues(self):
        with pytest.raises(ValueError):
            validate_pvalue(-0.1)  # negative
        with pytest.raises(ValueError):
            validate_pvalue(1.5)  # > 1


class TestValidateChromosome:
    """Tests for chromosome validation."""
    
    def test_valid_chromosomes(self):
        assert validate_chromosome("1") == "1"
        assert validate_chromosome("chr1") == "1"
        assert validate_chromosome("X") == "X"
        assert validate_chromosome("chrX") == "X"
        assert validate_chromosome("MT") == "MT"
    
    def test_invalid_chromosomes(self):
        with pytest.raises(ValueError):
            validate_chromosome("chr25")
        with pytest.raises(ValueError):
            validate_chromosome("invalid")


class TestValidateThreshold:
    """Tests for numeric threshold validation."""
    
    def test_valid_threshold(self):
        assert validate_threshold(0.5, "test", 0, 1) == 0.5
        assert validate_threshold(100, "test", 0, 1000) == 100
    
    def test_threshold_out_of_range(self):
        with pytest.raises(ValueError):
            validate_threshold(-1, "test", 0, 1)
        with pytest.raises(ValueError):
            validate_threshold(2, "test", 0, 1)
