"""
Pytest fixtures for GWAS MCP tests.
"""

import pytest
import tempfile
import os
from pathlib import Path


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_vcf(temp_dir):
    """Create a sample VCF file for testing."""
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100000\trs123\tA\tG\t100\tPASS\tAF=0.25
1\t200000\trs456\tC\tT\t100\tPASS\tAF=0.15
2\t100000\trs789\tG\tA\t100\tPASS\tAF=0.30
"""
    vcf_path = temp_dir / "test.vcf"
    vcf_path.write_text(vcf_content)
    return vcf_path


@pytest.fixture
def sample_sumstats(temp_dir):
    """Create sample GWAS summary statistics for testing."""
    sumstats_content = """CHR\tBP\tSNP\tA1\tA2\tBETA\tSE\tP\tN
1\t100000\trs123\tA\tG\t0.05\t0.01\t1e-6\t10000
1\t200000\trs456\tC\tT\t0.03\t0.01\t0.01\t10000
2\t100000\trs789\tG\tA\t0.08\t0.01\t5e-10\t10000
2\t200000\trs101\tT\tC\t-0.02\t0.01\t0.5\t10000
3\t100000\trs102\tA\tC\t0.04\t0.01\t1e-4\t10000
"""
    sumstats_path = temp_dir / "sumstats.txt"
    sumstats_path.write_text(sumstats_content)
    return sumstats_path


@pytest.fixture
def sample_phenotype(temp_dir):
    """Create sample phenotype file for testing."""
    pheno_content = """FID\tIID\tPHENO
1\t1\t0.5
2\t2\t1.2
3\t3\t-0.3
4\t4\t0.8
5\t5\t-0.1
"""
    pheno_path = temp_dir / "phenotype.txt"
    pheno_path.write_text(pheno_content)
    return pheno_path
