"""
Tests for analysis tools.
"""

import pytest
import json
from pathlib import Path
import sys

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


class TestGenomicInflation:
    """Tests for genomic inflation calculation."""
    
    @pytest.mark.asyncio
    async def test_calculate_lambda_gc(self, sample_sumstats):
        from gwas_mcp.tools.analysis_tools import calculate_genomic_inflation
        
        result = await calculate_genomic_inflation(str(sample_sumstats))
        result_dict = json.loads(result)
        
        assert "lambda_gc" in result_dict
        assert "n_variants" in result_dict
        assert result_dict["n_variants"] == 5


class TestSignificantSnps:
    """Tests for significant SNP identification."""
    
    @pytest.mark.asyncio
    async def test_identify_significant(self, sample_sumstats):
        from gwas_mcp.tools.analysis_tools import identify_significant_snps
        
        result = await identify_significant_snps(
            str(sample_sumstats),
            pvalue_threshold=5e-8,
            suggestive_threshold=1e-5
        )
        result_dict = json.loads(result)
        
        assert "genome_wide_significant" in result_dict
        assert "suggestive_significant" in result_dict
        # rs789 should be genome-wide significant (p=5e-10)
        assert result_dict["genome_wide_significant"]["count"] == 1
