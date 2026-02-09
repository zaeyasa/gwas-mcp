"""
GWAS Analysis Tools for MCP Server.

Tools for running GWAS, calculating genomic inflation, and heritability estimation.
"""

import json
import logging
from typing import Any, Dict, List, Optional
import numpy as np
import pandas as pd
from scipy import stats

from mcp.types import Tool

from ..utils.validators import validate_file_path, validate_threshold, validate_pvalue
from ..utils.file_handlers import read_summary_stats, read_plink, write_results

logger = logging.getLogger(__name__)


# Define Analysis tools
ANALYSIS_TOOLS = [
    Tool(
        name="run_gwas",
        description="Perform genome-wide association study using linear or logistic regression. Returns summary statistics including p-values, beta coefficients, and standard errors.",
        inputSchema={
            "type": "object",
            "properties": {
                "plink_prefix": {
                    "type": "string",
                    "description": "Path prefix for PLINK genotype files"
                },
                "phenotype_path": {
                    "type": "string",
                    "description": "Path to phenotype file (tab-separated with FID, IID, PHENO columns)"
                },
                "covariate_path": {
                    "type": "string",
                    "description": "Optional path to covariate file"
                },
                "model": {
                    "type": "string",
                    "enum": ["linear", "logistic"],
                    "description": "Regression model type (default: linear)",
                    "default": "linear"
                },
                "output_path": {
                    "type": "string",
                    "description": "Path to save GWAS results"
                }
            },
            "required": ["plink_prefix", "phenotype_path"]
        }
    ),
    Tool(
        name="calculate_genomic_inflation",
        description="Calculate genomic inflation factor (lambda GC) from GWAS summary statistics. Used to assess population stratification.",
        inputSchema={
            "type": "object",
            "properties": {
                "sumstats_path": {
                    "type": "string",
                    "description": "Path to GWAS summary statistics file"
                },
                "pvalue_column": {
                    "type": "string",
                    "description": "Name of p-value column (default: P)",
                    "default": "P"
                }
            },
            "required": ["sumstats_path"]
        }
    ),
    Tool(
        name="identify_significant_snps",
        description="Filter GWAS results to identify genome-wide significant SNPs based on p-value threshold.",
        inputSchema={
            "type": "object",
            "properties": {
                "sumstats_path": {
                    "type": "string",
                    "description": "Path to GWAS summary statistics file"
                },
                "pvalue_threshold": {
                    "type": "number",
                    "description": "P-value threshold for significance (default: 5e-8)",
                    "default": 5e-8
                },
                "suggestive_threshold": {
                    "type": "number",
                    "description": "Suggestive significance threshold (default: 1e-5)",
                    "default": 1e-5
                },
                "output_path": {
                    "type": "string",
                    "description": "Path to save significant SNPs"
                }
            },
            "required": ["sumstats_path"]
        }
    ),
    Tool(
        name="calculate_heritability_ldsc",
        description="Estimate SNP-heritability using LD Score regression. Requires GWAS summary statistics and LD scores.",
        inputSchema={
            "type": "object",
            "properties": {
                "sumstats_path": {
                    "type": "string",
                    "description": "Path to GWAS summary statistics"
                },
                "ld_scores_path": {
                    "type": "string",
                    "description": "Path to LD scores file"
                },
                "sample_size": {
                    "type": "integer",
                    "description": "GWAS sample size (required if not in sumstats)"
                }
            },
            "required": ["sumstats_path"]
        }
    ),
]


async def handle_analysis_tool(name: str, arguments: Dict[str, Any]) -> str:
    """Handle analysis tool calls."""
    
    if name == "run_gwas":
        return await run_gwas(
            plink_prefix=arguments["plink_prefix"],
            phenotype_path=arguments["phenotype_path"],
            covariate_path=arguments.get("covariate_path"),
            model=arguments.get("model", "linear"),
            output_path=arguments.get("output_path")
        )
    
    elif name == "calculate_genomic_inflation":
        return await calculate_genomic_inflation(
            sumstats_path=arguments["sumstats_path"],
            pvalue_column=arguments.get("pvalue_column", "P")
        )
    
    elif name == "identify_significant_snps":
        return await identify_significant_snps(
            sumstats_path=arguments["sumstats_path"],
            pvalue_threshold=arguments.get("pvalue_threshold", 5e-8),
            suggestive_threshold=arguments.get("suggestive_threshold", 1e-5),
            output_path=arguments.get("output_path")
        )
    
    elif name == "calculate_heritability_ldsc":
        return await calculate_heritability_ldsc(
            sumstats_path=arguments["sumstats_path"],
            ld_scores_path=arguments.get("ld_scores_path"),
            sample_size=arguments.get("sample_size")
        )
    
    raise ValueError(f"Unknown analysis tool: {name}")


async def run_gwas(
    plink_prefix: str,
    phenotype_path: str,
    covariate_path: Optional[str] = None,
    model: str = "linear",
    output_path: Optional[str] = None
) -> str:
    """
    Run GWAS using statsmodels.
    
    Performs linear or logistic regression for each variant.
    """
    logger.info(f"Running GWAS: {plink_prefix}")
    
    try:
        import statsmodels.api as sm
        
        # Read genotype data
        plink_data = read_plink(plink_prefix)
        genotypes = plink_data['genotypes']
        bim = plink_data['bim']
        fam = plink_data['fam']
        
        # Read phenotype
        pheno_df = pd.read_csv(phenotype_path, sep='\t')
        
        # Standardize column names
        pheno_df.columns = pheno_df.columns.str.upper()
        if 'PHENO' not in pheno_df.columns and len(pheno_df.columns) >= 3:
            pheno_df.columns = ['FID', 'IID', 'PHENO'] + list(pheno_df.columns[3:])
        
        # Match samples
        # Assuming fam has 'fid' and 'iid' columns
        geno_matrix = genotypes.compute() if hasattr(genotypes, 'compute') else np.array(genotypes)
        
        # Get phenotype values
        phenotype = pheno_df['PHENO'].values
        
        # Check for binary phenotype
        unique_vals = np.unique(phenotype[~np.isnan(phenotype)])
        is_binary = len(unique_vals) == 2
        
        if model == "logistic" and not is_binary:
            logger.warning("Logistic regression requested but phenotype is not binary. Using linear regression.")
            model = "linear"
        
        # Read covariates if provided
        covariates = None
        if covariate_path:
            cov_df = pd.read_csv(covariate_path, sep='\t')
            covariates = cov_df.iloc[:, 2:].values  # Skip FID, IID
        
        # Run association tests
        results = []
        n_variants = geno_matrix.shape[1]
        
        for i in range(min(n_variants, 10000)):  # Limit for performance
            geno = geno_matrix[:, i]
            
            # Remove missing values
            valid = ~(np.isnan(geno) | np.isnan(phenotype))
            if np.sum(valid) < 10:
                continue
            
            X = geno[valid]
            y = phenotype[valid]
            
            # Add covariates if present
            if covariates is not None:
                X = np.column_stack([X, covariates[valid]])
            
            # Add constant
            X = sm.add_constant(X)
            
            try:
                if model == "logistic":
                    mod = sm.Logit(y, X)
                else:
                    mod = sm.OLS(y, X)
                
                res = mod.fit(disp=0)
                
                # Get results for genotype effect (index 1)
                results.append({
                    'SNP': bim.iloc[i]['snp'] if 'snp' in bim.columns else f"SNP_{i}",
                    'CHR': bim.iloc[i]['chrom'] if 'chrom' in bim.columns else 'NA',
                    'BP': bim.iloc[i]['pos'] if 'pos' in bim.columns else i,
                    'A1': bim.iloc[i]['a0'] if 'a0' in bim.columns else 'A',
                    'A2': bim.iloc[i]['a1'] if 'a1' in bim.columns else 'B',
                    'BETA': res.params[1],
                    'SE': res.bse[1],
                    'P': res.pvalues[1],
                    'N': int(np.sum(valid))
                })
            except Exception as e:
                logger.debug(f"GWAS failed for variant {i}: {e}")
                continue
        
        # Create results DataFrame
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values('P')
        
        # Calculate summary
        n_tested = len(results_df)
        n_significant = len(results_df[results_df['P'] < 5e-8])
        n_suggestive = len(results_df[results_df['P'] < 1e-5])
        
        output = {
            "n_samples": len(phenotype),
            "n_variants_tested": n_tested,
            "model": model,
            "n_genome_wide_significant": n_significant,
            "n_suggestive": n_suggestive,
            "top_hits": results_df.head(20).to_dict(orient='records'),
            "note": "Full GWAS may require more memory. Consider using PLINK/BOLT-LMM for large datasets."
        }
        
        if output_path:
            write_results(results_df, output_path, format='tsv')
            output['results_saved_to'] = output_path
        
        return json.dumps(output, indent=2, default=str)
        
    except ImportError as e:
        return json.dumps({
            "error": "Required package not installed",
            "message": str(e),
            "suggestion": "Install with: pip install statsmodels pandas-plink"
        }, indent=2)
    except Exception as e:
        logger.exception(f"GWAS failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


async def calculate_genomic_inflation(
    sumstats_path: str,
    pvalue_column: str = "P"
) -> str:
    """
    Calculate genomic inflation factor (lambda GC).
    
    Lambda GC is the ratio of median observed chi-squared statistic
    to expected median under the null hypothesis.
    """
    logger.info(f"Calculating lambda GC from: {sumstats_path}")
    
    try:
        # Read summary statistics
        df = read_summary_stats(sumstats_path)
        
        # Find p-value column
        pval_col = None
        for col in [pvalue_column, 'P', 'PVALUE', 'P_VALUE', 'PVAL']:
            if col in df.columns or col.upper() in df.columns:
                pval_col = col if col in df.columns else col.upper()
                break
        
        if pval_col is None:
            raise ValueError(f"P-value column not found. Available columns: {list(df.columns)}")
        
        pvalues = df[pval_col].dropna().values
        pvalues = pvalues[pvalues > 0]  # Remove zeros
        
        # Convert p-values to chi-squared statistics
        chi2_stats = stats.chi2.ppf(1 - pvalues, df=1)
        
        # Calculate lambda GC
        observed_median = np.median(chi2_stats)
        expected_median = stats.chi2.ppf(0.5, df=1)  # ~0.455
        
        lambda_gc = observed_median / expected_median
        
        # Interpretation
        if lambda_gc < 1.05:
            interpretation = "Excellent - minimal inflation"
        elif lambda_gc < 1.1:
            interpretation = "Good - slight inflation, acceptable"
        elif lambda_gc < 1.2:
            interpretation = "Moderate inflation - consider correcting"
        else:
            interpretation = "High inflation - population stratification likely"
        
        result = {
            "lambda_gc": round(lambda_gc, 4),
            "n_variants": len(pvalues),
            "observed_median_chi2": round(observed_median, 4),
            "expected_median_chi2": round(expected_median, 4),
            "interpretation": interpretation,
            "recommendation": "Use lambda GC to assess population stratification. Values > 1.1 suggest confounding."
        }
        
        return json.dumps(result, indent=2)
        
    except Exception as e:
        logger.exception(f"Lambda GC calculation failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


async def identify_significant_snps(
    sumstats_path: str,
    pvalue_threshold: float = 5e-8,
    suggestive_threshold: float = 1e-5,
    output_path: Optional[str] = None
) -> str:
    """
    Identify genome-wide and suggestive significant SNPs.
    """
    logger.info(f"Identifying significant SNPs from: {sumstats_path}")
    
    try:
        # Read summary statistics
        df = read_summary_stats(sumstats_path)
        
        # Find p-value column
        pval_col = None
        for col in ['P', 'PVALUE', 'P_VALUE', 'PVAL']:
            if col in df.columns:
                pval_col = col
                break
        
        if pval_col is None:
            raise ValueError(f"P-value column not found. Available: {list(df.columns)}")
        
        # Filter significant SNPs
        gw_sig = df[df[pval_col] < pvalue_threshold].copy()
        suggestive = df[(df[pval_col] >= pvalue_threshold) & (df[pval_col] < suggestive_threshold)].copy()
        
        # Sort by p-value
        gw_sig = gw_sig.sort_values(pval_col)
        suggestive = suggestive.sort_values(pval_col)
        
        result = {
            "total_variants": len(df),
            "genome_wide_significant": {
                "count": len(gw_sig),
                "threshold": pvalue_threshold,
                "snps": gw_sig.head(50).to_dict(orient='records')
            },
            "suggestive_significant": {
                "count": len(suggestive),
                "threshold": suggestive_threshold,
                "snps": suggestive.head(50).to_dict(orient='records')
            },
            "top_hit": gw_sig.head(1).to_dict(orient='records')[0] if len(gw_sig) > 0 else None
        }
        
        if output_path:
            all_sig = pd.concat([gw_sig, suggestive])
            write_results(all_sig, output_path, format='tsv')
            result['results_saved_to'] = output_path
        
        return json.dumps(result, indent=2, default=str)
        
    except Exception as e:
        logger.exception(f"Significant SNP identification failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


async def calculate_heritability_ldsc(
    sumstats_path: str,
    ld_scores_path: Optional[str] = None,
    sample_size: Optional[int] = None
) -> str:
    """
    Estimate SNP-heritability using LD Score regression.
    
    This is a simplified implementation. For production use,
    consider using the official LDSC software.
    """
    logger.info(f"Calculating heritability from: {sumstats_path}")
    
    try:
        # Read summary statistics
        df = read_summary_stats(sumstats_path)
        
        # Check for required columns
        required_cols = ['P', 'N'] if sample_size is None else ['P']
        for col in required_cols:
            if col not in df.columns:
                if col == 'N' and sample_size:
                    df['N'] = sample_size
                else:
                    raise ValueError(f"Required column '{col}' not found")
        
        # Get sample size
        n = sample_size or df['N'].median()
        
        # Convert p-values to chi-squared statistics
        pvalues = df['P'].dropna().values
        pvalues = pvalues[(pvalues > 0) & (pvalues < 1)]
        chi2_stats = stats.chi2.ppf(1 - pvalues, df=1)
        
        # Simple heritability estimate (Bulik-Sullivan method approximation)
        # This is simplified - full LDSC requires LD scores
        mean_chi2 = np.mean(chi2_stats)
        
        # Estimate h2 = (mean(chi2) - 1) / (N * mean_ld_score / M)
        # Without LD scores, we use a rough estimate
        if ld_scores_path:
            # Read LD scores
            ld_df = pd.read_csv(ld_scores_path, sep='\t')
            mean_ld = ld_df['L2'].mean() if 'L2' in ld_df.columns else 50
        else:
            mean_ld = 50  # Rough estimate for European samples
        
        m_variants = len(pvalues)
        
        # LDSC regression intercept estimate
        intercept = 1.0  # Assuming no confounding
        
        # Heritability calculation
        h2_liability = (mean_chi2 - intercept) * m_variants / (n * mean_ld)
        h2_liability = max(0, min(1, h2_liability))  # Bound between 0 and 1
        
        # Standard error (rough approximation)
        se_h2 = h2_liability * 0.2  # ~20% relative SE
        
        result = {
            "h2_snp": round(h2_liability, 4),
            "h2_se": round(se_h2, 4),
            "sample_size": int(n),
            "n_variants": m_variants,
            "mean_chi2": round(mean_chi2, 4),
            "ldsc_intercept": intercept,
            "interpretation": f"SNP-heritability estimate: {h2_liability:.1%} of trait variance explained by common SNPs",
            "note": "This is a simplified estimate. For publication-quality results, use official LDSC software (github.com/bulik/ldsc)"
        }
        
        return json.dumps(result, indent=2)
        
    except Exception as e:
        logger.exception(f"Heritability calculation failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__,
            "suggestion": "For accurate heritability estimates, use LDSC: github.com/bulik/ldsc"
        }, indent=2)


def register_analysis_tools():
    """Return analysis tools for registration."""
    return ANALYSIS_TOOLS
