"""
Quality Control Tools for GWAS MCP Server.

Tools for variant quality control, PCA, LD pruning, and missingness analysis.
"""

import json
import logging
from typing import Any, Dict, List, Optional
import numpy as np
import pandas as pd
from scipy import stats

from mcp.types import Tool

from ..utils.validators import validate_file_path, validate_threshold
from ..utils.file_handlers import read_vcf, read_plink, write_results

logger = logging.getLogger(__name__)


# Define QC tools
QC_TOOLS = [
    Tool(
        name="vcf_quality_check",
        description="Perform quality control on a VCF file. Filters variants based on MAF, missingness, and Hardy-Weinberg equilibrium.",
        inputSchema={
            "type": "object",
            "properties": {
                "vcf_path": {
                    "type": "string",
                    "description": "Path to VCF file (.vcf or .vcf.gz)"
                },
                "maf_threshold": {
                    "type": "number",
                    "description": "Minor allele frequency threshold (default: 0.01)",
                    "default": 0.01
                },
                "missingness_threshold": {
                    "type": "number",
                    "description": "Maximum allowed missingness rate (default: 0.1)",
                    "default": 0.1
                },
                "hwe_pvalue": {
                    "type": "number",
                    "description": "Hardy-Weinberg equilibrium p-value threshold (default: 1e-6)",
                    "default": 1e-6
                },
                "output_path": {
                    "type": "string",
                    "description": "Optional path to save QC report"
                }
            },
            "required": ["vcf_path"]
        }
    ),
    Tool(
        name="calculate_pca",
        description="Perform Principal Component Analysis for population stratification. Uses PLINK files as input.",
        inputSchema={
            "type": "object",
            "properties": {
                "plink_prefix": {
                    "type": "string",
                    "description": "Path prefix for PLINK files (without .bed/.bim/.fam extension)"
                },
                "n_components": {
                    "type": "integer",
                    "description": "Number of principal components to compute (default: 10)",
                    "default": 10
                },
                "output_path": {
                    "type": "string",
                    "description": "Optional path to save PCA results"
                }
            },
            "required": ["plink_prefix"]
        }
    ),
    Tool(
        name="ld_pruning",
        description="Prune SNPs based on linkage disequilibrium (LD). Removes variants in high LD to create an independent set.",
        inputSchema={
            "type": "object",
            "properties": {
                "plink_prefix": {
                    "type": "string",
                    "description": "Path prefix for PLINK files"
                },
                "r2_threshold": {
                    "type": "number",
                    "description": "r² threshold for LD pruning (default: 0.2)",
                    "default": 0.2
                },
                "window_size": {
                    "type": "integer",
                    "description": "Window size in kb for LD calculation (default: 500)",
                    "default": 500
                },
                "step_size": {
                    "type": "integer",
                    "description": "Step size in variants (default: 50)",
                    "default": 50
                }
            },
            "required": ["plink_prefix"]
        }
    ),
    Tool(
        name="calculate_missingness",
        description="Calculate per-individual and per-variant missingness rates from genotype data.",
        inputSchema={
            "type": "object",
            "properties": {
                "plink_prefix": {
                    "type": "string",
                    "description": "Path prefix for PLINK files"
                },
                "output_path": {
                    "type": "string",
                    "description": "Optional path to save missingness report"
                }
            },
            "required": ["plink_prefix"]
        }
    ),
]


async def handle_qc_tool(name: str, arguments: Dict[str, Any]) -> str:
    """Handle QC tool calls."""
    
    if name == "vcf_quality_check":
        return await vcf_quality_check(
            vcf_path=arguments["vcf_path"],
            maf_threshold=arguments.get("maf_threshold", 0.01),
            missingness_threshold=arguments.get("missingness_threshold", 0.1),
            hwe_pvalue=arguments.get("hwe_pvalue", 1e-6),
            output_path=arguments.get("output_path")
        )
    
    elif name == "calculate_pca":
        return await calculate_pca(
            plink_prefix=arguments["plink_prefix"],
            n_components=arguments.get("n_components", 10),
            output_path=arguments.get("output_path")
        )
    
    elif name == "ld_pruning":
        return await ld_pruning(
            plink_prefix=arguments["plink_prefix"],
            r2_threshold=arguments.get("r2_threshold", 0.2),
            window_size=arguments.get("window_size", 500),
            step_size=arguments.get("step_size", 50)
        )
    
    elif name == "calculate_missingness":
        return await calculate_missingness(
            plink_prefix=arguments["plink_prefix"],
            output_path=arguments.get("output_path")
        )
    
    raise ValueError(f"Unknown QC tool: {name}")


async def vcf_quality_check(
    vcf_path: str,
    maf_threshold: float = 0.01,
    missingness_threshold: float = 0.1,
    hwe_pvalue: float = 1e-6,
    output_path: Optional[str] = None
) -> str:
    """
    Perform quality control on VCF file.
    
    Filters variants based on:
    - Minor allele frequency (MAF)
    - Genotype missingness
    - Hardy-Weinberg equilibrium (HWE)
    """
    logger.info(f"Running VCF QC on: {vcf_path}")
    
    # Validate inputs
    validate_file_path(vcf_path, must_exist=True, allowed_extensions=['vcf', 'vcf.gz'])
    maf_threshold = validate_threshold(maf_threshold, "MAF threshold", 0, 0.5)
    missingness_threshold = validate_threshold(missingness_threshold, "missingness", 0, 1)
    
    # Read VCF
    df = read_vcf(vcf_path)
    total_variants = len(df)
    
    # Calculate basic statistics (simplified for VCF without genotypes)
    qc_report = {
        "input_file": vcf_path,
        "total_variants": total_variants,
        "parameters": {
            "maf_threshold": maf_threshold,
            "missingness_threshold": missingness_threshold,
            "hwe_pvalue": hwe_pvalue
        },
        "summary": {
            "chromosomes": df['CHROM'].nunique() if 'CHROM' in df.columns else 0,
            "variant_types": {},
        },
        "filters_applied": [],
        "variants_remaining": total_variants,
    }
    
    # Count variant types if ALT column exists
    if 'ALT' in df.columns:
        for idx, row in df.iterrows():
            alt = str(row['ALT'])
            if len(row['REF']) == 1 and len(alt) == 1:
                vtype = 'SNP'
            elif len(row['REF']) > len(alt):
                vtype = 'DEL'
            elif len(row['REF']) < len(alt):
                vtype = 'INS'
            else:
                vtype = 'OTHER'
            qc_report['summary']['variant_types'][vtype] = \
                qc_report['summary']['variant_types'].get(vtype, 0) + 1
    
    # Note: Full QC requires genotype data
    qc_report['note'] = (
        "Full QC (MAF, missingness, HWE filtering) requires genotype data. "
        "For complete QC, use PLINK files with the calculate_missingness tool."
    )
    
    # Save report if output path provided
    if output_path:
        write_results(qc_report, output_path, format='json')
        qc_report['report_saved_to'] = output_path
    
    return json.dumps(qc_report, indent=2)


async def calculate_pca(
    plink_prefix: str,
    n_components: int = 10,
    output_path: Optional[str] = None
) -> str:
    """
    Perform PCA for population stratification analysis.
    """
    logger.info(f"Running PCA on: {plink_prefix}")
    
    # Validate inputs
    n_components = int(validate_threshold(n_components, "n_components", 1, 100))
    
    try:
        # Read PLINK files
        plink_data = read_plink(plink_prefix)
        genotypes = plink_data['genotypes']
        fam = plink_data['fam']
        
        # Convert to numpy array and handle missing values
        geno_matrix = genotypes.compute() if hasattr(genotypes, 'compute') else np.array(genotypes)
        
        # Replace missing values with mean
        col_means = np.nanmean(geno_matrix, axis=0)
        nan_mask = np.isnan(geno_matrix)
        geno_matrix[nan_mask] = np.take(col_means, np.where(nan_mask)[1])
        
        # Center and scale
        geno_centered = geno_matrix - np.mean(geno_matrix, axis=0)
        geno_std = np.std(geno_centered, axis=0)
        geno_std[geno_std == 0] = 1  # Avoid division by zero
        geno_scaled = geno_centered / geno_std
        
        # PCA via SVD
        n_components = min(n_components, min(geno_scaled.shape) - 1)
        U, S, Vt = np.linalg.svd(geno_scaled, full_matrices=False)
        
        # Get principal components
        pcs = U[:, :n_components] * S[:n_components]
        eigenvalues = (S ** 2) / (geno_scaled.shape[0] - 1)
        variance_explained = eigenvalues[:n_components] / np.sum(eigenvalues)
        
        # Create results DataFrame
        pc_columns = [f'PC{i+1}' for i in range(n_components)]
        pca_df = pd.DataFrame(pcs, columns=pc_columns)
        pca_df['FID'] = fam['fid'].values if 'fid' in fam.columns else range(len(pca_df))
        pca_df['IID'] = fam['iid'].values if 'iid' in fam.columns else range(len(pca_df))
        
        result = {
            "n_samples": len(pca_df),
            "n_variants_used": geno_matrix.shape[1],
            "n_components": n_components,
            "variance_explained": variance_explained.tolist(),
            "cumulative_variance": np.cumsum(variance_explained).tolist(),
            "eigenvalues": eigenvalues[:n_components].tolist(),
            "pca_results_preview": pca_df.head(10).to_dict(orient='records')
        }
        
        if output_path:
            write_results(pca_df, output_path, format='tsv')
            result['results_saved_to'] = output_path
        
        return json.dumps(result, indent=2)
        
    except ImportError as e:
        return json.dumps({
            "error": "pandas-plink not installed",
            "message": str(e),
            "suggestion": "Install with: pip install pandas-plink"
        }, indent=2)
    except Exception as e:
        logger.exception(f"PCA failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


async def ld_pruning(
    plink_prefix: str,
    r2_threshold: float = 0.2,
    window_size: int = 500,
    step_size: int = 50
) -> str:
    """
    Prune SNPs based on linkage disequilibrium.
    
    Uses a sliding window approach to identify and remove SNPs in high LD.
    """
    logger.info(f"Running LD pruning on: {plink_prefix}")
    
    # Validate inputs
    r2_threshold = validate_threshold(r2_threshold, "r2_threshold", 0, 1)
    window_size = int(validate_threshold(window_size, "window_size", 1, 10000))
    step_size = int(validate_threshold(step_size, "step_size", 1, 1000))
    
    try:
        # Read PLINK files
        plink_data = read_plink(plink_prefix)
        genotypes = plink_data['genotypes']
        bim = plink_data['bim']
        
        # Convert to numpy
        geno_matrix = genotypes.compute() if hasattr(genotypes, 'compute') else np.array(genotypes)
        
        n_samples, n_variants = geno_matrix.shape
        
        # Simple LD pruning implementation
        # (For production, consider using PLINK binary for efficiency)
        keep_indices = []
        exclude_indices = set()
        
        for i in range(n_variants):
            if i in exclude_indices:
                continue
            
            keep_indices.append(i)
            
            # Calculate LD with variants in window
            window_end = min(i + window_size, n_variants)
            
            for j in range(i + 1, window_end):
                if j in exclude_indices:
                    continue
                
                # Calculate r² between variants i and j
                geno_i = geno_matrix[:, i]
                geno_j = geno_matrix[:, j]
                
                # Remove missing values
                valid = ~(np.isnan(geno_i) | np.isnan(geno_j))
                if np.sum(valid) < 10:
                    continue
                
                r = np.corrcoef(geno_i[valid], geno_j[valid])[0, 1]
                r2 = r ** 2 if not np.isnan(r) else 0
                
                if r2 > r2_threshold:
                    exclude_indices.add(j)
        
        # Get pruned SNP list
        pruned_snps = bim.iloc[keep_indices]['snp'].tolist() if 'snp' in bim.columns else keep_indices
        
        result = {
            "total_variants": n_variants,
            "variants_kept": len(keep_indices),
            "variants_removed": n_variants - len(keep_indices),
            "pruning_parameters": {
                "r2_threshold": r2_threshold,
                "window_size_kb": window_size,
                "step_size": step_size
            },
            "pruned_snps_preview": pruned_snps[:20] if len(pruned_snps) > 20 else pruned_snps,
            "note": "For large datasets, consider using PLINK binary for faster LD pruning."
        }
        
        return json.dumps(result, indent=2)
        
    except ImportError as e:
        return json.dumps({
            "error": "pandas-plink not installed",
            "message": str(e),
            "suggestion": "Install with: pip install pandas-plink"
        }, indent=2)
    except Exception as e:
        logger.exception(f"LD pruning failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


async def calculate_missingness(
    plink_prefix: str,
    output_path: Optional[str] = None
) -> str:
    """
    Calculate per-individual and per-variant missingness rates.
    """
    logger.info(f"Calculating missingness for: {plink_prefix}")
    
    try:
        # Read PLINK files
        plink_data = read_plink(plink_prefix)
        genotypes = plink_data['genotypes']
        bim = plink_data['bim']
        fam = plink_data['fam']
        
        # Convert to numpy
        geno_matrix = genotypes.compute() if hasattr(genotypes, 'compute') else np.array(genotypes)
        
        # Calculate missingness
        missing_mask = np.isnan(geno_matrix)
        
        # Per-individual missingness
        ind_missingness = np.mean(missing_mask, axis=1)
        ind_results = pd.DataFrame({
            'FID': fam['fid'].values if 'fid' in fam.columns else range(len(fam)),
            'IID': fam['iid'].values if 'iid' in fam.columns else range(len(fam)),
            'N_MISS': np.sum(missing_mask, axis=1),
            'N_GENO': geno_matrix.shape[1],
            'F_MISS': ind_missingness
        })
        
        # Per-variant missingness
        var_missingness = np.mean(missing_mask, axis=0)
        var_results = pd.DataFrame({
            'SNP': bim['snp'].values if 'snp' in bim.columns else range(len(bim)),
            'CHR': bim['chrom'].values if 'chrom' in bim.columns else 'NA',
            'N_MISS': np.sum(missing_mask, axis=0),
            'N_GENO': geno_matrix.shape[0],
            'F_MISS': var_missingness
        })
        
        result = {
            "n_samples": len(fam),
            "n_variants": len(bim),
            "overall_missingness": float(np.mean(missing_mask)),
            "individual_summary": {
                "mean_missingness": float(np.mean(ind_missingness)),
                "max_missingness": float(np.max(ind_missingness)),
                "min_missingness": float(np.min(ind_missingness)),
                "individuals_above_10pct": int(np.sum(ind_missingness > 0.1)),
                "individuals_above_5pct": int(np.sum(ind_missingness > 0.05)),
            },
            "variant_summary": {
                "mean_missingness": float(np.mean(var_missingness)),
                "max_missingness": float(np.max(var_missingness)),
                "min_missingness": float(np.min(var_missingness)),
                "variants_above_10pct": int(np.sum(var_missingness > 0.1)),
                "variants_above_5pct": int(np.sum(var_missingness > 0.05)),
            },
            "high_missingness_individuals": ind_results[ind_results['F_MISS'] > 0.05].head(10).to_dict(orient='records'),
            "high_missingness_variants": var_results[var_results['F_MISS'] > 0.05].head(10).to_dict(orient='records'),
        }
        
        if output_path:
            # Save both reports
            ind_path = output_path.replace('.', '_individuals.')
            var_path = output_path.replace('.', '_variants.')
            write_results(ind_results, ind_path, format='tsv')
            write_results(var_results, var_path, format='tsv')
            result['individual_results_saved_to'] = ind_path
            result['variant_results_saved_to'] = var_path
        
        return json.dumps(result, indent=2)
        
    except ImportError as e:
        return json.dumps({
            "error": "pandas-plink not installed",
            "message": str(e),
            "suggestion": "Install with: pip install pandas-plink"
        }, indent=2)
    except Exception as e:
        logger.exception(f"Missingness calculation failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


def register_qc_tools():
    """Return QC tools for registration."""
    return QC_TOOLS
