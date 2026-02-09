"""
Post-GWAS Analysis Tools for MCP Server.

Tools for polygenic risk scores, LD clumping, and locus extraction.
"""

import json
import logging
from typing import Any, Dict, List, Optional
import numpy as np
import pandas as pd

from mcp.types import Tool

from ..utils.validators import validate_file_path, validate_threshold
from ..utils.file_handlers import read_summary_stats, read_plink, write_results

logger = logging.getLogger(__name__)


# Define Post-GWAS tools
POST_GWAS_TOOLS = [
    Tool(
        name="calculate_prs",
        description="Calculate Polygenic Risk Scores (PRS) for individuals using GWAS summary statistics and target genotypes.",
        inputSchema={
            "type": "object",
            "properties": {
                "sumstats_path": {
                    "type": "string",
                    "description": "Path to GWAS summary statistics (base data)"
                },
                "target_plink_prefix": {
                    "type": "string",
                    "description": "Path prefix for target PLINK files"
                },
                "pvalue_threshold": {
                    "type": "number",
                    "description": "P-value threshold for SNP inclusion (default: 1.0 - all SNPs)",
                    "default": 1.0
                },
                "clump_r2": {
                    "type": "number",
                    "description": "LD r² threshold for clumping (default: 0.1)",
                    "default": 0.1
                },
                "output_path": {
                    "type": "string",
                    "description": "Path to save PRS results"
                }
            },
            "required": ["sumstats_path", "target_plink_prefix"]
        }
    ),
    Tool(
        name="clump_snps",
        description="Perform LD-based clumping of GWAS results to identify independent significant signals.",
        inputSchema={
            "type": "object",
            "properties": {
                "sumstats_path": {
                    "type": "string",
                    "description": "Path to GWAS summary statistics"
                },
                "plink_prefix": {
                    "type": "string",
                    "description": "Path prefix for PLINK files (LD reference)"
                },
                "p1": {
                    "type": "number",
                    "description": "P-value threshold for index SNPs (default: 5e-8)",
                    "default": 5e-8
                },
                "p2": {
                    "type": "number",
                    "description": "P-value threshold for clumped SNPs (default: 1e-4)",
                    "default": 1e-4
                },
                "r2": {
                    "type": "number",
                    "description": "LD r² threshold (default: 0.5)",
                    "default": 0.5
                },
                "kb": {
                    "type": "integer",
                    "description": "Clumping window in kb (default: 250)",
                    "default": 250
                },
                "output_path": {
                    "type": "string",
                    "description": "Path to save clumped results"
                }
            },
            "required": ["sumstats_path"]
        }
    ),
    Tool(
        name="extract_locus",
        description="Extract genotype data for a specific genomic region around a lead SNP for fine-mapping.",
        inputSchema={
            "type": "object",
            "properties": {
                "plink_prefix": {
                    "type": "string",
                    "description": "Path prefix for PLINK files"
                },
                "lead_snp": {
                    "type": "string",
                    "description": "Lead SNP identifier (rsID or CHR:POS)"
                },
                "window_kb": {
                    "type": "integer",
                    "description": "Window size in kb (default: 500)",
                    "default": 500
                },
                "output_prefix": {
                    "type": "string",
                    "description": "Output prefix for extracted files"
                }
            },
            "required": ["plink_prefix", "lead_snp"]
        }
    ),
]


async def handle_post_gwas_tool(name: str, arguments: Dict[str, Any]) -> str:
    """Handle post-GWAS tool calls."""
    
    if name == "calculate_prs":
        return await calculate_prs(
            sumstats_path=arguments["sumstats_path"],
            target_plink_prefix=arguments["target_plink_prefix"],
            pvalue_threshold=arguments.get("pvalue_threshold", 1.0),
            clump_r2=arguments.get("clump_r2", 0.1),
            output_path=arguments.get("output_path")
        )
    
    elif name == "clump_snps":
        return await clump_snps(
            sumstats_path=arguments["sumstats_path"],
            plink_prefix=arguments.get("plink_prefix"),
            p1=arguments.get("p1", 5e-8),
            p2=arguments.get("p2", 1e-4),
            r2=arguments.get("r2", 0.5),
            kb=arguments.get("kb", 250),
            output_path=arguments.get("output_path")
        )
    
    elif name == "extract_locus":
        return await extract_locus(
            plink_prefix=arguments["plink_prefix"],
            lead_snp=arguments["lead_snp"],
            window_kb=arguments.get("window_kb", 500),
            output_prefix=arguments.get("output_prefix")
        )
    
    raise ValueError(f"Unknown post-GWAS tool: {name}")


async def calculate_prs(
    sumstats_path: str,
    target_plink_prefix: str,
    pvalue_threshold: float = 1.0,
    clump_r2: float = 0.1,
    output_path: Optional[str] = None
) -> str:
    """
    Calculate Polygenic Risk Scores.
    
    Uses a simple weighted sum approach:
    PRS = sum(beta_i * dosage_i) for each SNP i
    """
    logger.info(f"Calculating PRS from: {sumstats_path}")
    
    try:
        # Read summary statistics
        sumstats = read_summary_stats(sumstats_path)
        
        # Filter by p-value
        sumstats = sumstats[sumstats['P'] <= pvalue_threshold]
        logger.info(f"SNPs after p-value filter: {len(sumstats)}")
        
        if len(sumstats) == 0:
            return json.dumps({
                "error": "No SNPs passed p-value threshold",
                "pvalue_threshold": pvalue_threshold
            }, indent=2)
        
        # Read target genotypes
        plink_data = read_plink(target_plink_prefix)
        genotypes = plink_data['genotypes']
        bim = plink_data['bim']
        fam = plink_data['fam']
        
        geno_matrix = genotypes.compute() if hasattr(genotypes, 'compute') else np.array(genotypes)
        
        # Match SNPs between sumstats and genotype data
        if 'snp' in bim.columns and 'SNP' in sumstats.columns:
            common_snps = set(bim['snp']) & set(sumstats['SNP'])
            logger.info(f"Common SNPs: {len(common_snps)}")
            
            if len(common_snps) == 0:
                return json.dumps({
                    "error": "No common SNPs between summary stats and target genotypes",
                    "sumstats_snps_sample": sumstats['SNP'].head(10).tolist(),
                    "genotype_snps_sample": bim['snp'].head(10).tolist()
                }, indent=2)
            
            # Filter to common SNPs
            sumstats_matched = sumstats[sumstats['SNP'].isin(common_snps)].copy()
            bim_matched = bim[bim['snp'].isin(common_snps)].copy()
            
            # Create SNP to index mapping
            snp_to_idx = {snp: idx for idx, snp in enumerate(bim['snp'])}
            
            # Calculate PRS for each individual
            n_samples = len(fam)
            prs_scores = np.zeros(n_samples)
            n_snps_used = 0
            
            for _, row in sumstats_matched.iterrows():
                snp = row['SNP']
                beta = row.get('BETA', 0)
                
                if snp in snp_to_idx and beta != 0:
                    idx = snp_to_idx[snp]
                    dosages = geno_matrix[:, idx]
                    
                    # Handle missing values
                    valid = ~np.isnan(dosages)
                    prs_scores[valid] += beta * dosages[valid]
                    n_snps_used += 1
            
            # Standardize PRS (Z-score)
            prs_mean = np.mean(prs_scores)
            prs_std = np.std(prs_scores)
            prs_standardized = (prs_scores - prs_mean) / prs_std if prs_std > 0 else prs_scores
            
            # Create results DataFrame
            prs_df = pd.DataFrame({
                'FID': fam['fid'].values if 'fid' in fam.columns else range(n_samples),
                'IID': fam['iid'].values if 'iid' in fam.columns else range(n_samples),
                'PRS_raw': prs_scores,
                'PRS_standardized': prs_standardized
            })
            
            result = {
                "n_individuals": n_samples,
                "n_snps_in_sumstats": len(sumstats),
                "n_snps_used": n_snps_used,
                "pvalue_threshold": pvalue_threshold,
                "prs_statistics": {
                    "mean": float(prs_mean),
                    "std": float(prs_std),
                    "min": float(prs_scores.min()),
                    "max": float(prs_scores.max())
                },
                "prs_preview": prs_df.head(10).to_dict(orient='records')
            }
            
            if output_path:
                write_results(prs_df, output_path, format='tsv')
                result["results_saved_to"] = output_path
            
            return json.dumps(result, indent=2)
        else:
            return json.dumps({
                "error": "SNP column not found in data",
                "bim_columns": list(bim.columns),
                "sumstats_columns": list(sumstats.columns)
            }, indent=2)
            
    except ImportError as e:
        return json.dumps({
            "error": "pandas-plink not installed",
            "message": str(e)
        }, indent=2)
    except Exception as e:
        logger.exception(f"PRS calculation failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


async def clump_snps(
    sumstats_path: str,
    plink_prefix: Optional[str] = None,
    p1: float = 5e-8,
    p2: float = 1e-4,
    r2: float = 0.5,
    kb: int = 250,
    output_path: Optional[str] = None
) -> str:
    """
    Perform LD-based clumping of GWAS results.
    
    Identifies independent signals by grouping SNPs in LD.
    """
    logger.info(f"Clumping SNPs from: {sumstats_path}")
    
    try:
        # Read summary statistics
        sumstats = read_summary_stats(sumstats_path)
        
        # Filter to p1 threshold
        index_snps = sumstats[sumstats['P'] <= p1].copy()
        
        if len(index_snps) == 0:
            return json.dumps({
                "n_index_snps": 0,
                "p1_threshold": p1,
                "message": "No SNPs passed p1 threshold"
            }, indent=2)
        
        # Sort by p-value
        index_snps = index_snps.sort_values('P')
        
        # If no LD reference, just return index SNPs
        if not plink_prefix:
            result = {
                "note": "No LD reference provided - returning all SNPs passing p1 threshold",
                "n_index_snps": len(index_snps),
                "clumped_snps": index_snps.to_dict(orient='records')
            }
            
            if output_path:
                write_results(index_snps, output_path, format='tsv')
                result["results_saved_to"] = output_path
            
            return json.dumps(result, indent=2, default=str)
        
        # With LD reference - perform actual clumping
        plink_data = read_plink(plink_prefix)
        genotypes = plink_data['genotypes']
        bim = plink_data['bim']
        
        geno_matrix = genotypes.compute() if hasattr(genotypes, 'compute') else np.array(genotypes)
        
        # Create position lookup
        if 'snp' in bim.columns:
            snp_to_idx = {snp: idx for idx, snp in enumerate(bim['snp'])}
        else:
            snp_to_idx = {}
        
        # Clumping algorithm
        clumped = []
        used_snps = set()
        
        for _, index_row in index_snps.iterrows():
            snp = index_row.get('SNP')
            
            if snp in used_snps:
                continue
            
            # This is an index SNP
            clumped.append(index_row.to_dict())
            used_snps.add(snp)
            
            # Find SNPs in LD within window
            if snp in snp_to_idx:
                idx = snp_to_idx[snp]
                chrom = index_row.get('CHR')
                pos = index_row.get('BP')
                
                # Find nearby SNPs
                for _, nearby_row in sumstats.iterrows():
                    nearby_snp = nearby_row.get('SNP')
                    
                    if nearby_snp in used_snps:
                        continue
                    
                    if nearby_row.get('CHR') != chrom:
                        continue
                    
                    nearby_pos = nearby_row.get('BP')
                    if abs(nearby_pos - pos) > kb * 1000:
                        continue
                    
                    # Check LD
                    if nearby_snp in snp_to_idx and nearby_row['P'] <= p2:
                        nearby_idx = snp_to_idx[nearby_snp]
                        
                        geno1 = geno_matrix[:, idx]
                        geno2 = geno_matrix[:, nearby_idx]
                        
                        valid = ~(np.isnan(geno1) | np.isnan(geno2))
                        if np.sum(valid) > 10:
                            r = np.corrcoef(geno1[valid], geno2[valid])[0, 1]
                            r2_calc = r ** 2 if not np.isnan(r) else 0
                            
                            if r2_calc >= r2:
                                used_snps.add(nearby_snp)
        
        clumped_df = pd.DataFrame(clumped)
        
        result = {
            "n_input_snps": len(sumstats),
            "n_index_snps_p1": len(index_snps),
            "n_independent_signals": len(clumped),
            "parameters": {
                "p1": p1,
                "p2": p2,
                "r2": r2,
                "kb": kb
            },
            "clumped_snps": clumped_df.to_dict(orient='records')
        }
        
        if output_path:
            write_results(clumped_df, output_path, format='tsv')
            result["results_saved_to"] = output_path
        
        return json.dumps(result, indent=2, default=str)
        
    except Exception as e:
        logger.exception(f"Clumping failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


async def extract_locus(
    plink_prefix: str,
    lead_snp: str,
    window_kb: int = 500,
    output_prefix: Optional[str] = None
) -> str:
    """
    Extract genotype data for a specific locus.
    """
    logger.info(f"Extracting locus around: {lead_snp}")
    
    try:
        # Read PLINK files
        plink_data = read_plink(plink_prefix)
        bim = plink_data['bim']
        fam = plink_data['fam']
        genotypes = plink_data['genotypes']
        
        geno_matrix = genotypes.compute() if hasattr(genotypes, 'compute') else np.array(genotypes)
        
        # Find lead SNP
        if lead_snp.startswith('rs') and 'snp' in bim.columns:
            lead_row = bim[bim['snp'] == lead_snp]
        else:
            # Try CHR:POS format
            parts = lead_snp.replace('chr', '').split(':')
            if len(parts) == 2:
                lead_row = bim[(bim['chrom'].astype(str) == parts[0]) & 
                              (bim['pos'] == int(parts[1]))]
            else:
                lead_row = pd.DataFrame()
        
        if len(lead_row) == 0:
            return json.dumps({
                "error": f"Lead SNP not found: {lead_snp}",
                "snps_sample": bim['snp'].head(10).tolist() if 'snp' in bim.columns else []
            }, indent=2)
        
        lead_row = lead_row.iloc[0]
        chrom = lead_row['chrom']
        pos = lead_row['pos']
        
        # Extract variants in window
        window = window_kb * 1000
        region_mask = (bim['chrom'] == chrom) & \
                     (bim['pos'] >= pos - window) & \
                     (bim['pos'] <= pos + window)
        
        region_bim = bim[region_mask].copy()
        region_indices = np.where(region_mask)[0]
        region_geno = geno_matrix[:, region_indices]
        
        result = {
            "lead_snp": lead_snp,
            "chromosome": str(chrom),
            "position": int(pos),
            "window_kb": window_kb,
            "n_samples": len(fam),
            "n_variants_extracted": len(region_bim),
            "region": f"chr{chrom}:{int(pos-window)}-{int(pos+window)}",
            "variants_preview": region_bim.head(20).to_dict(orient='records')
        }
        
        if output_prefix:
            # Save region BIM
            bim_path = f"{output_prefix}.bim"
            region_bim.to_csv(bim_path, sep='\t', index=False, header=False)
            
            # Save genotypes as TSV
            geno_path = f"{output_prefix}.geno.tsv"
            geno_df = pd.DataFrame(
                region_geno,
                columns=region_bim['snp'].values if 'snp' in region_bim.columns else range(len(region_bim))
            )
            geno_df.insert(0, 'IID', fam['iid'].values if 'iid' in fam.columns else range(len(fam)))
            geno_df.insert(0, 'FID', fam['fid'].values if 'fid' in fam.columns else range(len(fam)))
            write_results(geno_df, geno_path, format='tsv')
            
            result["bim_saved_to"] = bim_path
            result["genotypes_saved_to"] = geno_path
        
        return json.dumps(result, indent=2)
        
    except Exception as e:
        logger.exception(f"Locus extraction failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


def register_post_gwas_tools():
    """Return post-GWAS tools for registration."""
    return POST_GWAS_TOOLS
