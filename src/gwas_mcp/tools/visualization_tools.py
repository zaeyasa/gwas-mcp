"""
Visualization Tools for GWAS MCP Server.

Tools for creating Manhattan plots, QQ plots, regional association plots, and PCA plots.
"""

import json
import logging
import base64
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, List, Optional
import numpy as np
import pandas as pd

from mcp.types import Tool

from ..utils.validators import validate_file_path
from ..utils.file_handlers import read_summary_stats

logger = logging.getLogger(__name__)


# Define Visualization tools
VISUALIZATION_TOOLS = [
    Tool(
        name="create_manhattan_plot",
        description="Create a Manhattan plot from GWAS summary statistics. Shows -log10(p-value) across chromosomes.",
        inputSchema={
            "type": "object",
            "properties": {
                "sumstats_path": {
                    "type": "string",
                    "description": "Path to GWAS summary statistics file"
                },
                "significance_threshold": {
                    "type": "number",
                    "description": "Genome-wide significance threshold (default: 5e-8)",
                    "default": 5e-8
                },
                "suggestive_threshold": {
                    "type": "number",
                    "description": "Suggestive significance threshold (default: 1e-5)",
                    "default": 1e-5
                },
                "output_path": {
                    "type": "string",
                    "description": "Path to save the plot (PNG format)"
                },
                "title": {
                    "type": "string",
                    "description": "Plot title",
                    "default": "Manhattan Plot"
                }
            },
            "required": ["sumstats_path"]
        }
    ),
    Tool(
        name="create_qq_plot",
        description="Create a Quantile-Quantile (QQ) plot from GWAS p-values. Shows observed vs expected p-values.",
        inputSchema={
            "type": "object",
            "properties": {
                "sumstats_path": {
                    "type": "string",
                    "description": "Path to GWAS summary statistics file"
                },
                "output_path": {
                    "type": "string",
                    "description": "Path to save the plot"
                },
                "title": {
                    "type": "string",
                    "description": "Plot title",
                    "default": "QQ Plot"
                }
            },
            "required": ["sumstats_path"]
        }
    ),
    Tool(
        name="create_regional_plot",
        description="Create a regional association plot (LocusZoom-style) for a specific genomic region around a lead SNP.",
        inputSchema={
            "type": "object",
            "properties": {
                "sumstats_path": {
                    "type": "string",
                    "description": "Path to GWAS summary statistics"
                },
                "lead_snp": {
                    "type": "string",
                    "description": "Lead SNP rsID or CHR:POS"
                },
                "window_kb": {
                    "type": "integer",
                    "description": "Window size in kilobases (default: 500)",
                    "default": 500
                },
                "output_path": {
                    "type": "string",
                    "description": "Path to save the plot"
                }
            },
            "required": ["sumstats_path", "lead_snp"]
        }
    ),
    Tool(
        name="create_pca_plot",
        description="Create a PCA scatter plot showing population structure. Plots PC1 vs PC2.",
        inputSchema={
            "type": "object",
            "properties": {
                "pca_results_path": {
                    "type": "string",
                    "description": "Path to PCA results file (with PC1, PC2 columns)"
                },
                "population_labels_path": {
                    "type": "string",
                    "description": "Optional path to population labels file"
                },
                "output_path": {
                    "type": "string",
                    "description": "Path to save the plot"
                },
                "title": {
                    "type": "string",
                    "description": "Plot title",
                    "default": "PCA Plot"
                }
            },
            "required": ["pca_results_path"]
        }
    ),
]


async def handle_visualization_tool(name: str, arguments: Dict[str, Any]) -> str:
    """Handle visualization tool calls."""
    
    if name == "create_manhattan_plot":
        return await create_manhattan_plot(
            sumstats_path=arguments["sumstats_path"],
            significance_threshold=arguments.get("significance_threshold", 5e-8),
            suggestive_threshold=arguments.get("suggestive_threshold", 1e-5),
            output_path=arguments.get("output_path"),
            title=arguments.get("title", "Manhattan Plot")
        )
    
    elif name == "create_qq_plot":
        return await create_qq_plot(
            sumstats_path=arguments["sumstats_path"],
            output_path=arguments.get("output_path"),
            title=arguments.get("title", "QQ Plot")
        )
    
    elif name == "create_regional_plot":
        return await create_regional_plot(
            sumstats_path=arguments["sumstats_path"],
            lead_snp=arguments["lead_snp"],
            window_kb=arguments.get("window_kb", 500),
            output_path=arguments.get("output_path")
        )
    
    elif name == "create_pca_plot":
        return await create_pca_plot(
            pca_results_path=arguments["pca_results_path"],
            population_labels_path=arguments.get("population_labels_path"),
            output_path=arguments.get("output_path"),
            title=arguments.get("title", "PCA Plot")
        )
    
    raise ValueError(f"Unknown visualization tool: {name}")


async def create_manhattan_plot(
    sumstats_path: str,
    significance_threshold: float = 5e-8,
    suggestive_threshold: float = 1e-5,
    output_path: Optional[str] = None,
    title: str = "Manhattan Plot"
) -> str:
    """
    Create a Manhattan plot from GWAS summary statistics.
    """
    logger.info(f"Creating Manhattan plot from: {sumstats_path}")
    
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        
        # Read summary statistics
        df = read_summary_stats(sumstats_path)
        
        # Ensure required columns
        if 'CHR' not in df.columns or 'BP' not in df.columns or 'P' not in df.columns:
            raise ValueError("Summary stats must have CHR, BP, and P columns")
        
        # Clean data
        df = df.dropna(subset=['CHR', 'BP', 'P'])
        df = df[df['P'] > 0]  # Remove zero p-values
        
        # Convert chromosome to numeric
        df['CHR'] = pd.to_numeric(df['CHR'].astype(str).str.replace('chr', '', case=False), errors='coerce')
        df = df.dropna(subset=['CHR'])
        df['CHR'] = df['CHR'].astype(int)
        
        # Sort by chromosome and position
        df = df.sort_values(['CHR', 'BP'])
        
        # Calculate cumulative position
        df['cumulative_pos'] = 0
        offset = 0
        chr_centers = {}
        
        for chrom in range(1, 23):
            chr_data = df[df['CHR'] == chrom]
            if len(chr_data) > 0:
                df.loc[df['CHR'] == chrom, 'cumulative_pos'] = chr_data['BP'] + offset
                chr_centers[chrom] = offset + chr_data['BP'].median()
                offset = df[df['CHR'] == chrom]['cumulative_pos'].max() + 1e6
        
        # Calculate -log10(p)
        df['neg_log_p'] = -np.log10(df['P'])
        
        # Create plot
        fig, ax = plt.subplots(figsize=(14, 6))
        
        # Color by chromosome
        colors = ['#1f77b4', '#aec7e8']  # Alternating colors
        
        for chrom in range(1, 23):
            chr_data = df[df['CHR'] == chrom]
            color = colors[chrom % 2]
            ax.scatter(chr_data['cumulative_pos'], chr_data['neg_log_p'],
                      c=color, s=3, alpha=0.6)
        
        # Add significance lines
        ax.axhline(y=-np.log10(significance_threshold), color='red', 
                  linestyle='--', linewidth=1, label=f'GW sig (p={significance_threshold})')
        ax.axhline(y=-np.log10(suggestive_threshold), color='blue',
                  linestyle='--', linewidth=0.5, label=f'Suggestive (p={suggestive_threshold})')
        
        # Labels
        ax.set_xlabel('Chromosome')
        ax.set_ylabel('-log₁₀(p-value)')
        ax.set_title(title)
        
        # Set x-axis ticks
        ax.set_xticks(list(chr_centers.values()))
        ax.set_xticklabels(list(chr_centers.keys()))
        
        ax.legend(loc='upper right')
        ax.set_xlim([0, df['cumulative_pos'].max()])
        ax.set_ylim([0, max(df['neg_log_p'].max() * 1.1, -np.log10(significance_threshold) + 2)])
        
        plt.tight_layout()
        
        # Save or encode
        result = {
            "title": title,
            "n_variants_plotted": len(df),
            "n_significant": len(df[df['P'] < significance_threshold]),
            "n_suggestive": len(df[(df['P'] >= significance_threshold) & (df['P'] < suggestive_threshold)]),
            "max_neg_log_p": float(df['neg_log_p'].max())
        }
        
        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            result["plot_saved_to"] = output_path
        else:
            # Return as base64
            buffer = BytesIO()
            plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
            buffer.seek(0)
            result["plot_base64"] = base64.b64encode(buffer.getvalue()).decode('utf-8')
        
        plt.close()
        
        return json.dumps(result, indent=2)
        
    except ImportError as e:
        return json.dumps({
            "error": "matplotlib not installed",
            "message": str(e),
            "suggestion": "Install with: pip install matplotlib seaborn"
        }, indent=2)
    except Exception as e:
        logger.exception(f"Manhattan plot failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


async def create_qq_plot(
    sumstats_path: str,
    output_path: Optional[str] = None,
    title: str = "QQ Plot"
) -> str:
    """
    Create a QQ plot from GWAS p-values.
    """
    logger.info(f"Creating QQ plot from: {sumstats_path}")
    
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')
        from scipy import stats
        
        # Read summary statistics
        df = read_summary_stats(sumstats_path)
        
        # Get p-values
        pvalues = df['P'].dropna().values
        pvalues = pvalues[(pvalues > 0) & (pvalues <= 1)]
        
        # Sort p-values
        observed = np.sort(pvalues)
        n = len(observed)
        
        # Expected p-values under null hypothesis
        expected = np.arange(1, n + 1) / (n + 1)
        
        # Transform to -log10
        observed_log = -np.log10(observed)
        expected_log = -np.log10(expected)
        
        # Calculate lambda GC
        chi2_stats = stats.chi2.ppf(1 - pvalues, df=1)
        lambda_gc = np.median(chi2_stats) / stats.chi2.ppf(0.5, df=1)
        
        # Create plot
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Plot points
        ax.scatter(expected_log, observed_log, c='#1f77b4', s=3, alpha=0.5)
        
        # Plot diagonal line
        max_val = max(expected_log.max(), observed_log.max())
        ax.plot([0, max_val], [0, max_val], 'r--', linewidth=1)
        
        # Labels
        ax.set_xlabel('Expected -log₁₀(p)')
        ax.set_ylabel('Observed -log₁₀(p)')
        ax.set_title(f'{title}\nλGC = {lambda_gc:.3f}')
        
        ax.set_xlim([0, max_val * 1.05])
        ax.set_ylim([0, max_val * 1.05])
        
        plt.tight_layout()
        
        result = {
            "title": title,
            "n_variants": len(pvalues),
            "lambda_gc": round(lambda_gc, 4),
            "max_observed_neg_log_p": float(observed_log.max())
        }
        
        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            result["plot_saved_to"] = output_path
        else:
            buffer = BytesIO()
            plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
            buffer.seek(0)
            result["plot_base64"] = base64.b64encode(buffer.getvalue()).decode('utf-8')
        
        plt.close()
        
        return json.dumps(result, indent=2)
        
    except ImportError as e:
        return json.dumps({
            "error": "Required packages not installed",
            "message": str(e)
        }, indent=2)
    except Exception as e:
        logger.exception(f"QQ plot failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


async def create_regional_plot(
    sumstats_path: str,
    lead_snp: str,
    window_kb: int = 500,
    output_path: Optional[str] = None
) -> str:
    """
    Create a regional association plot around a lead SNP.
    """
    logger.info(f"Creating regional plot for: {lead_snp}")
    
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')
        
        # Read summary statistics
        df = read_summary_stats(sumstats_path)
        
        # Find lead SNP
        if lead_snp.startswith('rs'):
            lead_row = df[df['SNP'] == lead_snp]
        else:
            # CHR:POS format
            parts = lead_snp.replace('chr', '').split(':')
            lead_row = df[(df['CHR'].astype(str) == parts[0]) & (df['BP'] == int(parts[1]))]
        
        if len(lead_row) == 0:
            return json.dumps({
                "error": f"Lead SNP not found: {lead_snp}",
                "available_snps_sample": df['SNP'].head(10).tolist() if 'SNP' in df.columns else []
            }, indent=2)
        
        lead_row = lead_row.iloc[0]
        chrom = lead_row['CHR']
        pos = lead_row['BP']
        
        # Filter to region
        window = window_kb * 1000
        region_df = df[(df['CHR'] == chrom) & 
                       (df['BP'] >= pos - window) & 
                       (df['BP'] <= pos + window)].copy()
        
        if len(region_df) == 0:
            return json.dumps({
                "error": "No variants found in region"
            }, indent=2)
        
        region_df['neg_log_p'] = -np.log10(region_df['P'])
        
        # Create plot
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Plot all SNPs
        ax.scatter(region_df['BP'] / 1e6, region_df['neg_log_p'],
                  c='#1f77b4', s=20, alpha=0.7)
        
        # Highlight lead SNP
        ax.scatter([pos / 1e6], [-np.log10(lead_row['P'])],
                  c='red', s=100, marker='D', zorder=10, label='Lead SNP')
        
        # Add significance line
        ax.axhline(y=-np.log10(5e-8), color='red', linestyle='--', linewidth=1)
        
        ax.set_xlabel(f'Position on chromosome {chrom} (Mb)')
        ax.set_ylabel('-log₁₀(p-value)')
        ax.set_title(f'Regional Plot: {lead_snp}')
        ax.legend()
        
        plt.tight_layout()
        
        result = {
            "lead_snp": lead_snp,
            "chromosome": str(chrom),
            "position": int(pos),
            "lead_pvalue": float(lead_row['P']),
            "n_variants_in_region": len(region_df),
            "region": f"chr{chrom}:{pos-window}-{pos+window}"
        }
        
        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            result["plot_saved_to"] = output_path
        else:
            buffer = BytesIO()
            plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
            buffer.seek(0)
            result["plot_base64"] = base64.b64encode(buffer.getvalue()).decode('utf-8')
        
        plt.close()
        
        return json.dumps(result, indent=2)
        
    except Exception as e:
        logger.exception(f"Regional plot failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


async def create_pca_plot(
    pca_results_path: str,
    population_labels_path: Optional[str] = None,
    output_path: Optional[str] = None,
    title: str = "PCA Plot"
) -> str:
    """
    Create a PCA scatter plot.
    """
    logger.info(f"Creating PCA plot from: {pca_results_path}")
    
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')
        
        # Read PCA results
        pca_df = pd.read_csv(pca_results_path, sep='\t')
        
        # Find PC columns
        pc1_col = 'PC1' if 'PC1' in pca_df.columns else pca_df.columns[2]
        pc2_col = 'PC2' if 'PC2' in pca_df.columns else pca_df.columns[3]
        
        # Read population labels if provided
        if population_labels_path:
            pop_df = pd.read_csv(population_labels_path, sep='\t')
            pca_df = pca_df.merge(pop_df, on=['FID', 'IID'], how='left')
            pop_col = [c for c in pop_df.columns if c not in ['FID', 'IID']][0]
        else:
            pop_col = None
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        if pop_col and pop_col in pca_df.columns:
            # Color by population
            populations = pca_df[pop_col].unique()
            colors = plt.cm.tab10(np.linspace(0, 1, len(populations)))
            
            for pop, color in zip(populations, colors):
                pop_data = pca_df[pca_df[pop_col] == pop]
                ax.scatter(pop_data[pc1_col], pop_data[pc2_col],
                          c=[color], s=20, alpha=0.7, label=pop)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            ax.scatter(pca_df[pc1_col], pca_df[pc2_col],
                      c='#1f77b4', s=20, alpha=0.7)
        
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_title(title)
        
        plt.tight_layout()
        
        result = {
            "title": title,
            "n_samples": len(pca_df),
            "n_populations": len(pca_df[pop_col].unique()) if pop_col else 1
        }
        
        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            result["plot_saved_to"] = output_path
        else:
            buffer = BytesIO()
            plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
            buffer.seek(0)
            result["plot_base64"] = base64.b64encode(buffer.getvalue()).decode('utf-8')
        
        plt.close()
        
        return json.dumps(result, indent=2)
        
    except Exception as e:
        logger.exception(f"PCA plot failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


def register_visualization_tools():
    """Return visualization tools for registration."""
    return VISUALIZATION_TOOLS
