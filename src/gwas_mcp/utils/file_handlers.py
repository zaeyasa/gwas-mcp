"""
File handlers for VCF, PLINK, and summary statistics files.
"""

import os
import json
import logging
from pathlib import Path
from typing import Optional, Dict, Any, Union
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def read_vcf(vcf_path: str, region: Optional[str] = None) -> pd.DataFrame:
    """
    Read a VCF file and return as DataFrame.
    
    Args:
        vcf_path: Path to VCF file (.vcf or .vcf.gz)
        region: Optional genomic region (chr:start-end)
    
    Returns:
        DataFrame with variant information
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    
    try:
        # Try using pysam for efficient VCF parsing
        import pysam
        
        vcf = pysam.VariantFile(str(vcf_path))
        records = []
        
        if region:
            iterator = vcf.fetch(region=region)
        else:
            iterator = vcf.fetch()
        
        for record in iterator:
            rec_dict = {
                'CHROM': record.chrom,
                'POS': record.pos,
                'ID': record.id or '.',
                'REF': record.ref,
                'ALT': ','.join(str(a) for a in record.alts) if record.alts else '.',
                'QUAL': record.qual,
                'FILTER': ','.join(record.filter.keys()) if record.filter else 'PASS',
            }
            # Add INFO fields
            for key in record.info.keys():
                rec_dict[f'INFO_{key}'] = record.info[key]
            records.append(rec_dict)
        
        vcf.close()
        return pd.DataFrame(records)
        
    except ImportError:
        logger.warning("pysam not available, using fallback VCF parser")
        return _read_vcf_fallback(vcf_path)


def _read_vcf_fallback(vcf_path: Path) -> pd.DataFrame:
    """Fallback VCF parser using pandas (less efficient)."""
    import gzip
    
    opener = gzip.open if str(vcf_path).endswith('.gz') else open
    
    records = []
    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    headers = line.strip().split('\t')
                continue
            
            fields = line.strip().split('\t')
            rec = {
                'CHROM': fields[0],
                'POS': int(fields[1]),
                'ID': fields[2],
                'REF': fields[3],
                'ALT': fields[4],
                'QUAL': fields[5],
                'FILTER': fields[6],
                'INFO': fields[7] if len(fields) > 7 else '',
            }
            records.append(rec)
    
    return pd.DataFrame(records)


def read_plink(prefix: str) -> Dict[str, Any]:
    """
    Read PLINK binary files (bed/bim/fam).
    
    Args:
        prefix: Path prefix for PLINK files (without .bed/.bim/.fam extension)
    
    Returns:
        Dictionary with 'bim', 'fam', and 'genotypes' DataFrames
    """
    try:
        from pandas_plink import read_plink as pandas_read_plink
        
        (bim, fam, genotypes) = pandas_read_plink(prefix)
        
        return {
            'bim': bim,
            'fam': fam,
            'genotypes': genotypes,
            'n_samples': len(fam),
            'n_variants': len(bim),
        }
        
    except ImportError:
        logger.error("pandas-plink not installed. Please install: pip install pandas-plink")
        raise ImportError("pandas-plink required for reading PLINK files")
    except Exception as e:
        logger.error(f"Error reading PLINK files: {e}")
        raise


def read_summary_stats(
    file_path: str,
    sep: str = '\t',
    column_mapping: Optional[Dict[str, str]] = None
) -> pd.DataFrame:
    """
    Read GWAS summary statistics file.
    
    Args:
        file_path: Path to summary statistics file
        sep: Field separator (default: tab)
        column_mapping: Optional mapping of file columns to standard names
    
    Returns:
        DataFrame with standardized column names
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Summary stats file not found: {file_path}")
    
    # Detect compression
    if str(file_path).endswith('.gz'):
        df = pd.read_csv(file_path, sep=sep, compression='gzip')
    else:
        df = pd.read_csv(file_path, sep=sep)
    
    # Standard column names
    standard_columns = {
        'chr': 'CHR',
        'chromosome': 'CHR',
        'chrom': 'CHR',
        'pos': 'BP',
        'position': 'BP',
        'bp': 'BP',
        'snp': 'SNP',
        'rsid': 'SNP',
        'marker': 'SNP',
        'variant_id': 'SNP',
        'p': 'P',
        'pvalue': 'P',
        'p_value': 'P',
        'pval': 'P',
        'beta': 'BETA',
        'effect': 'BETA',
        'b': 'BETA',
        'se': 'SE',
        'stderr': 'SE',
        'standard_error': 'SE',
        'a1': 'A1',
        'effect_allele': 'A1',
        'alt': 'A1',
        'a2': 'A2',
        'ref': 'A2',
        'other_allele': 'A2',
        'maf': 'MAF',
        'eaf': 'MAF',
        'freq': 'MAF',
    }
    
    # Apply default column mapping
    df.columns = df.columns.str.lower()
    df = df.rename(columns=standard_columns)
    
    # Apply custom mapping if provided
    if column_mapping:
        df = df.rename(columns=column_mapping)
    
    return df


def write_results(
    data: Union[pd.DataFrame, Dict, list],
    output_path: str,
    format: str = 'tsv'
) -> str:
    """
    Write results to file.
    
    Args:
        data: Data to write (DataFrame, dict, or list)
        output_path: Output file path
        format: Output format ('tsv', 'csv', 'json')
    
    Returns:
        Path to written file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if isinstance(data, pd.DataFrame):
        if format == 'tsv':
            data.to_csv(output_path, sep='\t', index=False)
        elif format == 'csv':
            data.to_csv(output_path, index=False)
        elif format == 'json':
            data.to_json(output_path, orient='records', indent=2)
        else:
            raise ValueError(f"Unsupported format: {format}")
    else:
        # Dict or list
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2, default=str)
    
    logger.info(f"Results written to: {output_path}")
    return str(output_path)
