"""
Input validation utilities for GWAS MCP Server.
"""

import os
import re
import logging
from pathlib import Path
from typing import Optional, Tuple, List

logger = logging.getLogger(__name__)


def validate_file_path(
    file_path: str,
    must_exist: bool = True,
    allowed_extensions: Optional[List[str]] = None,
    base_dir: Optional[str] = None
) -> Path:
    """
    Validate and sanitize a file path.
    
    Args:
        file_path: Path to validate
        must_exist: Whether the file must exist
        allowed_extensions: List of allowed file extensions
        base_dir: Optional base directory to restrict access
    
    Returns:
        Validated Path object
    
    Raises:
        ValueError: If path is invalid
        FileNotFoundError: If file doesn't exist and must_exist=True
    """
    path = Path(file_path).resolve()
    
    # Prevent directory traversal attacks
    if base_dir:
        base = Path(base_dir).resolve()
        if not str(path).startswith(str(base)):
            raise ValueError(f"Access denied: path outside allowed directory")
    
    # Check for suspicious patterns
    suspicious_patterns = ['..', '~', '$', '|', ';', '&']
    for pattern in suspicious_patterns:
        if pattern in str(file_path):
            raise ValueError(f"Invalid path: contains suspicious pattern '{pattern}'")
    
    # Check extension
    if allowed_extensions:
        ext = path.suffix.lower().lstrip('.')
        # Handle double extensions like .vcf.gz
        double_ext = ''.join(path.suffixes).lower().lstrip('.')
        
        if ext not in allowed_extensions and double_ext not in allowed_extensions:
            raise ValueError(
                f"Invalid file extension: {path.suffix}. "
                f"Allowed: {allowed_extensions}"
            )
    
    # Check existence
    if must_exist and not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    
    return path


def validate_rsid(rsid: str) -> str:
    """
    Validate an rsID (SNP identifier).
    
    Args:
        rsid: rsID to validate (e.g., 'rs12345')
    
    Returns:
        Normalized rsID
    
    Raises:
        ValueError: If rsID format is invalid
    """
    rsid = rsid.strip().lower()
    
    # Standard rsID format
    if re.match(r'^rs\d+$', rsid):
        return rsid
    
    # Some databases use RS prefix
    if re.match(r'^RS\d+$', rsid, re.IGNORECASE):
        return rsid.lower()
    
    raise ValueError(
        f"Invalid rsID format: {rsid}. "
        "Expected format: rs followed by numbers (e.g., rs12345)"
    )


def validate_genomic_region(region: str) -> Tuple[str, int, int]:
    """
    Validate and parse a genomic region string.
    
    Args:
        region: Region string (e.g., 'chr1:1000000-2000000' or '1:1000000-2000000')
    
    Returns:
        Tuple of (chromosome, start, end)
    
    Raises:
        ValueError: If region format is invalid
    """
    # Pattern: chr1:1000-2000 or 1:1000-2000
    pattern = r'^(chr)?(\d{1,2}|X|Y|MT?):(\d+)-(\d+)$'
    match = re.match(pattern, region, re.IGNORECASE)
    
    if not match:
        raise ValueError(
            f"Invalid genomic region: {region}. "
            "Expected format: chr1:1000000-2000000 or 1:1000000-2000000"
        )
    
    chrom = match.group(2).upper()
    if chrom == 'M':
        chrom = 'MT'
    
    start = int(match.group(3))
    end = int(match.group(4))
    
    if start >= end:
        raise ValueError(f"Invalid region: start ({start}) must be less than end ({end})")
    
    if start < 0:
        raise ValueError(f"Invalid region: start position must be non-negative")
    
    return (chrom, start, end)


def validate_pvalue(
    pvalue: float,
    allow_zero: bool = False,
    min_value: float = 0.0,
    max_value: float = 1.0
) -> float:
    """
    Validate a p-value.
    
    Args:
        pvalue: P-value to validate
        allow_zero: Whether to allow exactly 0
        min_value: Minimum allowed value
        max_value: Maximum allowed value
    
    Returns:
        Validated p-value
    
    Raises:
        ValueError: If p-value is invalid
    """
    try:
        pvalue = float(pvalue)
    except (TypeError, ValueError):
        raise ValueError(f"Invalid p-value: {pvalue}. Must be a number.")
    
    if pvalue < 0:
        raise ValueError(f"Invalid p-value: {pvalue}. Cannot be negative.")
    
    if pvalue > max_value:
        raise ValueError(f"Invalid p-value: {pvalue}. Cannot exceed {max_value}.")
    
    if not allow_zero and pvalue == 0:
        raise ValueError("P-value cannot be exactly 0. Use a very small value instead.")
    
    return pvalue


def validate_threshold(
    value: float,
    name: str,
    min_value: Optional[float] = None,
    max_value: Optional[float] = None
) -> float:
    """
    Validate a numeric threshold parameter.
    
    Args:
        value: Value to validate
        name: Parameter name (for error messages)
        min_value: Minimum allowed value
        max_value: Maximum allowed value
    
    Returns:
        Validated value
    """
    try:
        value = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"Invalid {name}: {value}. Must be a number.")
    
    if min_value is not None and value < min_value:
        raise ValueError(f"Invalid {name}: {value}. Must be >= {min_value}.")
    
    if max_value is not None and value > max_value:
        raise ValueError(f"Invalid {name}: {value}. Must be <= {max_value}.")
    
    return value


def validate_chromosome(chrom: str) -> str:
    """
    Validate and normalize a chromosome name.
    
    Args:
        chrom: Chromosome name (e.g., '1', 'chr1', 'X', 'chrX')
    
    Returns:
        Normalized chromosome name
    """
    chrom = str(chrom).strip().upper()
    
    # Remove 'CHR' prefix if present
    if chrom.startswith('CHR'):
        chrom = chrom[3:]
    
    # Valid chromosomes
    valid_chroms = set(str(i) for i in range(1, 23)) | {'X', 'Y', 'MT', 'M'}
    
    if chrom == 'M':
        chrom = 'MT'
    
    if chrom not in valid_chroms:
        raise ValueError(
            f"Invalid chromosome: {chrom}. "
            f"Valid values: 1-22, X, Y, MT"
        )
    
    return chrom
