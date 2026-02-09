"""
Database Resources for GWAS MCP Server.

Provides access to GWAS Catalog, Ensembl, and GTEx data via MCP resources.
"""

import json
import logging
import re
from typing import Optional
from urllib.parse import urlparse, unquote

import aiohttp

from mcp.types import Resource

logger = logging.getLogger(__name__)

# API endpoints
ENSEMBL_REST_API = "https://rest.ensembl.org"
GWAS_CATALOG_API = "https://www.ebi.ac.uk/gwas/rest/api"
GTEX_API = "https://gtexportal.org/api/v2"


# Define Resources
# Only include static resources that work without user input
# Dynamic lookups (genes, proteins, SNPs) are handled via Claude's natural language
RESOURCES = [
    Resource(
        uri="gwas://catalog/traits",
        name="GWAS Catalog Traits",
        description="List of all traits/diseases in the GWAS Catalog database",
        mimeType="application/json"
    ),
]


async def handle_resource(uri: str) -> str:
    """Handle resource read requests."""
    logger.info(f"Reading resource: {uri}")
    
    parsed = urlparse(str(uri))
    scheme = parsed.scheme
    path = unquote(parsed.netloc + parsed.path)
    
    if scheme == "gwas":
        return await handle_gwas_resource(path)
    elif scheme == "ensembl":
        return await handle_ensembl_resource(path)
    elif scheme == "gtex":
        return await handle_gtex_resource(path)
    else:
        raise ValueError(f"Unknown resource scheme: {scheme}")


async def handle_gwas_resource(path: str) -> str:
    """Handle GWAS Catalog resource requests."""
    
    parts = path.strip('/').split('/')
    
    if parts[0] == "catalog":
        if len(parts) == 2 and parts[1] == "traits":
            # List all traits
            return await fetch_gwas_traits()
        elif len(parts) == 3 and parts[1] == "snp":
            # SNP lookup
            rsid = parts[2]
            return await fetch_gwas_snp(rsid)
    
    raise ValueError(f"Invalid GWAS Catalog resource path: {path}")


async def handle_ensembl_resource(path: str) -> str:
    """Handle Ensembl resource requests."""
    
    parts = path.strip('/').split('/')
    
    if len(parts) >= 2:
        if parts[0] == "gene":
            gene_name = parts[1]
            return await fetch_ensembl_gene(gene_name)
        elif parts[0] == "variant":
            rsid = parts[1]
            return await fetch_ensembl_variant(rsid)
    
    raise ValueError(f"Invalid Ensembl resource path: {path}")


async def handle_gtex_resource(path: str) -> str:
    """Handle GTEx resource requests."""
    
    parts = path.strip('/').split('/')
    
    if len(parts) >= 2 and parts[0] == "eqtl":
        snp = parts[1]
        tissue = parts[2] if len(parts) > 2 else None
        return await fetch_gtex_eqtl(snp, tissue)
    
    raise ValueError(f"Invalid GTEx resource path: {path}")


async def fetch_gwas_traits() -> str:
    """Fetch list of traits from GWAS Catalog."""
    
    try:
        async with aiohttp.ClientSession() as session:
            url = f"{GWAS_CATALOG_API}/efoTraits?size=100"
            
            async with session.get(url, headers={"Accept": "application/json"}) as response:
                if response.status == 200:
                    data = await response.json()
                    
                    traits = []
                    if "_embedded" in data and "efoTraits" in data["_embedded"]:
                        for trait in data["_embedded"]["efoTraits"]:
                            traits.append({
                                "trait": trait.get("trait"),
                                "uri": trait.get("uri"),
                                "shortForm": trait.get("shortForm")
                            })
                    
                    return json.dumps({
                        "source": "GWAS Catalog",
                        "n_traits": len(traits),
                        "traits": traits
                    }, indent=2)
                else:
                    return json.dumps({
                        "error": f"API returned status {response.status}"
                    }, indent=2)
                    
    except Exception as e:
        logger.exception(f"GWAS Catalog fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def fetch_gwas_snp(rsid: str) -> str:
    """Fetch GWAS associations for a SNP."""
    
    try:
        async with aiohttp.ClientSession() as session:
            url = f"{GWAS_CATALOG_API}/singleNucleotidePolymorphisms/{rsid}/associations"
            
            async with session.get(url, headers={"Accept": "application/json"}) as response:
                if response.status == 200:
                    data = await response.json()
                    
                    associations = []
                    if "_embedded" in data and "associations" in data["_embedded"]:
                        for assoc in data["_embedded"]["associations"][:20]:
                            associations.append({
                                "pvalue": assoc.get("pvalue"),
                                "pvalue_text": assoc.get("pvalueText"),
                                "beta": assoc.get("betaNum"),
                                "odds_ratio": assoc.get("orPerCopyNum"),
                                "risk_frequency": assoc.get("riskFrequency")
                            })
                    
                    return json.dumps({
                        "rsid": rsid,
                        "source": "GWAS Catalog",
                        "n_associations": len(associations),
                        "associations": associations
                    }, indent=2)
                elif response.status == 404:
                    return json.dumps({
                        "rsid": rsid,
                        "source": "GWAS Catalog",
                        "n_associations": 0,
                        "message": "No associations found"
                    }, indent=2)
                else:
                    return json.dumps({
                        "error": f"API returned status {response.status}"
                    }, indent=2)
                    
    except Exception as e:
        logger.exception(f"GWAS Catalog SNP fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def fetch_ensembl_gene(gene_name: str) -> str:
    """Fetch gene information from Ensembl."""
    
    try:
        async with aiohttp.ClientSession() as session:
            # Search for gene
            url = f"{ENSEMBL_REST_API}/lookup/symbol/homo_sapiens/{gene_name}"
            headers = {"Content-Type": "application/json"}
            
            async with session.get(url, headers=headers) as response:
                if response.status == 200:
                    data = await response.json()
                    
                    gene_info = {
                        "symbol": data.get("display_name"),
                        "ensembl_id": data.get("id"),
                        "description": data.get("description"),
                        "biotype": data.get("biotype"),
                        "chromosome": data.get("seq_region_name"),
                        "start": data.get("start"),
                        "end": data.get("end"),
                        "strand": data.get("strand"),
                        "source": "Ensembl"
                    }
                    
                    return json.dumps(gene_info, indent=2)
                elif response.status == 404:
                    return json.dumps({
                        "gene": gene_name,
                        "error": "Gene not found"
                    }, indent=2)
                else:
                    return json.dumps({
                        "error": f"API returned status {response.status}"
                    }, indent=2)
                    
    except Exception as e:
        logger.exception(f"Ensembl gene fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def fetch_ensembl_variant(rsid: str) -> str:
    """Fetch variant information from Ensembl."""
    
    try:
        async with aiohttp.ClientSession() as session:
            url = f"{ENSEMBL_REST_API}/variation/human/{rsid}"
            headers = {"Content-Type": "application/json"}
            
            async with session.get(url, headers=headers) as response:
                if response.status == 200:
                    data = await response.json()
                    
                    # Get mapping info
                    mapping = data.get("mappings", [{}])[0]
                    
                    variant_info = {
                        "rsid": rsid,
                        "source": "Ensembl",
                        "chromosome": mapping.get("seq_region_name"),
                        "position": mapping.get("start"),
                        "alleles": mapping.get("allele_string"),
                        "ancestral_allele": data.get("ancestral_allele"),
                        "minor_allele": data.get("minor_allele"),
                        "maf": data.get("MAF"),
                        "consequence": data.get("most_severe_consequence"),
                        "clinical_significance": data.get("clinical_significance", [])
                    }
                    
                    return json.dumps(variant_info, indent=2)
                elif response.status == 404:
                    return json.dumps({
                        "rsid": rsid,
                        "error": "Variant not found"
                    }, indent=2)
                else:
                    return json.dumps({
                        "error": f"API returned status {response.status}"
                    }, indent=2)
                    
    except Exception as e:
        logger.exception(f"Ensembl variant fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def fetch_gtex_eqtl(snp: str, tissue: Optional[str] = None) -> str:
    """Fetch eQTL data from GTEx."""
    
    try:
        async with aiohttp.ClientSession() as session:
            params = {"snpId": snp, "datasetId": "gtex_v8"}
            if tissue:
                params["tissueSiteDetailId"] = tissue
            
            url = f"{GTEX_API}/association/singleTissueEqtl"
            
            async with session.get(url, params=params) as response:
                if response.status == 200:
                    data = await response.json()
                    
                    eqtls = []
                    for eqtl in data.get("data", [])[:30]:
                        eqtls.append({
                            "gene_symbol": eqtl.get("geneSymbol"),
                            "gene_id": eqtl.get("geneId"),
                            "tissue": eqtl.get("tissueSiteDetailId"),
                            "pvalue": eqtl.get("pValue"),
                            "nes": eqtl.get("nes"),
                            "slope": eqtl.get("slope")
                        })
                    
                    return json.dumps({
                        "snp": snp,
                        "tissue_filter": tissue,
                        "source": "GTEx v8",
                        "n_eqtls": len(eqtls),
                        "eqtls": eqtls
                    }, indent=2)
                elif response.status == 404:
                    return json.dumps({
                        "snp": snp,
                        "source": "GTEx",
                        "n_eqtls": 0,
                        "message": "No eQTL data found"
                    }, indent=2)
                else:
                    return json.dumps({
                        "error": f"API returned status {response.status}"
                    }, indent=2)
                    
    except Exception as e:
        logger.exception(f"GTEx fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


def register_resources():
    """Return resources for registration."""
    return RESOURCES
