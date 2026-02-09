"""
Annotation Tools for GWAS MCP Server.

Tools for SNP annotation, GWAS Catalog queries, eQTL lookups, and gene enrichment.
"""

import json
import logging
from typing import Any, Dict, List, Optional
import asyncio

import aiohttp
import requests
import pandas as pd

from mcp.types import Tool

from ..utils.validators import validate_rsid, validate_genomic_region

logger = logging.getLogger(__name__)

# API endpoints
ENSEMBL_REST_API = "https://rest.ensembl.org"
GWAS_CATALOG_API = "https://www.ebi.ac.uk/gwas/rest/api"
GTEX_API = "https://gtexportal.org/api/v2"


# Define Annotation tools
ANNOTATION_TOOLS = [
    Tool(
        name="annotate_snps",
        description="Annotate SNPs with gene names, functional consequences, and allele frequencies from Ensembl VEP and gnomAD.",
        inputSchema={
            "type": "object",
            "properties": {
                "rsids": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of rsIDs to annotate (e.g., ['rs12345', 'rs67890'])"
                },
                "include_frequencies": {
                    "type": "boolean",
                    "description": "Include population allele frequencies (default: true)",
                    "default": True
                }
            },
            "required": ["rsids"]
        }
    ),
    Tool(
        name="query_gwas_catalog",
        description="Query the NHGRI-EBI GWAS Catalog for previously reported associations for a SNP or trait.",
        inputSchema={
            "type": "object",
            "properties": {
                "rsid": {
                    "type": "string",
                    "description": "rsID to query (e.g., 'rs12345')"
                },
                "trait": {
                    "type": "string",
                    "description": "Trait/disease to search for (e.g., 'diabetes')"
                },
                "gene": {
                    "type": "string",
                    "description": "Gene name to search for (e.g., 'BRCA1')"
                }
            }
        }
    ),
    Tool(
        name="get_eqtl_data",
        description="Get expression quantitative trait loci (eQTL) data from GTEx for a SNP and tissue.",
        inputSchema={
            "type": "object",
            "properties": {
                "rsid": {
                    "type": "string",
                    "description": "rsID to query"
                },
                "tissue": {
                    "type": "string",
                    "description": "GTEx tissue type (e.g., 'Whole_Blood', 'Brain_Cortex')"
                },
                "gene": {
                    "type": "string",
                    "description": "Optional gene to filter results"
                }
            },
            "required": ["rsid"]
        }
    ),
    Tool(
        name="gene_set_enrichment",
        description="Perform gene set enrichment analysis using GO and KEGG pathways via Enrichr.",
        inputSchema={
            "type": "object",
            "properties": {
                "genes": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of gene symbols (e.g., ['BRCA1', 'TP53', 'EGFR'])"
                },
                "databases": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Databases to query (default: GO_Biological_Process, KEGG)",
                    "default": ["GO_Biological_Process_2021", "KEGG_2021_Human"]
                }
            },
            "required": ["genes"]
        }
    ),
]


async def handle_annotation_tool(name: str, arguments: Dict[str, Any]) -> str:
    """Handle annotation tool calls."""
    
    if name == "annotate_snps":
        return await annotate_snps(
            rsids=arguments["rsids"],
            include_frequencies=arguments.get("include_frequencies", True)
        )
    
    elif name == "query_gwas_catalog":
        return await query_gwas_catalog(
            rsid=arguments.get("rsid"),
            trait=arguments.get("trait"),
            gene=arguments.get("gene")
        )
    
    elif name == "get_eqtl_data":
        return await get_eqtl_data(
            rsid=arguments["rsid"],
            tissue=arguments.get("tissue"),
            gene=arguments.get("gene")
        )
    
    elif name == "gene_set_enrichment":
        return await gene_set_enrichment(
            genes=arguments["genes"],
            databases=arguments.get("databases", ["GO_Biological_Process_2021", "KEGG_2021_Human"])
        )
    
    raise ValueError(f"Unknown annotation tool: {name}")


async def annotate_snps(
    rsids: List[str],
    include_frequencies: bool = True
) -> str:
    """
    Annotate SNPs using Ensembl VEP API.
    """
    logger.info(f"Annotating {len(rsids)} SNPs")
    
    # Limit batch size
    rsids = rsids[:100]  # API limit
    
    results = []
    
    async with aiohttp.ClientSession() as session:
        for rsid in rsids:
            try:
                # Validate rsID
                rsid = validate_rsid(rsid)
                
                # Query Ensembl VEP
                url = f"{ENSEMBL_REST_API}/variation/human/{rsid}?content-type=application/json"
                
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        
                        annotation = {
                            "rsid": rsid,
                            "chromosome": data.get("mappings", [{}])[0].get("seq_region_name"),
                            "position": data.get("mappings", [{}])[0].get("start"),
                            "alleles": data.get("mappings", [{}])[0].get("allele_string"),
                            "ancestral_allele": data.get("ancestral_allele"),
                            "minor_allele": data.get("minor_allele"),
                            "maf": data.get("MAF"),
                            "clinical_significance": data.get("clinical_significance", []),
                            "most_severe_consequence": data.get("most_severe_consequence"),
                        }
                        
                        # Get mapped genes
                        if "mappings" in data:
                            annotation["genes"] = []
                            for mapping in data["mappings"]:
                                if "consequence_type" in mapping:
                                    annotation["consequence"] = mapping["consequence_type"]
                        
                        # Get population frequencies if requested
                        if include_frequencies and "populations" in data:
                            annotation["population_frequencies"] = {}
                            for pop in data.get("populations", []):
                                pop_name = pop.get("population", "Unknown")
                                annotation["population_frequencies"][pop_name] = pop.get("frequency")
                        
                        results.append(annotation)
                    
                    elif response.status == 429:
                        # Rate limited, wait and retry
                        await asyncio.sleep(1)
                    else:
                        results.append({
                            "rsid": rsid,
                            "error": f"API returned status {response.status}"
                        })
                
                # Rate limiting
                await asyncio.sleep(0.1)
                
            except Exception as e:
                results.append({
                    "rsid": rsid,
                    "error": str(e)
                })
    
    output = {
        "n_queried": len(rsids),
        "n_annotated": len([r for r in results if "error" not in r]),
        "annotations": results
    }
    
    return json.dumps(output, indent=2)


async def query_gwas_catalog(
    rsid: Optional[str] = None,
    trait: Optional[str] = None,
    gene: Optional[str] = None
) -> str:
    """
    Query NHGRI-EBI GWAS Catalog.
    """
    logger.info(f"Querying GWAS Catalog: rsid={rsid}, trait={trait}, gene={gene}")
    
    if not any([rsid, trait, gene]):
        return json.dumps({
            "error": "Must provide at least one of: rsid, trait, or gene"
        }, indent=2)
    
    results = {
        "query": {"rsid": rsid, "trait": trait, "gene": gene},
        "associations": []
    }
    
    try:
        async with aiohttp.ClientSession() as session:
            if rsid:
                rsid = validate_rsid(rsid)
                url = f"{GWAS_CATALOG_API}/singleNucleotidePolymorphisms/{rsid}/associations"
                
                async with session.get(url, headers={"Accept": "application/json"}) as response:
                    if response.status == 200:
                        data = await response.json()
                        
                        if "_embedded" in data and "associations" in data["_embedded"]:
                            for assoc in data["_embedded"]["associations"]:
                                results["associations"].append({
                                    "pvalue": assoc.get("pvalue"),
                                    "pvalue_text": assoc.get("pvalueText"),
                                    "risk_allele_frequency": assoc.get("riskFrequency"),
                                    "beta": assoc.get("betaNum"),
                                    "beta_unit": assoc.get("betaUnit"),
                                    "odds_ratio": assoc.get("orPerCopyNum"),
                                    "trait": assoc.get("_links", {}).get("efoTraits", {}).get("href"),
                                    "study": assoc.get("_links", {}).get("study", {}).get("href"),
                                })
            
            if trait:
                url = f"{GWAS_CATALOG_API}/efoTraits/search/findByEfoUri?uri={trait}"
                # Simplified trait search
                search_url = f"{GWAS_CATALOG_API}/efoTraits/search?query={trait}"
                
                async with session.get(search_url, headers={"Accept": "application/json"}) as response:
                    if response.status == 200:
                        data = await response.json()
                        if "_embedded" in data:
                            results["trait_matches"] = [
                                {"trait": t.get("trait"), "uri": t.get("uri")}
                                for t in data["_embedded"].get("efoTraits", [])[:10]
                            ]
            
            if gene:
                url = f"{GWAS_CATALOG_API}/singleNucleotidePolymorphisms/search/findByGene?geneName={gene}"
                
                async with session.get(url, headers={"Accept": "application/json"}) as response:
                    if response.status == 200:
                        data = await response.json()
                        if "_embedded" in data:
                            results["snps_in_gene"] = [
                                snp.get("rsId")
                                for snp in data["_embedded"].get("singleNucleotidePolymorphisms", [])[:20]
                            ]
        
        results["n_associations"] = len(results["associations"])
        return json.dumps(results, indent=2)
        
    except Exception as e:
        logger.exception(f"GWAS Catalog query failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


async def get_eqtl_data(
    rsid: str,
    tissue: Optional[str] = None,
    gene: Optional[str] = None
) -> str:
    """
    Get eQTL data from GTEx.
    """
    logger.info(f"Querying GTEx eQTLs for: {rsid}")
    
    try:
        rsid = validate_rsid(rsid)
        
        # GTEx API v2
        params = {"snpId": rsid, "datasetId": "gtex_v8"}
        if tissue:
            params["tissueSiteDetailId"] = tissue
        
        async with aiohttp.ClientSession() as session:
            # First, get available tissues
            tissues_url = f"{GTEX_API}/dataset/tissueSiteDetail"
            
            async with session.get(tissues_url) as response:
                available_tissues = []
                if response.status == 200:
                    data = await response.json()
                    available_tissues = [t["tissueSiteDetailId"] for t in data.get("data", [])]
            
            # Query eQTLs
            eqtl_url = f"{GTEX_API}/association/singleTissueEqtl"
            
            results = {
                "rsid": rsid,
                "tissue_queried": tissue,
                "eqtls": [],
                "available_tissues_sample": available_tissues[:20] if available_tissues else []
            }
            
            async with session.get(eqtl_url, params=params) as response:
                if response.status == 200:
                    data = await response.json()
                    
                    for eqtl in data.get("data", [])[:50]:
                        results["eqtls"].append({
                            "gene_symbol": eqtl.get("geneSymbol"),
                            "gene_id": eqtl.get("geneId"),
                            "tissue": eqtl.get("tissueSiteDetailId"),
                            "pvalue": eqtl.get("pValue"),
                            "nes": eqtl.get("nes"),  # Normalized effect size
                            "slope": eqtl.get("slope"),
                        })
                elif response.status == 404:
                    results["note"] = "No eQTL data found for this variant"
                else:
                    results["error"] = f"GTEx API returned status {response.status}"
            
            results["n_eqtls"] = len(results["eqtls"])
            
            # Filter by gene if specified
            if gene and results["eqtls"]:
                results["eqtls"] = [e for e in results["eqtls"] 
                                   if e.get("gene_symbol", "").upper() == gene.upper()]
                results["n_eqtls_filtered"] = len(results["eqtls"])
            
            return json.dumps(results, indent=2)
            
    except Exception as e:
        logger.exception(f"GTEx query failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


async def gene_set_enrichment(
    genes: List[str],
    databases: List[str] = ["GO_Biological_Process_2021", "KEGG_2021_Human"]
) -> str:
    """
    Perform gene set enrichment analysis using Enrichr.
    """
    logger.info(f"Running enrichment analysis for {len(genes)} genes")
    
    ENRICHR_URL = "https://maayanlab.cloud/Enrichr"
    
    try:
        # Limit gene list size
        genes = genes[:500]
        
        async with aiohttp.ClientSession() as session:
            # Step 1: Submit gene list
            genes_str = "\n".join(genes)
            payload = {"list": genes_str, "description": "GWAS gene list"}
            
            async with session.post(f"{ENRICHR_URL}/addList", data=payload) as response:
                if response.status != 200:
                    return json.dumps({
                        "error": f"Enrichr submission failed: {response.status}"
                    }, indent=2)
                
                add_result = await response.json()
                user_list_id = add_result.get("userListId")
            
            # Step 2: Get enrichment results for each database
            results = {
                "n_genes": len(genes),
                "databases": {},
            }
            
            for db in databases:
                async with session.get(
                    f"{ENRICHR_URL}/enrich",
                    params={"userListId": user_list_id, "backgroundType": db}
                ) as response:
                    if response.status == 200:
                        data = await response.json()
                        
                        db_results = []
                        for term in data.get(db, [])[:20]:
                            db_results.append({
                                "term": term[1],
                                "pvalue": term[2],
                                "adjusted_pvalue": term[6],
                                "odds_ratio": term[3],
                                "combined_score": term[4],
                                "genes": term[5],
                                "n_genes_overlap": len(term[5]) if isinstance(term[5], list) else 0
                            })
                        
                        results["databases"][db] = {
                            "n_significant_terms": len([r for r in db_results if r["adjusted_pvalue"] < 0.05]),
                            "top_terms": db_results
                        }
            
            return json.dumps(results, indent=2)
            
    except Exception as e:
        logger.exception(f"Enrichment analysis failed: {e}")
        return json.dumps({
            "error": str(e),
            "type": type(e).__name__
        }, indent=2)


def register_annotation_tools():
    """Return annotation tools for registration."""
    return ANNOTATION_TOOLS
