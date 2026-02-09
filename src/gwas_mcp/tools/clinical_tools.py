"""
Clinical and Network Tools for GWAS MCP Server.

Tools for ClinVar clinical variant data and STRING protein-protein interactions.
"""

import json
import logging
from typing import Any, Dict, List, Optional
from functools import lru_cache
import time

import aiohttp

from mcp.types import Tool

logger = logging.getLogger(__name__)

# API endpoints
CLINVAR_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
STRING_API = "https://string-db.org/api"

# Simple in-memory cache
_cache: Dict[str, tuple] = {}  # {key: (result, timestamp)}
CACHE_TTL = 3600  # 1 hour


def get_cached(key: str) -> Optional[str]:
    """Get cached result if not expired."""
    if key in _cache:
        result, timestamp = _cache[key]
        if time.time() - timestamp < CACHE_TTL:
            logger.info(f"Cache hit for: {key[:50]}...")
            return result
        else:
            del _cache[key]
    return None


def set_cached(key: str, result: str):
    """Cache a result."""
    _cache[key] = (result, time.time())
    # Keep cache size reasonable
    if len(_cache) > 100:
        oldest_key = min(_cache.keys(), key=lambda k: _cache[k][1])
        del _cache[oldest_key]


# Define Clinical/Network tools
CLINICAL_TOOLS = [
    Tool(
        name="search_clinvar",
        description="Search ClinVar for clinical variant interpretations. Find pathogenic/benign classifications for genetic variants.",
        inputSchema={
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Search query - can be gene name (e.g., 'BRCA1'), rsID (e.g., 'rs80357906'), or condition (e.g., 'breast cancer')"
                },
                "limit": {
                    "type": "integer",
                    "description": "Maximum number of results (default: 10)",
                    "default": 10
                }
            },
            "required": ["query"]
        }
    ),
    Tool(
        name="get_clinvar_variant",
        description="Get detailed ClinVar information for a specific variant by rsID or ClinVar ID.",
        inputSchema={
            "type": "object",
            "properties": {
                "variant_id": {
                    "type": "string",
                    "description": "Variant identifier - rsID (e.g., 'rs80357906') or ClinVar ID"
                }
            },
            "required": ["variant_id"]
        }
    ),
    Tool(
        name="get_protein_interactions",
        description="Get protein-protein interaction network from STRING database. Shows what proteins interact with your protein of interest.",
        inputSchema={
            "type": "object",
            "properties": {
                "protein": {
                    "type": "string",
                    "description": "Protein or gene name (e.g., 'TP53', 'BRCA1')"
                },
                "species": {
                    "type": "integer",
                    "description": "NCBI species taxonomy ID (default: 9606 for human)",
                    "default": 9606
                },
                "limit": {
                    "type": "integer",
                    "description": "Maximum number of interaction partners to return (default: 20)",
                    "default": 20
                },
                "score_threshold": {
                    "type": "number",
                    "description": "Minimum interaction confidence score 0-1 (default: 0.4 = medium confidence)",
                    "default": 0.4
                }
            },
            "required": ["protein"]
        }
    ),
    Tool(
        name="get_interaction_network",
        description="Get interaction network between a list of proteins from STRING. Shows how multiple proteins interact with each other.",
        inputSchema={
            "type": "object",
            "properties": {
                "proteins": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of protein/gene names (e.g., ['TP53', 'BRCA1', 'ATM'])"
                },
                "species": {
                    "type": "integer",
                    "description": "NCBI species taxonomy ID (default: 9606 for human)",
                    "default": 9606
                }
            },
            "required": ["proteins"]
        }
    ),
    Tool(
        name="get_functional_enrichment",
        description="Get functional enrichment analysis for a list of proteins using STRING. Returns enriched GO terms, pathways, and domains.",
        inputSchema={
            "type": "object",
            "properties": {
                "proteins": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of protein/gene names"
                },
                "species": {
                    "type": "integer",
                    "description": "NCBI species taxonomy ID (default: 9606 for human)",
                    "default": 9606
                }
            },
            "required": ["proteins"]
        }
    ),
]


async def handle_clinical_tool(name: str, arguments: Dict[str, Any]) -> str:
    """Handle clinical/network tool calls."""
    
    if name == "search_clinvar":
        return await search_clinvar(
            query=arguments["query"],
            limit=arguments.get("limit", 10)
        )
    elif name == "get_clinvar_variant":
        return await get_clinvar_variant(arguments["variant_id"])
    elif name == "get_protein_interactions":
        return await get_protein_interactions(
            protein=arguments["protein"],
            species=arguments.get("species", 9606),
            limit=arguments.get("limit", 20),
            score_threshold=arguments.get("score_threshold", 0.4)
        )
    elif name == "get_interaction_network":
        return await get_interaction_network(
            proteins=arguments["proteins"],
            species=arguments.get("species", 9606)
        )
    elif name == "get_functional_enrichment":
        return await get_functional_enrichment(
            proteins=arguments["proteins"],
            species=arguments.get("species", 9606)
        )
    else:
        raise ValueError(f"Unknown clinical tool: {name}")


async def search_clinvar(query: str, limit: int = 10) -> str:
    """
    Search ClinVar for clinical variant information.
    """
    cache_key = f"clinvar_search:{query}:{limit}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            # Search ClinVar
            search_url = f"{CLINVAR_API}/esearch.fcgi"
            search_params = {
                "db": "clinvar",
                "term": query,
                "retmode": "json",
                "retmax": limit
            }
            
            async with session.get(search_url, params=search_params) as response:
                if response.status != 200:
                    return json.dumps({"error": f"ClinVar search failed: {response.status}"}, indent=2)
                
                search_data = await response.json()
                ids = search_data.get("esearchresult", {}).get("idlist", [])
                
                if not ids:
                    return json.dumps({
                        "query": query,
                        "n_results": 0,
                        "message": "No ClinVar entries found"
                    }, indent=2)
                
                # Get summaries
                summary_url = f"{CLINVAR_API}/esummary.fcgi"
                summary_params = {
                    "db": "clinvar",
                    "id": ",".join(ids),
                    "retmode": "json"
                }
                
                async with session.get(summary_url, params=summary_params) as sum_response:
                    if sum_response.status != 200:
                        return json.dumps({"error": f"ClinVar summary failed: {sum_response.status}"}, indent=2)
                    
                    summary_data = await sum_response.json()
                    
                    results = []
                    for var_id in ids:
                        var_info = summary_data.get("result", {}).get(var_id, {})
                        if var_info:
                            # Parse clinical significance
                            clin_sig = var_info.get("clinical_significance", {})
                            
                            results.append({
                                "clinvar_id": var_id,
                                "title": var_info.get("title"),
                                "gene": var_info.get("genes", [{}])[0].get("symbol") if var_info.get("genes") else None,
                                "clinical_significance": clin_sig.get("description") if isinstance(clin_sig, dict) else str(clin_sig),
                                "review_status": var_info.get("clinical_significance", {}).get("review_status") if isinstance(clin_sig, dict) else None,
                                "condition": var_info.get("trait_set", [{}])[0].get("trait_name") if var_info.get("trait_set") else None,
                                "variant_type": var_info.get("variant_type")
                            })
                    
                    result = json.dumps({
                        "query": query,
                        "n_results": len(results),
                        "results": results
                    }, indent=2)
                    
                    set_cached(cache_key, result)
                    return result
                    
    except Exception as e:
        logger.exception(f"ClinVar search failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def get_clinvar_variant(variant_id: str) -> str:
    """
    Get detailed ClinVar information for a variant.
    """
    cache_key = f"clinvar_variant:{variant_id}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            # If it's an rsID, search first
            if variant_id.startswith("rs"):
                search_url = f"{CLINVAR_API}/esearch.fcgi"
                search_params = {
                    "db": "clinvar",
                    "term": f"{variant_id}[Variant ID]",
                    "retmode": "json",
                    "retmax": 5
                }
                
                async with session.get(search_url, params=search_params) as response:
                    if response.status == 200:
                        data = await response.json()
                        ids = data.get("esearchresult", {}).get("idlist", [])
                        if ids:
                            variant_id = ids[0]
                        else:
                            return json.dumps({
                                "variant_id": variant_id,
                                "error": "Variant not found in ClinVar"
                            }, indent=2)
            
            # Get detailed info
            summary_url = f"{CLINVAR_API}/esummary.fcgi"
            summary_params = {
                "db": "clinvar",
                "id": variant_id,
                "retmode": "json"
            }
            
            async with session.get(summary_url, params=summary_params) as response:
                if response.status != 200:
                    return json.dumps({"error": f"ClinVar API error: {response.status}"}, indent=2)
                
                data = await response.json()
                var_info = data.get("result", {}).get(str(variant_id), {})
                
                if not var_info:
                    return json.dumps({
                        "variant_id": variant_id,
                        "error": "Variant not found"
                    }, indent=2)
                
                clin_sig = var_info.get("clinical_significance", {})
                
                result = json.dumps({
                    "clinvar_id": variant_id,
                    "title": var_info.get("title"),
                    "gene": var_info.get("genes", [{}])[0].get("symbol") if var_info.get("genes") else None,
                    "clinical_significance": clin_sig.get("description") if isinstance(clin_sig, dict) else str(clin_sig),
                    "review_status": clin_sig.get("review_status") if isinstance(clin_sig, dict) else None,
                    "conditions": [t.get("trait_name") for t in var_info.get("trait_set", [])],
                    "variant_type": var_info.get("variant_type"),
                    "chromosome": var_info.get("chr"),
                    "start": var_info.get("start"),
                    "stop": var_info.get("stop"),
                    "source": "ClinVar"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"ClinVar variant fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def get_protein_interactions(
    protein: str,
    species: int = 9606,
    limit: int = 20,
    score_threshold: float = 0.4
) -> str:
    """
    Get protein-protein interactions from STRING.
    """
    cache_key = f"string_interactions:{protein}:{species}:{limit}:{score_threshold}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            url = f"{STRING_API}/json/network"
            params = {
                "identifiers": protein,
                "species": species,
                "limit": limit,
                "required_score": int(score_threshold * 1000),  # STRING uses 0-1000 scale
                "caller_identity": "gwas_mcp"
            }
            
            async with session.get(url, params=params) as response:
                if response.status != 200:
                    return json.dumps({"error": f"STRING API error: {response.status}"}, indent=2)
                
                data = await response.json()
                
                if not data:
                    return json.dumps({
                        "protein": protein,
                        "n_interactions": 0,
                        "message": "No interactions found"
                    }, indent=2)
                
                # Parse interactions
                interactions = []
                seen_partners = set()
                
                for edge in data:
                    partner = edge.get("preferredName_B") if edge.get("preferredName_A") == protein.upper() else edge.get("preferredName_A")
                    
                    if partner and partner not in seen_partners:
                        seen_partners.add(partner)
                        interactions.append({
                            "partner": partner,
                            "score": edge.get("score", 0) / 1000,  # Convert to 0-1
                            "string_id": edge.get("stringId_B") if edge.get("preferredName_A") == protein.upper() else edge.get("stringId_A")
                        })
                
                # Sort by score
                interactions.sort(key=lambda x: x["score"], reverse=True)
                
                result = json.dumps({
                    "protein": protein,
                    "species": species,
                    "score_threshold": score_threshold,
                    "n_interactions": len(interactions),
                    "interacting_proteins": interactions[:limit],
                    "source": "STRING"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"STRING interaction fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def get_interaction_network(proteins: List[str], species: int = 9606) -> str:
    """
    Get interaction network between a list of proteins.
    """
    cache_key = f"string_network:{','.join(sorted(proteins))}:{species}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            url = f"{STRING_API}/json/network"
            params = {
                "identifiers": "%0d".join(proteins),  # Newline-separated
                "species": species,
                "required_score": 400,  # Medium confidence
                "caller_identity": "gwas_mcp"
            }
            
            async with session.get(url, params=params) as response:
                if response.status != 200:
                    return json.dumps({"error": f"STRING API error: {response.status}"}, indent=2)
                
                data = await response.json()
                
                # Parse network
                edges = []
                nodes = set()
                
                for edge in data:
                    protein_a = edge.get("preferredName_A")
                    protein_b = edge.get("preferredName_B")
                    
                    if protein_a and protein_b:
                        nodes.add(protein_a)
                        nodes.add(protein_b)
                        edges.append({
                            "protein_a": protein_a,
                            "protein_b": protein_b,
                            "score": edge.get("score", 0) / 1000
                        })
                
                result = json.dumps({
                    "query_proteins": proteins,
                    "n_nodes": len(nodes),
                    "n_edges": len(edges),
                    "nodes": list(nodes),
                    "edges": edges,
                    "source": "STRING"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"STRING network fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def get_functional_enrichment(proteins: List[str], species: int = 9606) -> str:
    """
    Get functional enrichment analysis from STRING.
    """
    cache_key = f"string_enrichment:{','.join(sorted(proteins))}:{species}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            url = f"{STRING_API}/json/enrichment"
            params = {
                "identifiers": "%0d".join(proteins),
                "species": species,
                "caller_identity": "gwas_mcp"
            }
            
            async with session.get(url, params=params) as response:
                if response.status != 200:
                    return json.dumps({"error": f"STRING API error: {response.status}"}, indent=2)
                
                data = await response.json()
                
                if not data:
                    return json.dumps({
                        "proteins": proteins,
                        "n_enrichments": 0,
                        "message": "No enrichments found"
                    }, indent=2)
                
                # Group by category
                categories = {}
                for term in data[:50]:  # Limit to top 50
                    category = term.get("category", "Other")
                    if category not in categories:
                        categories[category] = []
                    
                    categories[category].append({
                        "term": term.get("term"),
                        "description": term.get("description"),
                        "p_value": term.get("p_value"),
                        "fdr": term.get("fdr"),
                        "genes_in_term": term.get("number_of_genes"),
                        "genes_matched": term.get("inputGenes", "").split(",") if term.get("inputGenes") else []
                    })
                
                result = json.dumps({
                    "proteins": proteins,
                    "n_categories": len(categories),
                    "enrichment_by_category": categories,
                    "source": "STRING"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"STRING enrichment failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


def register_clinical_tools():
    """Return clinical tools for registration."""
    return CLINICAL_TOOLS
