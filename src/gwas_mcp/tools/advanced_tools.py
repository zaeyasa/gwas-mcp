"""
Advanced Bioinformatics Tools for GWAS MCP Server.

AlphaFold AI structures, Open Targets drug discovery, and OMIM genetic diseases.
"""

import json
import logging
from typing import Any, Dict, List, Optional
import time

import aiohttp

from mcp.types import Tool

logger = logging.getLogger(__name__)

# API endpoints
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api"
OPENTARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"
OMIM_API = "https://api.omim.org/api"  # Note: requires API key for full access

# Cache
_cache: Dict[str, tuple] = {}
CACHE_TTL = 3600


def get_cached(key: str) -> Optional[str]:
    if key in _cache:
        result, timestamp = _cache[key]
        if time.time() - timestamp < CACHE_TTL:
            return result
        del _cache[key]
    return None


def set_cached(key: str, result: str):
    _cache[key] = (result, time.time())
    if len(_cache) > 100:
        oldest = min(_cache.keys(), key=lambda k: _cache[k][1])
        del _cache[oldest]


# Define tools
ADVANCED_TOOLS = [
    # AlphaFold Tools
    Tool(
        name="get_alphafold_structure",
        description="Get AlphaFold AI-predicted protein structure by UniProt ID. Returns structure confidence and download links.",
        inputSchema={
            "type": "object",
            "properties": {
                "uniprot_id": {
                    "type": "string",
                    "description": "UniProt accession ID (e.g., 'P53_HUMAN' or 'P04637')"
                }
            },
            "required": ["uniprot_id"]
        }
    ),
    Tool(
        name="search_alphafold",
        description="Search AlphaFold database for predicted structures by gene name or protein name.",
        inputSchema={
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Gene symbol or protein name (e.g., 'TP53', 'BRCA1')"
                },
                "organism": {
                    "type": "string",
                    "description": "Organism (default: 'Homo sapiens')",
                    "default": "Homo sapiens"
                }
            },
            "required": ["query"]
        }
    ),
    
    # Open Targets Tools
    Tool(
        name="get_drug_targets",
        description="Get drug target information from Open Targets for a gene. Shows drugs in development and approved drugs targeting this gene.",
        inputSchema={
            "type": "object",
            "properties": {
                "gene": {
                    "type": "string",
                    "description": "Gene symbol (e.g., 'EGFR', 'BRAF')"
                }
            },
            "required": ["gene"]
        }
    ),
    Tool(
        name="get_disease_associations",
        description="Get disease associations for a gene from Open Targets with evidence scores.",
        inputSchema={
            "type": "object",
            "properties": {
                "gene": {
                    "type": "string",
                    "description": "Gene symbol"
                },
                "limit": {
                    "type": "integer",
                    "description": "Maximum number of diseases (default: 20)",
                    "default": 20
                }
            },
            "required": ["gene"]
        }
    ),
    Tool(
        name="search_open_targets",
        description="Search Open Targets Platform for genes, diseases, or drugs.",
        inputSchema={
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Search term (gene, disease, or drug name)"
                }
            },
            "required": ["query"]
        }
    ),
    
    # OMIM Tools
    Tool(
        name="search_omim",
        description="Search OMIM (Online Mendelian Inheritance in Man) for genetic diseases and phenotypes.",
        inputSchema={
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Gene symbol or disease name"
                },
                "limit": {
                    "type": "integer",
                    "description": "Maximum results (default: 10)",
                    "default": 10
                }
            },
            "required": ["query"]
        }
    ),
    Tool(
        name="get_gene_diseases",
        description="Get genetic diseases associated with a gene from OMIM and other sources.",
        inputSchema={
            "type": "object",
            "properties": {
                "gene": {
                    "type": "string",
                    "description": "Gene symbol (e.g., 'BRCA1', 'CFTR')"
                }
            },
            "required": ["gene"]
        }
    ),
]


async def handle_advanced_tool(name: str, arguments: Dict[str, Any]) -> str:
    """Handle advanced bioinformatics tool calls."""
    
    if name == "get_alphafold_structure":
        return await get_alphafold_structure(arguments["uniprot_id"])
    elif name == "search_alphafold":
        return await search_alphafold(
            query=arguments["query"],
            organism=arguments.get("organism", "Homo sapiens")
        )
    elif name == "get_drug_targets":
        return await get_drug_targets(arguments["gene"])
    elif name == "get_disease_associations":
        return await get_disease_associations(
            gene=arguments["gene"],
            limit=arguments.get("limit", 20)
        )
    elif name == "search_open_targets":
        return await search_open_targets(arguments["query"])
    elif name == "search_omim":
        return await search_omim(
            query=arguments["query"],
            limit=arguments.get("limit", 10)
        )
    elif name == "get_gene_diseases":
        return await get_gene_diseases(arguments["gene"])
    else:
        raise ValueError(f"Unknown advanced tool: {name}")


async def get_alphafold_structure(uniprot_id: str) -> str:
    """Get AlphaFold predicted structure for a protein."""
    cache_key = f"alphafold:{uniprot_id}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            # Clean up UniProt ID format
            clean_id = uniprot_id.split("_")[0] if "_" in uniprot_id else uniprot_id
            
            url = f"{ALPHAFOLD_API}/prediction/{clean_id}"
            
            async with session.get(url) as response:
                if response.status == 404:
                    return json.dumps({
                        "uniprot_id": uniprot_id,
                        "error": "No AlphaFold prediction found for this protein"
                    }, indent=2)
                if response.status != 200:
                    return json.dumps({"error": f"AlphaFold API error: {response.status}"}, indent=2)
                
                data = await response.json()
                
                if isinstance(data, list) and len(data) > 0:
                    entry = data[0]
                else:
                    entry = data
                
                result = json.dumps({
                    "uniprot_id": entry.get("uniprotAccession"),
                    "entry_name": entry.get("uniprotId"),
                    "gene": entry.get("gene"),
                    "organism": entry.get("organismScientificName"),
                    "sequence_length": entry.get("uniprotSequence", ""),
                    "model_created": entry.get("modelCreatedDate"),
                    "latest_version": entry.get("latestVersion"),
                    "average_plddt": entry.get("globalMetricValue"),
                    "confidence_description": "pLDDT > 90: Very high confidence, 70-90: Confident, 50-70: Low confidence, < 50: Very low",
                    "structure_url": entry.get("cifUrl"),
                    "pdb_url": entry.get("pdbUrl"),
                    "viewer_url": f"https://alphafold.ebi.ac.uk/entry/{entry.get('uniprotAccession')}",
                    "source": "AlphaFold DB"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"AlphaFold fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def search_alphafold(query: str, organism: str = "Homo sapiens") -> str:
    """Search AlphaFold by gene name (via UniProt mapping)."""
    cache_key = f"alphafold_search:{query}:{organism}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            # First search UniProt to get IDs
            uniprot_url = "https://rest.uniprot.org/uniprotkb/search"
            params = {
                "query": f'(gene:{query} OR protein_name:{query}) AND organism_name:"{organism}"',
                "format": "json",
                "size": 5,
                "fields": "accession,id,gene_names,protein_name,organism_name"
            }
            
            async with session.get(uniprot_url, params=params) as response:
                if response.status != 200:
                    return json.dumps({"error": f"UniProt search failed: {response.status}"}, indent=2)
                
                data = await response.json()
                
                results = []
                for entry in data.get("results", []):
                    uniprot_id = entry.get("primaryAccession")
                    
                    # Check if AlphaFold has this structure
                    af_url = f"{ALPHAFOLD_API}/prediction/{uniprot_id}"
                    async with session.get(af_url) as af_response:
                        has_structure = af_response.status == 200
                        
                        results.append({
                            "uniprot_id": uniprot_id,
                            "gene": [g.get("geneName", {}).get("value") for g in entry.get("genes", []) if g.get("geneName")][:1],
                            "protein_name": entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"),
                            "organism": entry.get("organism", {}).get("scientificName"),
                            "has_alphafold": has_structure,
                            "alphafold_url": f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}" if has_structure else None
                        })
                
                result = json.dumps({
                    "query": query,
                    "organism": organism,
                    "n_results": len(results),
                    "results": results,
                    "source": "AlphaFold + UniProt"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"AlphaFold search failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def get_drug_targets(gene: str) -> str:
    """Get drug target information from Open Targets."""
    cache_key = f"opentargets_drugs:{gene}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            # GraphQL query
            query = """
            query DrugTargets($ensemblId: String!) {
                target(ensemblId: $ensemblId) {
                    id
                    approvedSymbol
                    approvedName
                    knownDrugs {
                        uniqueDrugs
                        rows {
                            drug {
                                id
                                name
                                drugType
                                maximumClinicalTrialPhase
                                isApproved
                            }
                            mechanismOfAction
                            disease {
                                id
                                name
                            }
                        }
                    }
                }
            }
            """
            
            # First get Ensembl ID for the gene
            search_query = """
            query SearchGene($queryString: String!) {
                search(queryString: $queryString, entityNames: ["target"], page: {size: 1, index: 0}) {
                    hits {
                        id
                        name
                        entity
                    }
                }
            }
            """
            
            # Search for gene
            async with session.post(
                OPENTARGETS_API,
                json={"query": search_query, "variables": {"queryString": gene}},
                headers={"Content-Type": "application/json"}
            ) as search_response:
                if search_response.status != 200:
                    return json.dumps({"error": f"Open Targets search failed: {search_response.status}"}, indent=2)
                
                search_data = await search_response.json()
                hits = search_data.get("data", {}).get("search", {}).get("hits", [])
                
                if not hits:
                    return json.dumps({
                        "gene": gene,
                        "error": "Gene not found in Open Targets"
                    }, indent=2)
                
                ensembl_id = hits[0]["id"]
            
            # Get drug targets
            async with session.post(
                OPENTARGETS_API,
                json={"query": query, "variables": {"ensemblId": ensembl_id}},
                headers={"Content-Type": "application/json"}
            ) as response:
                if response.status != 200:
                    return json.dumps({"error": f"Open Targets API error: {response.status}"}, indent=2)
                
                data = await response.json()
                target = data.get("data", {}).get("target", {})
                
                if not target:
                    return json.dumps({
                        "gene": gene,
                        "error": "No target data found"
                    }, indent=2)
                
                known_drugs = target.get("knownDrugs", {})
                drugs = []
                
                for row in known_drugs.get("rows", [])[:15]:
                    drug_info = row.get("drug", {})
                    drugs.append({
                        "drug_name": drug_info.get("name"),
                        "drug_type": drug_info.get("drugType"),
                        "max_phase": drug_info.get("maximumClinicalTrialPhase"),
                        "is_approved": drug_info.get("isApproved"),
                        "mechanism": row.get("mechanismOfAction"),
                        "indication": row.get("disease", {}).get("name")
                    })
                
                result = json.dumps({
                    "gene": gene,
                    "ensembl_id": ensembl_id,
                    "approved_name": target.get("approvedName"),
                    "total_drugs": known_drugs.get("uniqueDrugs", 0),
                    "drugs": drugs,
                    "open_targets_url": f"https://platform.opentargets.org/target/{ensembl_id}",
                    "source": "Open Targets Platform"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"Open Targets drug fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def get_disease_associations(gene: str, limit: int = 20) -> str:
    """Get disease associations from Open Targets."""
    cache_key = f"opentargets_diseases:{gene}:{limit}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            # Search for gene first
            search_query = """
            query SearchGene($queryString: String!) {
                search(queryString: $queryString, entityNames: ["target"], page: {size: 1, index: 0}) {
                    hits { id }
                }
            }
            """
            
            async with session.post(
                OPENTARGETS_API,
                json={"query": search_query, "variables": {"queryString": gene}},
                headers={"Content-Type": "application/json"}
            ) as search_response:
                search_data = await search_response.json()
                hits = search_data.get("data", {}).get("search", {}).get("hits", [])
                
                if not hits:
                    return json.dumps({"gene": gene, "error": "Gene not found"}, indent=2)
                
                ensembl_id = hits[0]["id"]
            
            # Get disease associations
            disease_query = """
            query DiseaseAssociations($ensemblId: String!, $size: Int!) {
                target(ensemblId: $ensemblId) {
                    id
                    approvedSymbol
                    associatedDiseases(page: {size: $size, index: 0}) {
                        count
                        rows {
                            disease {
                                id
                                name
                            }
                            score
                            datatypeScores {
                                id
                                score
                            }
                        }
                    }
                }
            }
            """
            
            async with session.post(
                OPENTARGETS_API,
                json={"query": disease_query, "variables": {"ensemblId": ensembl_id, "size": limit}},
                headers={"Content-Type": "application/json"}
            ) as response:
                data = await response.json()
                target = data.get("data", {}).get("target", {})
                
                assoc_data = target.get("associatedDiseases", {})
                diseases = []
                
                for row in assoc_data.get("rows", []):
                    disease = row.get("disease", {})
                    diseases.append({
                        "disease": disease.get("name"),
                        "disease_id": disease.get("id"),
                        "overall_score": round(row.get("score", 0), 3),
                        "evidence_types": {d["id"]: round(d["score"], 3) for d in row.get("datatypeScores", [])}
                    })
                
                result = json.dumps({
                    "gene": gene,
                    "ensembl_id": ensembl_id,
                    "total_associations": assoc_data.get("count", 0),
                    "diseases": diseases,
                    "source": "Open Targets Platform"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"Open Targets disease fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def search_open_targets(query: str) -> str:
    """Search Open Targets Platform."""
    cache_key = f"opentargets_search:{query}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            search_query = """
            query Search($queryString: String!) {
                search(queryString: $queryString, page: {size: 15, index: 0}) {
                    total
                    hits {
                        id
                        name
                        entity
                        description
                    }
                }
            }
            """
            
            async with session.post(
                OPENTARGETS_API,
                json={"query": search_query, "variables": {"queryString": query}},
                headers={"Content-Type": "application/json"}
            ) as response:
                if response.status != 200:
                    return json.dumps({"error": f"Open Targets API error: {response.status}"}, indent=2)
                
                data = await response.json()
                search_results = data.get("data", {}).get("search", {})
                
                results = []
                for hit in search_results.get("hits", []):
                    entity_type = hit.get("entity")
                    url_type = "target" if entity_type == "target" else "disease" if entity_type == "disease" else "drug"
                    
                    results.append({
                        "id": hit.get("id"),
                        "name": hit.get("name"),
                        "type": entity_type,
                        "description": hit.get("description", "")[:200] + "..." if len(hit.get("description", "")) > 200 else hit.get("description"),
                        "url": f"https://platform.opentargets.org/{url_type}/{hit.get('id')}"
                    })
                
                result = json.dumps({
                    "query": query,
                    "total": search_results.get("total", 0),
                    "results": results,
                    "source": "Open Targets Platform"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"Open Targets search failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def search_omim(query: str, limit: int = 10) -> str:
    """Search OMIM via NCBI integration."""
    cache_key = f"omim_search:{query}:{limit}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            # Use NCBI's OMIM integration
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {
                "db": "omim",
                "term": query,
                "retmode": "json",
                "retmax": limit
            }
            
            async with session.get(search_url, params=params) as response:
                if response.status != 200:
                    return json.dumps({"error": f"OMIM search failed: {response.status}"}, indent=2)
                
                search_data = await response.json()
                ids = search_data.get("esearchresult", {}).get("idlist", [])
                
                if not ids:
                    return json.dumps({
                        "query": query,
                        "n_results": 0,
                        "message": "No OMIM entries found"
                    }, indent=2)
                
                # Get summaries
                summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                summary_params = {
                    "db": "omim",
                    "id": ",".join(ids),
                    "retmode": "json"
                }
                
                async with session.get(summary_url, params=summary_params) as sum_response:
                    summary_data = await sum_response.json()
                    
                    results = []
                    for omim_id in ids:
                        entry = summary_data.get("result", {}).get(omim_id, {})
                        if entry:
                            results.append({
                                "omim_id": omim_id,
                                "title": entry.get("title"),
                                "prefix": entry.get("prefix"),  # *, #, %, +, etc.
                                "mim_type": entry.get("mim_type"),
                                "omim_url": f"https://omim.org/entry/{omim_id}"
                            })
                    
                    result = json.dumps({
                        "query": query,
                        "n_results": len(results),
                        "results": results,
                        "prefix_legend": {
                            "*": "Gene",
                            "#": "Phenotype (disease)",
                            "%": "Phenotype (disease) or locus",
                            "+": "Gene and phenotype combined",
                            "none": "Other"
                        },
                        "source": "OMIM via NCBI"
                    }, indent=2)
                    
                    set_cached(cache_key, result)
                    return result
                    
    except Exception as e:
        logger.exception(f"OMIM search failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def get_gene_diseases(gene: str) -> str:
    """Get diseases associated with a gene from multiple sources."""
    cache_key = f"gene_diseases:{gene}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        # Combine OMIM and Open Targets data
        omim_result = await search_omim(gene, limit=10)
        ot_result = await get_disease_associations(gene, limit=10)
        
        omim_data = json.loads(omim_result)
        ot_data = json.loads(ot_result)
        
        result = json.dumps({
            "gene": gene,
            "omim_entries": omim_data.get("results", []) if "results" in omim_data else [],
            "open_targets_diseases": ot_data.get("diseases", []) if "diseases" in ot_data else [],
            "sources": ["OMIM", "Open Targets Platform"]
        }, indent=2)
        
        set_cached(cache_key, result)
        return result
        
    except Exception as e:
        logger.exception(f"Gene diseases fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


def register_advanced_tools():
    """Return advanced tools for registration."""
    return ADVANCED_TOOLS
