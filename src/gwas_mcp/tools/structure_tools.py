"""
Structure, Pathway, and Pharmacogenomics Tools for GWAS MCP Server.

Tools for PDB protein structures, KEGG pathways, and PharmGKB drug-gene interactions.
"""

import json
import logging
from typing import Any, Dict, List, Optional
import time

import aiohttp

from mcp.types import Tool

logger = logging.getLogger(__name__)

# API endpoints
PDB_API = "https://data.rcsb.org/rest/v1"
PDB_SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2"
KEGG_API = "https://rest.kegg.jp"
PHARMGKB_API = "https://api.pharmgkb.org/v1"

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
STRUCTURE_TOOLS = [
    Tool(
        name="search_pdb_structures",
        description="Search PDB for protein 3D structures by protein name, gene name, or UniProt ID.",
        inputSchema={
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Search query - protein name, gene symbol, or UniProt ID"
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
        name="get_pdb_structure",
        description="Get detailed information about a PDB structure by PDB ID.",
        inputSchema={
            "type": "object",
            "properties": {
                "pdb_id": {
                    "type": "string",
                    "description": "PDB ID (e.g., '1TUP', '6VXX')"
                }
            },
            "required": ["pdb_id"]
        }
    ),
    Tool(
        name="search_kegg_pathway",
        description="Search KEGG for metabolic and signaling pathways by name or gene.",
        inputSchema={
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Pathway name or gene symbol"
                },
                "organism": {
                    "type": "string",
                    "description": "Organism code (default: 'hsa' for human)",
                    "default": "hsa"
                }
            },
            "required": ["query"]
        }
    ),
    Tool(
        name="get_kegg_pathway",
        description="Get detailed KEGG pathway information including genes and description.",
        inputSchema={
            "type": "object",
            "properties": {
                "pathway_id": {
                    "type": "string",
                    "description": "KEGG pathway ID (e.g., 'hsa04110' for cell cycle)"
                }
            },
            "required": ["pathway_id"]
        }
    ),
    Tool(
        name="get_gene_pathways",
        description="Find all KEGG pathways that a gene participates in.",
        inputSchema={
            "type": "object",
            "properties": {
                "gene": {
                    "type": "string",
                    "description": "Gene symbol (e.g., 'TP53', 'BRCA1')"
                }
            },
            "required": ["gene"]
        }
    ),
    Tool(
        name="search_pharmgkb",
        description="Search PharmGKB for drug-gene interactions and pharmacogenomics data.",
        inputSchema={
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Gene name, drug name, or rsID"
                },
                "resource_type": {
                    "type": "string",
                    "description": "Type of resource: 'gene', 'drug', 'variant', or 'all' (default: 'all')",
                    "default": "all"
                }
            },
            "required": ["query"]
        }
    ),
    Tool(
        name="get_drug_gene_interactions",
        description="Get drug-gene interaction annotations from PharmGKB for a specific gene.",
        inputSchema={
            "type": "object",
            "properties": {
                "gene": {
                    "type": "string",
                    "description": "Gene symbol (e.g., 'CYP2D6', 'VKORC1')"
                }
            },
            "required": ["gene"]
        }
    ),
]


async def handle_structure_tool(name: str, arguments: Dict[str, Any]) -> str:
    """Handle structure/pathway/pharmacogenomics tool calls."""
    
    if name == "search_pdb_structures":
        return await search_pdb_structures(
            query=arguments["query"],
            limit=arguments.get("limit", 10)
        )
    elif name == "get_pdb_structure":
        return await get_pdb_structure(arguments["pdb_id"])
    elif name == "search_kegg_pathway":
        return await search_kegg_pathway(
            query=arguments["query"],
            organism=arguments.get("organism", "hsa")
        )
    elif name == "get_kegg_pathway":
        return await get_kegg_pathway(arguments["pathway_id"])
    elif name == "get_gene_pathways":
        return await get_gene_pathways(arguments["gene"])
    elif name == "search_pharmgkb":
        return await search_pharmgkb(
            query=arguments["query"],
            resource_type=arguments.get("resource_type", "all")
        )
    elif name == "get_drug_gene_interactions":
        return await get_drug_gene_interactions(arguments["gene"])
    else:
        raise ValueError(f"Unknown structure tool: {name}")


async def search_pdb_structures(query: str, limit: int = 10) -> str:
    """Search PDB for protein structures."""
    cache_key = f"pdb_search:{query}:{limit}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            # Use text search API
            search_query = {
                "query": {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": query
                    }
                },
                "return_type": "entry",
                "request_options": {
                    "paginate": {
                        "start": 0,
                        "rows": limit
                    },
                    "results_content_type": ["experimental"]
                }
            }
            
            url = f"{PDB_SEARCH_API}/query"
            headers = {"Content-Type": "application/json"}
            
            async with session.post(url, json=search_query, headers=headers) as response:
                if response.status != 200:
                    # Fallback to simpler search
                    return await _pdb_simple_search(session, query, limit)
                
                data = await response.json()
                
                pdb_ids = [hit["identifier"] for hit in data.get("result_set", [])]
                
                if not pdb_ids:
                    return json.dumps({
                        "query": query,
                        "n_results": 0,
                        "message": "No PDB structures found"
                    }, indent=2)
                
                # Get basic info for each
                structures = []
                for pdb_id in pdb_ids[:limit]:
                    info_url = f"{PDB_API}/core/entry/{pdb_id}"
                    async with session.get(info_url) as info_response:
                        if info_response.status == 200:
                            info = await info_response.json()
                            structures.append({
                                "pdb_id": pdb_id,
                                "title": info.get("struct", {}).get("title"),
                                "method": info.get("exptl", [{}])[0].get("method") if info.get("exptl") else None,
                                "resolution": info.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
                                "deposit_date": info.get("rcsb_accession_info", {}).get("deposit_date")
                            })
                
                result = json.dumps({
                    "query": query,
                    "n_results": len(structures),
                    "structures": structures,
                    "source": "PDB/RCSB"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"PDB search failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def _pdb_simple_search(session, query: str, limit: int) -> str:
    """Fallback simple PDB search."""
    try:
        url = f"https://www.rcsb.org/search/suggest/entry/{query}"
        async with session.get(url) as response:
            if response.status == 200:
                suggestions = await response.json()
                return json.dumps({
                    "query": query,
                    "n_results": len(suggestions[:limit]),
                    "suggestions": suggestions[:limit],
                    "source": "PDB/RCSB"
                }, indent=2)
    except:
        pass
    return json.dumps({"query": query, "n_results": 0, "message": "Search failed"}, indent=2)


async def get_pdb_structure(pdb_id: str) -> str:
    """Get detailed PDB structure information."""
    cache_key = f"pdb_structure:{pdb_id}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            url = f"{PDB_API}/core/entry/{pdb_id.upper()}"
            
            async with session.get(url) as response:
                if response.status == 404:
                    return json.dumps({"pdb_id": pdb_id, "error": "Structure not found"}, indent=2)
                if response.status != 200:
                    return json.dumps({"error": f"PDB API error: {response.status}"}, indent=2)
                
                data = await response.json()
                
                # Get polymer entities (proteins, DNA, etc.)
                entities_url = f"{PDB_API}/core/polymer_entity/{pdb_id.upper()}"
                entities = []
                try:
                    async with session.get(entities_url) as ent_response:
                        if ent_response.status == 200:
                            ent_data = await ent_response.json()
                            for ent in ent_data if isinstance(ent_data, list) else [ent_data]:
                                entities.append({
                                    "entity_id": ent.get("rcsb_id"),
                                    "type": ent.get("entity_poly", {}).get("type"),
                                    "description": ent.get("rcsb_polymer_entity", {}).get("pdbx_description")
                                })
                except:
                    pass
                
                result = json.dumps({
                    "pdb_id": pdb_id.upper(),
                    "title": data.get("struct", {}).get("title"),
                    "method": data.get("exptl", [{}])[0].get("method") if data.get("exptl") else None,
                    "resolution": data.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
                    "deposit_date": data.get("rcsb_accession_info", {}).get("deposit_date"),
                    "release_date": data.get("rcsb_accession_info", {}).get("initial_release_date"),
                    "keywords": data.get("struct_keywords", {}).get("pdbx_keywords"),
                    "organism": data.get("rcsb_entry_info", {}).get("polymer_entity_count_protein"),
                    "entities": entities[:5],
                    "pdb_url": f"https://www.rcsb.org/structure/{pdb_id.upper()}",
                    "source": "PDB/RCSB"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"PDB structure fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def search_kegg_pathway(query: str, organism: str = "hsa") -> str:
    """Search KEGG for pathways."""
    cache_key = f"kegg_search:{query}:{organism}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            # Search pathways
            url = f"{KEGG_API}/find/pathway/{query}"
            
            async with session.get(url) as response:
                if response.status != 200:
                    return json.dumps({"error": f"KEGG API error: {response.status}"}, indent=2)
                
                text = await response.text()
                
                if not text.strip():
                    return json.dumps({
                        "query": query,
                        "n_results": 0,
                        "message": "No pathways found"
                    }, indent=2)
                
                pathways = []
                for line in text.strip().split("\n"):
                    if line:
                        parts = line.split("\t")
                        if len(parts) >= 2:
                            pathway_id = parts[0].replace("path:", "")
                            # Filter for organism
                            if organism == "hsa" or pathway_id.startswith(organism):
                                pathways.append({
                                    "pathway_id": pathway_id,
                                    "name": parts[1],
                                    "url": f"https://www.kegg.jp/pathway/{pathway_id}"
                                })
                
                result = json.dumps({
                    "query": query,
                    "organism": organism,
                    "n_results": len(pathways),
                    "pathways": pathways[:20],
                    "source": "KEGG"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"KEGG search failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def get_kegg_pathway(pathway_id: str) -> str:
    """Get detailed KEGG pathway information."""
    cache_key = f"kegg_pathway:{pathway_id}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            url = f"{KEGG_API}/get/{pathway_id}"
            
            async with session.get(url) as response:
                if response.status == 404:
                    return json.dumps({"pathway_id": pathway_id, "error": "Pathway not found"}, indent=2)
                if response.status != 200:
                    return json.dumps({"error": f"KEGG API error: {response.status}"}, indent=2)
                
                text = await response.text()
                
                # Parse KEGG flat file format
                pathway_info = {
                    "pathway_id": pathway_id,
                    "name": None,
                    "description": None,
                    "genes": [],
                    "diseases": [],
                    "drugs": [],
                    "url": f"https://www.kegg.jp/pathway/{pathway_id}",
                    "source": "KEGG"
                }
                
                current_section = None
                for line in text.split("\n"):
                    if line.startswith("NAME"):
                        pathway_info["name"] = line[12:].strip()
                    elif line.startswith("DESCRIPTION"):
                        pathway_info["description"] = line[12:].strip()
                    elif line.startswith("GENE"):
                        current_section = "GENE"
                        gene_part = line[12:].strip()
                        if gene_part:
                            parts = gene_part.split(";")
                            if len(parts) >= 1:
                                gene_id = parts[0].strip().split()[0] if parts[0].strip() else None
                                gene_name = parts[0].strip().split()[1] if len(parts[0].strip().split()) > 1 else None
                                pathway_info["genes"].append({"id": gene_id, "symbol": gene_name})
                    elif line.startswith("DISEASE"):
                        current_section = "DISEASE"
                    elif line.startswith("DRUG"):
                        current_section = "DRUG"
                    elif line.startswith(" ") and current_section == "GENE":
                        gene_part = line.strip()
                        if gene_part:
                            parts = gene_part.split(";")
                            if len(parts) >= 1:
                                gene_split = parts[0].strip().split()
                                if len(gene_split) >= 2:
                                    pathway_info["genes"].append({
                                        "id": gene_split[0],
                                        "symbol": gene_split[1] if len(gene_split) > 1 else None
                                    })
                
                # Limit genes to avoid huge output
                pathway_info["genes"] = pathway_info["genes"][:30]
                pathway_info["n_genes"] = len(pathway_info["genes"])
                
                result = json.dumps(pathway_info, indent=2)
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"KEGG pathway fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def get_gene_pathways(gene: str) -> str:
    """Find KEGG pathways for a gene."""
    cache_key = f"kegg_gene_pathways:{gene}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            # First find the KEGG gene ID
            find_url = f"{KEGG_API}/find/genes/{gene}"
            
            async with session.get(find_url) as response:
                if response.status != 200:
                    return json.dumps({"error": f"KEGG API error: {response.status}"}, indent=2)
                
                text = await response.text()
                
                # Find human gene
                kegg_gene_id = None
                for line in text.strip().split("\n"):
                    if line and "hsa:" in line:
                        kegg_gene_id = line.split("\t")[0]
                        break
                
                if not kegg_gene_id:
                    return json.dumps({
                        "gene": gene,
                        "n_pathways": 0,
                        "message": "Gene not found in KEGG"
                    }, indent=2)
                
                # Get pathways for this gene
                link_url = f"{KEGG_API}/link/pathway/{kegg_gene_id}"
                
                async with session.get(link_url) as link_response:
                    if link_response.status != 200:
                        return json.dumps({"error": "Could not fetch pathways"}, indent=2)
                    
                    link_text = await link_response.text()
                    
                    pathway_ids = []
                    for line in link_text.strip().split("\n"):
                        if line and "\t" in line:
                            pathway_id = line.split("\t")[1].replace("path:", "")
                            pathway_ids.append(pathway_id)
                    
                    # Get pathway names
                    pathways = []
                    for pid in pathway_ids[:15]:
                        try:
                            info_url = f"{KEGG_API}/get/{pid}"
                            async with session.get(info_url) as info_response:
                                if info_response.status == 200:
                                    info_text = await info_response.text()
                                    for line in info_text.split("\n"):
                                        if line.startswith("NAME"):
                                            pathways.append({
                                                "pathway_id": pid,
                                                "name": line[12:].strip(),
                                                "url": f"https://www.kegg.jp/pathway/{pid}"
                                            })
                                            break
                        except:
                            pathways.append({"pathway_id": pid})
                    
                    result = json.dumps({
                        "gene": gene,
                        "kegg_gene_id": kegg_gene_id,
                        "n_pathways": len(pathways),
                        "pathways": pathways,
                        "source": "KEGG"
                    }, indent=2)
                    
                    set_cached(cache_key, result)
                    return result
                    
    except Exception as e:
        logger.exception(f"KEGG gene pathways fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def search_pharmgkb(query: str, resource_type: str = "all") -> str:
    """Search PharmGKB for pharmacogenomics data."""
    cache_key = f"pharmgkb_search:{query}:{resource_type}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            url = f"{PHARMGKB_API}/search/autocomplete"
            params = {"term": query}
            
            async with session.get(url, params=params) as response:
                if response.status != 200:
                    # Try alternative approach
                    return await _pharmgkb_gene_search(session, query)
                
                data = await response.json()
                
                results = []
                for item in data.get("data", []):
                    item_type = item.get("type", "").lower()
                    
                    if resource_type == "all" or resource_type in item_type:
                        results.append({
                            "id": item.get("id"),
                            "name": item.get("name"),
                            "type": item.get("type"),
                            "symbol": item.get("symbol")
                        })
                
                result = json.dumps({
                    "query": query,
                    "resource_type": resource_type,
                    "n_results": len(results),
                    "results": results[:20],
                    "source": "PharmGKB"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"PharmGKB search failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def _pharmgkb_gene_search(session, query: str) -> str:
    """Fallback PharmGKB gene search."""
    try:
        url = f"{PHARMGKB_API}/data/gene"
        params = {"symbol": query}
        
        async with session.get(url, params=params) as response:
            if response.status == 200:
                data = await response.json()
                if data.get("data"):
                    gene_data = data["data"]
                    return json.dumps({
                        "query": query,
                        "id": gene_data.get("id"),
                        "symbol": gene_data.get("symbol"),
                        "name": gene_data.get("name"),
                        "has_clinical_annotation": gene_data.get("hasClinicalAnnotation", False),
                        "has_cpic": gene_data.get("hasCpic", False),
                        "source": "PharmGKB"
                    }, indent=2)
    except:
        pass
    return json.dumps({"query": query, "n_results": 0, "message": "No results found"}, indent=2)


async def get_drug_gene_interactions(gene: str) -> str:
    """Get drug-gene interaction data from PharmGKB."""
    cache_key = f"pharmgkb_interactions:{gene}"
    cached = get_cached(cache_key)
    if cached:
        return cached
    
    try:
        async with aiohttp.ClientSession() as session:
            # Get gene info first
            gene_url = f"{PHARMGKB_API}/data/gene"
            params = {"symbol": gene}
            
            async with session.get(gene_url, params=params) as response:
                if response.status != 200:
                    return json.dumps({
                        "gene": gene,
                        "error": f"PharmGKB API error: {response.status}"
                    }, indent=2)
                
                gene_data = await response.json()
                
                if not gene_data.get("data"):
                    return json.dumps({
                        "gene": gene,
                        "error": "Gene not found in PharmGKB"
                    }, indent=2)
                
                gene_info = gene_data["data"]
                gene_id = gene_info.get("id")
                
                # Get clinical annotations
                annotations = []
                guidelines = []
                
                # Try to get related drugs
                try:
                    drugs_url = f"{PHARMGKB_API}/data/gene/{gene_id}/relatedDrugs"
                    async with session.get(drugs_url) as drugs_response:
                        if drugs_response.status == 200:
                            drugs_data = await drugs_response.json()
                            for drug in drugs_data.get("data", [])[:15]:
                                annotations.append({
                                    "drug": drug.get("name"),
                                    "drug_id": drug.get("id")
                                })
                except:
                    pass
                
                result = json.dumps({
                    "gene": gene,
                    "pharmgkb_id": gene_id,
                    "name": gene_info.get("name"),
                    "has_clinical_annotation": gene_info.get("hasClinicalAnnotation", False),
                    "has_cpic_guideline": gene_info.get("hasCpic", False),
                    "has_fda_label": gene_info.get("hasFdaLabel", False),
                    "n_related_drugs": len(annotations),
                    "related_drugs": annotations,
                    "pharmgkb_url": f"https://www.pharmgkb.org/gene/{gene_id}",
                    "source": "PharmGKB"
                }, indent=2)
                
                set_cached(cache_key, result)
                return result
                
    except Exception as e:
        logger.exception(f"PharmGKB interaction fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


def register_structure_tools():
    """Return structure/pathway/pharma tools for registration."""
    return STRUCTURE_TOOLS
