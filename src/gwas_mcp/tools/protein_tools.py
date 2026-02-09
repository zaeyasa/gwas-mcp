"""
Protein Lookup Tools for GWAS MCP Server.

Tools for searching proteins, genes, and domains via UniProt, NCBI, and InterPro APIs.
"""

import json
import logging
from typing import Any, Dict, List, Optional

import aiohttp

from mcp.types import Tool

logger = logging.getLogger(__name__)

# API endpoints
UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
NCBI_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"


# Define Protein Lookup tools
PROTEIN_TOOLS = [
    Tool(
        name="search_uniprot",
        description="Search UniProt for protein information by protein name, gene name, or UniProt ID. Returns protein function, sequence info, and associated genes.",
        inputSchema={
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Search query - can be protein name (e.g., 'hemoglobin'), gene symbol (e.g., 'BRCA1'), or UniProt ID (e.g., 'P38398')"
                },
                "organism": {
                    "type": "string",
                    "description": "Organism to filter results (default: 'human')",
                    "default": "human"
                },
                "limit": {
                    "type": "integer",
                    "description": "Maximum number of results to return (default: 10)",
                    "default": 10
                }
            },
            "required": ["query"]
        }
    ),
    Tool(
        name="get_protein_details",
        description="Get detailed protein information from UniProt by UniProt accession ID. Returns full protein details including function, domains, and GO annotations.",
        inputSchema={
            "type": "object",
            "properties": {
                "uniprot_id": {
                    "type": "string",
                    "description": "UniProt accession ID (e.g., 'P38398' for BRCA1)"
                }
            },
            "required": ["uniprot_id"]
        }
    ),
    Tool(
        name="search_ncbi_gene",
        description="Search NCBI Gene database for gene information by gene symbol, name, or ID.",
        inputSchema={
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Gene symbol (e.g., 'TP53'), name, or NCBI Gene ID"
                },
                "organism": {
                    "type": "string",
                    "description": "Organism to filter results (default: 'Homo sapiens')",
                    "default": "Homo sapiens"
                }
            },
            "required": ["query"]
        }
    ),
    Tool(
        name="get_interpro_domains",
        description="Get protein domain information from InterPro for a given UniProt ID or protein sequence.",
        inputSchema={
            "type": "object",
            "properties": {
                "uniprot_id": {
                    "type": "string",
                    "description": "UniProt accession ID to look up domains for"
                }
            },
            "required": ["uniprot_id"]
        }
    ),
    Tool(
        name="search_ensembl_gene",
        description="Search Ensembl for gene information by gene symbol or Ensembl ID. Returns gene location, biotype, and description.",
        inputSchema={
            "type": "object",
            "properties": {
                "gene": {
                    "type": "string",
                    "description": "Gene symbol (e.g., 'BRCA1') or Ensembl ID (e.g., 'ENSG00000012048')"
                },
                "species": {
                    "type": "string",
                    "description": "Species (default: 'homo_sapiens')",
                    "default": "homo_sapiens"
                }
            },
            "required": ["gene"]
        }
    ),
    Tool(
        name="get_variant_info",
        description="Get detailed variant/SNP information from Ensembl by rsID. Returns position, alleles, clinical significance, and consequences.",
        inputSchema={
            "type": "object",
            "properties": {
                "rsid": {
                    "type": "string",
                    "description": "dbSNP rsID (e.g., 'rs1234567')"
                }
            },
            "required": ["rsid"]
        }
    ),
]


async def handle_protein_tool(name: str, arguments: Dict[str, Any]) -> str:
    """Handle protein lookup tool calls."""
    
    if name == "search_uniprot":
        return await search_uniprot(
            query=arguments["query"],
            organism=arguments.get("organism", "human"),
            limit=arguments.get("limit", 10)
        )
    elif name == "get_protein_details":
        return await get_protein_details(arguments["uniprot_id"])
    elif name == "search_ncbi_gene":
        return await search_ncbi_gene(
            query=arguments["query"],
            organism=arguments.get("organism", "Homo sapiens")
        )
    elif name == "get_interpro_domains":
        return await get_interpro_domains(arguments["uniprot_id"])
    elif name == "search_ensembl_gene":
        return await search_ensembl_gene(
            gene=arguments["gene"],
            species=arguments.get("species", "homo_sapiens")
        )
    elif name == "get_variant_info":
        return await get_variant_info(arguments["rsid"])
    else:
        raise ValueError(f"Unknown protein tool: {name}")


async def search_uniprot(query: str, organism: str = "human", limit: int = 10) -> str:
    """
    Search UniProt for proteins.
    """
    try:
        async with aiohttp.ClientSession() as session:
            # Build search query
            search_query = f'({query}) AND (organism_name:{organism})'
            url = f"{UNIPROT_API}/search"
            params = {
                "query": search_query,
                "format": "json",
                "size": limit,
                "fields": "accession,id,gene_names,protein_name,organism_name,length,cc_function"
            }
            
            async with session.get(url, params=params) as response:
                if response.status == 200:
                    data = await response.json()
                    
                    results = []
                    for entry in data.get("results", []):
                        protein_info = {
                            "uniprot_id": entry.get("primaryAccession"),
                            "entry_name": entry.get("uniProtkbId"),
                            "protein_name": entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "N/A"),
                            "gene_names": [g.get("geneName", {}).get("value") for g in entry.get("genes", []) if g.get("geneName")],
                            "organism": entry.get("organism", {}).get("scientificName"),
                            "length": entry.get("sequence", {}).get("length"),
                            "function": None
                        }
                        
                        # Get function if available
                        for comment in entry.get("comments", []):
                            if comment.get("commentType") == "FUNCTION":
                                texts = comment.get("texts", [])
                                if texts:
                                    protein_info["function"] = texts[0].get("value")
                        
                        results.append(protein_info)
                    
                    return json.dumps({
                        "query": query,
                        "organism": organism,
                        "n_results": len(results),
                        "results": results
                    }, indent=2)
                else:
                    return json.dumps({
                        "error": f"UniProt API returned status {response.status}"
                    }, indent=2)
                    
    except Exception as e:
        logger.exception(f"UniProt search failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def get_protein_details(uniprot_id: str) -> str:
    """
    Get detailed protein information from UniProt.
    """
    try:
        async with aiohttp.ClientSession() as session:
            url = f"{UNIPROT_API}/{uniprot_id}.json"
            
            async with session.get(url) as response:
                if response.status == 200:
                    data = await response.json()
                    
                    # Extract key information
                    protein_info = {
                        "uniprot_id": data.get("primaryAccession"),
                        "entry_name": data.get("uniProtkbId"),
                        "protein_name": data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"),
                        "gene_names": [g.get("geneName", {}).get("value") for g in data.get("genes", []) if g.get("geneName")],
                        "organism": data.get("organism", {}).get("scientificName"),
                        "sequence_length": data.get("sequence", {}).get("length"),
                        "function": None,
                        "subcellular_location": [],
                        "go_terms": [],
                        "domains": []
                    }
                    
                    # Parse comments
                    for comment in data.get("comments", []):
                        if comment.get("commentType") == "FUNCTION":
                            texts = comment.get("texts", [])
                            if texts:
                                protein_info["function"] = texts[0].get("value")
                        elif comment.get("commentType") == "SUBCELLULAR LOCATION":
                            for loc in comment.get("subcellularLocations", []):
                                location = loc.get("location", {}).get("value")
                                if location:
                                    protein_info["subcellular_location"].append(location)
                    
                    # Parse cross-references for GO terms
                    for xref in data.get("uniProtKBCrossReferences", []):
                        if xref.get("database") == "GO":
                            go_id = xref.get("id")
                            props = {p.get("key"): p.get("value") for p in xref.get("properties", [])}
                            protein_info["go_terms"].append({
                                "id": go_id,
                                "term": props.get("GoTerm", "").split(":")[-1] if props.get("GoTerm") else None,
                                "category": props.get("GoTerm", "").split(":")[0] if props.get("GoTerm") else None
                            })
                        elif xref.get("database") == "InterPro":
                            props = {p.get("key"): p.get("value") for p in xref.get("properties", [])}
                            protein_info["domains"].append({
                                "id": xref.get("id"),
                                "name": props.get("EntryName")
                            })
                    
                    # Limit GO terms to avoid huge output
                    protein_info["go_terms"] = protein_info["go_terms"][:15]
                    
                    return json.dumps(protein_info, indent=2)
                elif response.status == 404:
                    return json.dumps({
                        "uniprot_id": uniprot_id,
                        "error": "Protein not found"
                    }, indent=2)
                else:
                    return json.dumps({
                        "error": f"UniProt API returned status {response.status}"
                    }, indent=2)
                    
    except Exception as e:
        logger.exception(f"UniProt details fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def search_ncbi_gene(query: str, organism: str = "Homo sapiens") -> str:
    """
    Search NCBI Gene database.
    """
    try:
        async with aiohttp.ClientSession() as session:
            # First search for gene IDs
            search_url = f"{NCBI_API}/esearch.fcgi"
            search_params = {
                "db": "gene",
                "term": f'{query}[Gene Name] AND "{organism}"[Organism]',
                "retmode": "json",
                "retmax": 5
            }
            
            async with session.get(search_url, params=search_params) as search_response:
                if search_response.status != 200:
                    return json.dumps({"error": f"NCBI search failed with status {search_response.status}"}, indent=2)
                
                search_data = await search_response.json()
                gene_ids = search_data.get("esearchresult", {}).get("idlist", [])
                
                if not gene_ids:
                    return json.dumps({
                        "query": query,
                        "organism": organism,
                        "n_results": 0,
                        "message": "No genes found"
                    }, indent=2)
                
                # Fetch gene summaries
                summary_url = f"{NCBI_API}/esummary.fcgi"
                summary_params = {
                    "db": "gene",
                    "id": ",".join(gene_ids),
                    "retmode": "json"
                }
                
                async with session.get(summary_url, params=summary_params) as summary_response:
                    if summary_response.status != 200:
                        return json.dumps({"error": f"NCBI summary failed with status {summary_response.status}"}, indent=2)
                    
                    summary_data = await summary_response.json()
                    
                    results = []
                    for gene_id in gene_ids:
                        gene_info = summary_data.get("result", {}).get(gene_id, {})
                        if gene_info:
                            results.append({
                                "ncbi_gene_id": gene_id,
                                "symbol": gene_info.get("name"),
                                "description": gene_info.get("description"),
                                "organism": gene_info.get("organism", {}).get("scientificname"),
                                "chromosome": gene_info.get("chromosome"),
                                "aliases": gene_info.get("otheraliases", "").split(", ") if gene_info.get("otheraliases") else [],
                                "summary": gene_info.get("summary", "")[:500] + "..." if len(gene_info.get("summary", "")) > 500 else gene_info.get("summary", "")
                            })
                    
                    return json.dumps({
                        "query": query,
                        "organism": organism,
                        "n_results": len(results),
                        "results": results
                    }, indent=2)
                    
    except Exception as e:
        logger.exception(f"NCBI Gene search failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def get_interpro_domains(uniprot_id: str) -> str:
    """
    Get protein domain information from InterPro.
    """
    try:
        async with aiohttp.ClientSession() as session:
            url = f"{INTERPRO_API}/protein/UniProt/{uniprot_id}"
            headers = {"Accept": "application/json"}
            
            async with session.get(url, headers=headers) as response:
                if response.status == 200:
                    data = await response.json()
                    
                    domains = []
                    for entry in data.get("results", []):
                        metadata = entry.get("metadata", {})
                        domains.append({
                            "interpro_id": metadata.get("accession"),
                            "name": metadata.get("name"),
                            "type": metadata.get("type"),
                            "description": metadata.get("description", {}).get("p", [{}])[0].get("text") if metadata.get("description") else None,
                            "source_database": metadata.get("source_database"),
                            "go_terms": [{"id": go.get("identifier"), "name": go.get("name")} for go in metadata.get("go_terms", [])][:5]
                        })
                    
                    return json.dumps({
                        "uniprot_id": uniprot_id,
                        "n_domains": len(domains),
                        "domains": domains
                    }, indent=2)
                elif response.status == 404:
                    return json.dumps({
                        "uniprot_id": uniprot_id,
                        "n_domains": 0,
                        "message": "No domain information found"
                    }, indent=2)
                else:
                    return json.dumps({
                        "error": f"InterPro API returned status {response.status}"
                    }, indent=2)
                    
    except Exception as e:
        logger.exception(f"InterPro fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def search_ensembl_gene(gene: str, species: str = "homo_sapiens") -> str:
    """
    Search Ensembl for gene information.
    """
    try:
        async with aiohttp.ClientSession() as session:
            url = f"https://rest.ensembl.org/lookup/symbol/{species}/{gene}"
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
                        "species": data.get("species"),
                        "source": "Ensembl"
                    }
                    
                    return json.dumps(gene_info, indent=2)
                elif response.status == 404:
                    return json.dumps({
                        "gene": gene,
                        "error": "Gene not found in Ensembl"
                    }, indent=2)
                else:
                    return json.dumps({
                        "error": f"Ensembl API returned status {response.status}"
                    }, indent=2)
                    
    except Exception as e:
        logger.exception(f"Ensembl gene search failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


async def get_variant_info(rsid: str) -> str:
    """
    Get variant information from Ensembl.
    """
    try:
        async with aiohttp.ClientSession() as session:
            url = f"https://rest.ensembl.org/variation/human/{rsid}"
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
                        "error": f"Ensembl API returned status {response.status}"
                    }, indent=2)
                    
    except Exception as e:
        logger.exception(f"Ensembl variant fetch failed: {e}")
        return json.dumps({"error": str(e)}, indent=2)


def register_protein_tools():
    """Return protein tools for registration."""
    return PROTEIN_TOOLS
