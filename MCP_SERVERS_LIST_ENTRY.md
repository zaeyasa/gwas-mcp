# MCP Servers List Entry

This file contains the entry for submitting to the official MCP Servers List.

## Submission Format

For the [Anthropic MCP Servers List](https://github.com/modelcontextprotocol/servers):

### Entry

```markdown
### gwas-mcp

Bioinformatics MCP server providing access to 30+ tools from 14 major biological databases including UniProt, ClinVar, AlphaFold, KEGG, Open Targets, and more. Perfect for GWAS research, variant analysis, drug discovery, and protein structure exploration.

**Features:**
- Protein & gene information (UniProt, Ensembl, NCBI)
- Clinical variants (ClinVar, GWAS Catalog)
- Protein interactions (STRING)
- AI structures (AlphaFold) & experimental structures (PDB)
- Pathways (KEGG)
- Drug targets (Open Targets, PharmGKB)
- Genetic diseases (OMIM)

[GitHub](https://github.com/zaeyasa/gwas-mcp) | [PyPI](https://pypi.org/project/gwas-mcp/)
```

## Categories

This server fits in the following categories:
- **Science & Research**
- **Healthcare & Medical**
- **Data & APIs**

## Installation Command

```bash
pip install gwas-mcp
```

## Configuration

```json
{
  "mcpServers": {
    "gwas-bioinformatics": {
      "command": "python",
      "args": ["-m", "gwas_mcp.server"]
    }
  }
}
```
