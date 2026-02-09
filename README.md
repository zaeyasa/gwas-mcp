# 🧬 GWAS-MCP: Bioinformatics MCP Server

[![PyPI version](https://badge.fury.io/py/gwas-mcp.svg)](https://badge.fury.io/py/gwas-mcp)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MCP](https://img.shields.io/badge/MCP-Compatible-blue)](https://modelcontextprotocol.io)

A powerful **Model Context Protocol (MCP)** server for GWAS and bioinformatics research. Seamlessly integrates with Claude Desktop and other MCP clients to provide AI-powered access to major biological databases.

<p align="center">
  <img src="https://img.shields.io/badge/Tools-30+-green" alt="30+ Tools">
  <img src="https://img.shields.io/badge/Databases-12+-blue" alt="12+ Databases">
  <img src="https://img.shields.io/badge/Python-3.10+-blue" alt="Python 3.10+">
</p>

---

## ✨ Features

### 🔬 Protein & Gene Lookup
- **UniProt** - Search proteins by name, gene, or ID
- **Ensembl** - Gene information and variant details
- **NCBI Gene** - Comprehensive gene database

### 🧪 Clinical & Variants
- **ClinVar** - Clinical variant interpretations (pathogenic/benign)
- **GWAS Catalog** - Genome-wide association studies
- **GTEx** - Expression quantitative trait loci (eQTL)

### 🔗 Protein Interactions & Networks
- **STRING** - Protein-protein interactions
- **InterPro** - Protein domains and families

### 🏗️ Structures & Pathways
- **AlphaFold** - AI-predicted protein structures
- **PDB** - Experimental 3D structures
- **KEGG** - Metabolic and signaling pathways

### 💊 Drug Discovery
- **Open Targets** - Drug target validation & disease associations
- **PharmGKB** - Pharmacogenomics & drug-gene interactions

### 🏥 Genetic Diseases
- **OMIM** - Online Mendelian Inheritance in Man

---

## 🚀 Quick Start

### Installation

```bash
pip install gwas-mcp
```

### Claude Desktop Configuration

Add to your `claude_desktop_config.json`:

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

**Config file location:**
- Windows: `%APPDATA%\Claude\claude_desktop_config.json`
- macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
- Linux: `~/.config/Claude/claude_desktop_config.json`

### Restart Claude Desktop

After adding the configuration, restart Claude Desktop to load the MCP server.

---

## 🛠️ Available Tools

### Protein & Gene Tools

| Tool | Description |
|------|-------------|
| `search_uniprot` | Search UniProt by protein name, gene, or ID |
| `get_protein_details` | Get detailed protein info (function, domains, GO terms) |
| `search_ncbi_gene` | Search NCBI Gene database |
| `search_ensembl_gene` | Get gene location and details from Ensembl |
| `get_variant_info` | Get SNP/variant info by rsID |
| `get_interpro_domains` | Get protein domain information |

### Clinical & Variant Tools

| Tool | Description |
|------|-------------|
| `search_clinvar` | Search ClinVar for clinical variants |
| `get_clinvar_variant` | Get clinical interpretation for a variant |
| `annotate_snps` | Annotate SNPs with functional consequences |
| `query_gwas_catalog` | Query GWAS Catalog for associations |
| `get_eqtl_data` | Get eQTL data from GTEx |

### Protein Interaction Tools

| Tool | Description |
|------|-------------|
| `get_protein_interactions` | Find interacting proteins (STRING) |
| `get_interaction_network` | Get network between multiple proteins |
| `get_functional_enrichment` | Pathway/GO enrichment analysis |

### Structure & Pathway Tools

| Tool | Description |
|------|-------------|
| `get_alphafold_structure` | Get AI-predicted structure |
| `search_alphafold` | Search AlphaFold database |
| `search_pdb_structures` | Search PDB for 3D structures |
| `get_pdb_structure` | Get PDB structure details |
| `search_kegg_pathway` | Search KEGG pathways |
| `get_kegg_pathway` | Get pathway genes and details |
| `get_gene_pathways` | Find pathways for a gene |

### Drug Discovery Tools

| Tool | Description |
|------|-------------|
| `get_drug_targets` | Find drugs targeting a gene (Open Targets) |
| `get_disease_associations` | Get disease associations with scores |
| `search_open_targets` | Search genes, diseases, or drugs |
| `search_pharmgkb` | Search PharmGKB database |
| `get_drug_gene_interactions` | Get drug-gene interactions |

### Genetic Disease Tools

| Tool | Description |
|------|-------------|
| `search_omim` | Search OMIM for genetic diseases |
| `get_gene_diseases` | Get all diseases for a gene |

---

## 💬 Example Prompts

Once configured, ask Claude naturally:

### Protein & Gene Queries
> "Get information about the BRCA1 gene"
> 
> "Search UniProt for hemoglobin"
> 
> "What protein has UniProt ID P53_HUMAN?"

### Clinical Variants
> "Is the BRCA1 variant rs80357906 pathogenic?"
> 
> "Search ClinVar for TP53 variants"

### Protein Interactions
> "What proteins interact with TP53?"
> 
> "Find functional enrichment for BRCA1, ATM, and CHEK2"

### Structures & Pathways
> "Get the AlphaFold structure for TP53"
> 
> "What pathways is BRCA1 involved in?"
> 
> "Search PDB for insulin structures"

### Drug Discovery
> "What drugs target EGFR?"
> 
> "What diseases is BRAF associated with?"

### Genetic Diseases
> "Search OMIM for cystic fibrosis"
> 
> "What diseases are linked to the CFTR gene?"

---

## ⚡ Performance Features

- **Smart Caching** - API responses cached for 1 hour to improve speed
- **Async Operations** - All API calls are non-blocking
- **Error Handling** - Graceful handling of API failures

---

## 🔧 Development

### From Source

```bash
# Clone the repository
git clone https://github.com/zaeyasa/gwas-mcp.git
cd gwas-mcp

# Install dependencies
pip install -e .

# Run the server
python -m gwas_mcp.server
```

### Project Structure

```
gwas-mcp/
├── src/
│   └── gwas_mcp/
│       ├── server.py           # Main MCP server
│       ├── tools/
│       │   ├── protein_tools.py     # UniProt, NCBI, Ensembl
│       │   ├── clinical_tools.py    # ClinVar, STRING
│       │   ├── structure_tools.py   # PDB, KEGG, PharmGKB
│       │   └── advanced_tools.py    # AlphaFold, Open Targets, OMIM
│       └── resources/
│           └── db_resources.py      # Database resources
├── pyproject.toml
├── README.md
└── LICENSE
```

---

## 📊 Supported Databases

| Database | Type | Description |
|----------|------|-------------|
| [UniProt](https://www.uniprot.org/) | Protein | Protein sequences and annotations |
| [Ensembl](https://www.ensembl.org/) | Gene/Variant | Genome browser and variant data |
| [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/) | Gene | Gene information database |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Clinical | Clinical variant interpretations |
| [GWAS Catalog](https://www.ebi.ac.uk/gwas/) | GWAS | Genome-wide association studies |
| [GTEx](https://gtexportal.org/) | Expression | Expression QTL data |
| [STRING](https://string-db.org/) | Interactions | Protein-protein interactions |
| [InterPro](https://www.ebi.ac.uk/interpro/) | Domains | Protein families and domains |
| [AlphaFold](https://alphafold.ebi.ac.uk/) | Structure | AI-predicted structures |
| [PDB](https://www.rcsb.org/) | Structure | Experimental 3D structures |
| [KEGG](https://www.kegg.jp/) | Pathways | Metabolic and signaling pathways |
| [Open Targets](https://platform.opentargets.org/) | Drug Discovery | Drug targets and disease associations |
| [PharmGKB](https://www.pharmgkb.org/) | Pharmacogenomics | Drug-gene interactions |
| [OMIM](https://omim.org/) | Diseases | Genetic disease database |

---

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- [Model Context Protocol](https://modelcontextprotocol.io/) - The MCP specification
- [Anthropic](https://www.anthropic.com/) - Claude AI and MCP development
- All the amazing bioinformatics databases that make this possible

---

## 📬 Contact

- GitHub: [@zaeyasa](https://github.com/zaeyasa)

---

<p align="center">Made with ❤️ for the bioinformatics community</p>
