# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-02-09

### Added
- **Protein & Gene Tools**
  - UniProt search and protein details
  - NCBI Gene database search
  - Ensembl gene and variant lookup
  - InterPro domain information

- **Clinical & Variant Tools**
  - ClinVar clinical variant search and details
  - GWAS Catalog query
  - GTEx eQTL data

- **Protein Interaction Tools**
  - STRING protein-protein interactions
  - Protein interaction networks
  - Functional enrichment analysis

- **Structure & Pathway Tools**
  - AlphaFold AI-predicted structure lookup
  - PDB experimental structure search
  - KEGG pathway search and details
  - Gene-to-pathway mapping

- **Drug Discovery Tools**
  - Open Targets drug target information
  - Open Targets disease associations
  - PharmGKB drug-gene interactions

- **Genetic Disease Tools**
  - OMIM genetic disease search
  - Combined gene-disease lookup

- **Performance Features**
  - Smart caching with 1-hour TTL
  - Async API operations
  - Graceful error handling

### Infrastructure
- Full MCP server implementation
- Claude Desktop integration
- PyPI package distribution ready
