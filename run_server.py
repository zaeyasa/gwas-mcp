#!/usr/bin/env python
"""
Wrapper script to run the GWAS MCP Server.
This script ensures the correct Python path is set before importing the server.
"""
import sys
import os

# Add the src directory to the Python path
src_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src")
sys.path.insert(0, src_dir)

# Now import and run the server
from gwas_mcp.server import main

if __name__ == "__main__":
    main()
