"""
GWAS MCP Server - Main Entry Point

A professional MCP server for GWAS and bioinformatics research,
compatible with Claude Desktop and other MCP clients.
"""

import asyncio
import logging
import os
import sys
from typing import Any
from pathlib import Path

from dotenv import load_dotenv
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import TextContent, Tool, Resource, TextResourceContents

# Load environment variables
load_dotenv()

# Configure logging
logging.basicConfig(
    level=getattr(logging, os.getenv('LOG_LEVEL', 'INFO')),
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stderr),
    ]
)
logger = logging.getLogger(__name__)

# Create the MCP server instance
server = Server("gwas-bioinformatics")

# Import and register tools
from gwas_mcp.tools.qc_tools import QC_TOOLS, handle_qc_tool
from gwas_mcp.tools.analysis_tools import ANALYSIS_TOOLS, handle_analysis_tool
from gwas_mcp.tools.annotation_tools import ANNOTATION_TOOLS, handle_annotation_tool
from gwas_mcp.tools.visualization_tools import VISUALIZATION_TOOLS, handle_visualization_tool
from gwas_mcp.tools.post_gwas_tools import POST_GWAS_TOOLS, handle_post_gwas_tool
from gwas_mcp.tools.protein_tools import PROTEIN_TOOLS, handle_protein_tool
from gwas_mcp.tools.clinical_tools import CLINICAL_TOOLS, handle_clinical_tool
from gwas_mcp.tools.structure_tools import STRUCTURE_TOOLS, handle_structure_tool
from gwas_mcp.tools.advanced_tools import ADVANCED_TOOLS, handle_advanced_tool
from gwas_mcp.resources.db_resources import RESOURCES, handle_resource


@server.list_tools()
async def list_tools() -> list[Tool]:
    """List all available GWAS tools."""
    all_tools = []
    all_tools.extend(QC_TOOLS)
    all_tools.extend(ANALYSIS_TOOLS)
    all_tools.extend(ANNOTATION_TOOLS)
    all_tools.extend(VISUALIZATION_TOOLS)
    all_tools.extend(POST_GWAS_TOOLS)
    all_tools.extend(PROTEIN_TOOLS)
    all_tools.extend(CLINICAL_TOOLS)
    all_tools.extend(STRUCTURE_TOOLS)
    all_tools.extend(ADVANCED_TOOLS)
    
    logger.info(f"Listing {len(all_tools)} available tools")
    return all_tools


@server.call_tool()
async def call_tool(name: str, arguments: dict[str, Any]) -> list[TextContent]:
    """Handle tool calls by routing to appropriate handler."""
    logger.info(f"Tool called: {name} with args: {arguments}")
    
    try:
        # Route to appropriate handler based on tool name
        if name in [t.name for t in QC_TOOLS]:
            result = await handle_qc_tool(name, arguments)
        elif name in [t.name for t in ANALYSIS_TOOLS]:
            result = await handle_analysis_tool(name, arguments)
        elif name in [t.name for t in ANNOTATION_TOOLS]:
            result = await handle_annotation_tool(name, arguments)
        elif name in [t.name for t in VISUALIZATION_TOOLS]:
            result = await handle_visualization_tool(name, arguments)
        elif name in [t.name for t in POST_GWAS_TOOLS]:
            result = await handle_post_gwas_tool(name, arguments)
        elif name in [t.name for t in PROTEIN_TOOLS]:
            result = await handle_protein_tool(name, arguments)
        elif name in [t.name for t in CLINICAL_TOOLS]:
            result = await handle_clinical_tool(name, arguments)
        elif name in [t.name for t in STRUCTURE_TOOLS]:
            result = await handle_structure_tool(name, arguments)
        elif name in [t.name for t in ADVANCED_TOOLS]:
            result = await handle_advanced_tool(name, arguments)
        else:
            raise ValueError(f"Unknown tool: {name}")
        
        return [TextContent(type="text", text=result)]
        
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        return [TextContent(type="text", text=f"Error: File not found - {e}")]
    except ValueError as e:
        logger.error(f"Validation error: {e}")
        return [TextContent(type="text", text=f"Error: Invalid input - {e}")]
    except Exception as e:
        logger.exception(f"Tool execution failed: {e}")
        return [TextContent(type="text", text=f"Error: {type(e).__name__} - {e}")]


@server.list_resources()
async def list_resources() -> list[Resource]:
    """List available database resources."""
    logger.info(f"Listing {len(RESOURCES)} available resources")
    return RESOURCES


@server.read_resource()
async def read_resource(uri: str) -> str:
    """Read data from a resource URI."""
    logger.info(f"Reading resource: {uri}")
    
    try:
        result = await handle_resource(str(uri))
        return result
    except Exception as e:
        logger.exception(f"Resource read failed: {e}")
        raise


async def run_server():
    """Run the MCP server using stdio transport."""
    logger.info("Starting GWAS MCP Server...")
    
    async with stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            server.create_initialization_options()
        )


def main():
    """Main entry point."""
    try:
        # Windows-specific fix for asyncio pipes
        if sys.platform == 'win32':
             asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
             
        asyncio.run(run_server())
    except KeyboardInterrupt:
        logger.info("Server shutdown requested")
    except Exception as e:
        logger.exception(f"Server error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
