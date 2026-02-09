
import sys
import os
import asyncio
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path("src").resolve()))

from gwas_mcp.resources.db_resources import handle_resource

async def test_resources():
    print("--- Testing Static Resource (Traits) ---")
    try:
        # This one should work
        result = await handle_resource("gwas://catalog/traits")
        print(f"SUCCESS (len={len(result)})")
    except Exception as e:
        print(f"FAILED: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()

    print("\n--- Testing Templated Resource (SNP) ---")
    try:
        # This simulates clicking the menu item directly which sends the literal template
        result = await handle_resource("gwas://catalog/snp/{rsid}")
        # It won't "work" in fetching real data, but it shouldn't crash with a 500
        print(f"Result (len={len(result)}): {result[:100]}...")
    except Exception as e:
        print(f"FAILED: {type(e).__name__}: {e}")

if __name__ == "__main__":
    # if sys.platform == 'win32':
    #     asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
        
    try:
        asyncio.run(test_resources())
    except KeyboardInterrupt:
        pass
