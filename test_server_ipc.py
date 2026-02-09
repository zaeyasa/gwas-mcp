
import subprocess
import sys
import json
import os
from pathlib import Path

# Path to the source directory
SRC_DIR = Path("I:/ZZZ_Projects/GWAS-MCP/src")

def test_server_resource_read():
    print("Starting server process...")
    
    # Run the server as a subprocess
    process = subprocess.Popen(
        [sys.executable, "-m", "gwas_mcp.server"],
        cwd=str(SRC_DIR),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=sys.stderr, # Let stderr flow through to console
        text=True,
        env={**os.environ, "PYTHONPATH": str(SRC_DIR)}
    )
    
    try:
        # 1. Initialize
        print("\nSending initialize request...")
        init_req = {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "initialize",
            "params": {
                "protocolVersion": "2024-11-05",
                "capabilities": {},
                "clientInfo": {"name": "test-script", "version": "1.0"}
            }
        }
        process.stdin.write(json.dumps(init_req) + "\n")
        process.stdin.flush()
        
        # Read init response
        response = process.stdout.readline()
        print(f"Init Response: {response.strip()}")
        
        # 2. Send initialized notification
        process.stdin.write(json.dumps({
            "jsonrpc": "2.0",
            "method": "notifications/initialized"
        }) + "\n")
        process.stdin.flush()
        
        # 3. Request logic: Read a resource
        # We'll try the 'GWAS Catalog Traits' resource which is static and should work
        print("\nSending resources/read request...")
        read_req = {
            "jsonrpc": "2.0",
            "id": 2,
            "method": "resources/read",
            "params": {
                "uri": "gwas://catalog/traits"
            }
        }
        process.stdin.write(json.dumps(read_req) + "\n")
        process.stdin.flush()
        
        # Read response
        response = process.stdout.readline()
        print(f"Read Response (first 200 chars): {response.strip()[:200]}...")
        
        # Parse and check structure
        resp_json = json.loads(response)
        if "error" in resp_json:
            print(f"\n❌ Error in response: {resp_json['error']}")
        elif "result" in resp_json:
            result = resp_json["result"]
            print("\n✅ Result received!")
            print(f"Keys in result: {result.keys()}")
            
            if "contents" in result:
                contents = result["contents"]
                print(f"Number of content items: {len(contents)}")
                if len(contents) > 0:
                    item = contents[0]
                    print(f"Item mimeType: {item.get('mimeType')}")
                    print(f"Item uri: {item.get('uri')}")
                    print(f"Item text length: {len(item.get('text', ''))}")
            else:
                print("❌ 'contents' key missing in result!")
        
    except Exception as e:
        print(f"\n❌ Test failed with exception: {e}")
    finally:
        print("\nClosing server...")
        process.terminate()
        process.wait()

if __name__ == "__main__":
    test_server_resource_read()
