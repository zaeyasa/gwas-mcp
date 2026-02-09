# MCP Registry Submission Guide

## The New Process (2025+)

The MCP Registry no longer uses pull requests. Instead, you use the `mcp-publisher` CLI tool.

---

## Step-by-Step Instructions

### Step 1: Publish to PyPI First

```bash
twine upload dist/*
```

The MCP registry only stores metadata - your package must be on PyPI first.

### Step 2: Install mcp-publisher CLI

**Windows PowerShell:**
```powershell
$arch = if ([System.Runtime.InteropServices.RuntimeInformation]::ProcessArchitecture -eq "Arm64") { "arm64" } else { "amd64" }
Invoke-WebRequest -Uri "https://github.com/modelcontextprotocol/registry/releases/latest/download/mcp-publisher_windows_$arch.tar.gz" -OutFile "mcp-publisher.tar.gz"
tar xf mcp-publisher.tar.gz mcp-publisher.exe
rm mcp-publisher.tar.gz
# Move mcp-publisher.exe to a directory in your PATH (e.g., C:\Windows\System32)
```

Verify installation:
```bash
mcp-publisher --help
```

### Step 3: Authenticate with GitHub

```bash
mcp-publisher login github
```

This will:
1. Give you a URL to visit
2. Show a code to enter
3. Ask you to authorize the app

After authorization, you'll see: `✓ Successfully logged in`

### Step 4: Publish to MCP Registry

Navigate to your project folder and run:

```bash
cd I:\ZZZ_Projects\GWAS-MCP
mcp-publisher publish
```

You should see:
```
Publishing to https://registry.modelcontextprotocol.io...
✓ Successfully published
✓ Server io.github.zaeyasa/gwas-mcp version 1.0.0
```

### Step 5: Verify

Check your server is listed:
```bash
curl "https://registry.modelcontextprotocol.io/v0.1/servers?search=gwas-mcp"
```

---

## Files I Created For You

| File | Purpose |
|------|---------|
| `server.json` | Registry metadata (name, description, package info) |
| `README.md` | Updated with `<!-- mcp-name: ... -->` comment for verification |

---

## Important Notes

- Your server name MUST start with `io.github.zaeyasa/` because you're using GitHub auth
- The `mcp-name` comment in README.md is required for PyPI package verification
- You don't need to fork the registry repo anymore - the CLI handles everything!

---

## Troubleshooting

| Error | Solution |
|-------|----------|
| "Registry validation failed" | Make sure `mcp-name` comment is in README.md |
| "Invalid token" | Run `mcp-publisher login github` again |
| "Permission denied" | Server name must match `io.github.YOUR-USERNAME/...` |
