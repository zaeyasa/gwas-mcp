
import sys
import os
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path("src").resolve()))

def main():
    print("Testing VCF Fallback Parser...")
    
    try:
        import pysam
        print("WARNING: pysam IS installed! Fallback not being tested.")
    except ImportError:
        print("CONFIRMED: pysam is NOT installed. Testing fallback...")
    
    from gwas_mcp.utils.file_handlers import read_vcf
    
    # Create valid dummy VCF
    with open("test_fallback.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("1\t100\trs1\tA\tG\t.\tPASS\t.\n")
        f.write("1\t200\trs2\tC\tT\t.\tPASS\t.\n")
        
    try:
        df = read_vcf("test_fallback.vcf")
        print("\nSUCCESS! Parsed VCF without pysam.")
        print("Data:")
        print(df)
        
        if len(df) == 2 and df.iloc[0]['ID'] == 'rs1':
             print("\nVerification: Data content is correct.")
        else:
             print("\nVerification: Data content is INCORRECT.")
             
    except Exception as e:
        print(f"\nFAILURE: {e}")
    finally:
        if os.path.exists("test_fallback.vcf"):
            os.remove("test_fallback.vcf")

if __name__ == "__main__":
    main()
