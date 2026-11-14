"""
Generate both standard signature plots (using plot.py) and RePrint plots (using plot_static.py).

Usage:
    reprint-all-plots <cosmic_file> [output_base_name]
    
Example:
    reprint-all-plots tests/data/COSMIC_v3.4_SBS_GRCh37.txt cosmic_v3.4_GRCh37
"""

import pandas as pd
import os
import sys
from reprint_tool.core import reprint
from reprint_tool.plot import save_all_signatures_to_single_pdf
from reprint_tool.plot_static import generate_all_reprints

def main():
    if len(sys.argv) < 2:
        print("Usage: reprint-all-plots <cosmic_file> [output_base_name]")
        print("\nExample:")
        print("  reprint-all-plots tests/data/COSMIC_v3.4_SBS_GRCh37.txt cosmic_v3.4_GRCh37")
        sys.exit(1)
    
    cosmic_file = sys.argv[1]
    
    # Determine base name for output files
    if len(sys.argv) > 2:
        base_name = sys.argv[2]
    else:
        base_name = os.path.basename(cosmic_file).replace('.txt', '').replace('COSMIC_', '')
    
    print(f"Processing: {cosmic_file}")
    print(f"Output base name: {base_name}\n")
    
    # Load COSMIC data
    print("=" * 60)
    print("STEP 1: Loading COSMIC data")
    print("=" * 60)
    df = pd.read_csv(cosmic_file, sep='\t', index_col=0)
    print(f"Loaded {len(df)} mutations and {len(df.columns)} signatures\n")
    
    # Save original signatures to CSV
    signatures_csv = f'signatures_{base_name}.csv'
    df.to_csv(signatures_csv)
    print(f"Saved signatures to: {signatures_csv}\n")
    
    # Generate RePrint data
    print("=" * 60)
    print("STEP 2: Generating RePrint data")
    print("=" * 60)
    df_reprint = reprint(df)
    print(f"Generated RePrint for {len(df_reprint.columns)} signatures")
    
    # Save RePrint to CSV
    reprints_csv = f'reprints_{base_name}.csv'
    df_reprint.to_csv(reprints_csv)
    print(f"Saved RePrint data to: {reprints_csv}\n")
    
    # Generate standard signature plots (using plot.py)
    print("=" * 60)
    print("STEP 3: Generating standard signature plots")
    print("=" * 60)
    signature_output_dir = f'signature_plots_{base_name}'
    os.makedirs(signature_output_dir, exist_ok=True)
    
    # Generate single PDF with all signatures
    signature_pdf = f'all_signatures_{base_name}.pdf'
    try:
        save_all_signatures_to_single_pdf(
            df, 
            output_pdf=signature_pdf,
            prefix="", 
            yaxis_title="Probabilities"
        )
        print(f"\n✓ Standard signature plots PDF: {signature_pdf}\n")
    except Exception as e:
        print(f"\n⚠ Warning: Could not generate standard signature plots PDF: {e}")
        print("  (This requires kaleido: pip install -U kaleido)")
        print("  Continuing with RePrint plots...\n")
        signature_pdf = None
    
    # Generate RePrint plots (using plot_static.py)
    print("=" * 60)
    print("STEP 4: Generating RePrint plots (seamless matplotlib)")
    print("=" * 60)
    reprint_output_dir = f'reprint_plots_{base_name}'
    os.makedirs(reprint_output_dir, exist_ok=True)
    
    reprint_pdf = f'all_reprints_{base_name}.pdf'
    generate_all_reprints(
        reprints_csv, 
        output_pdf=reprint_pdf, 
        output_dir=reprint_output_dir
    )
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"✓ Signatures CSV: {signatures_csv}")
    print(f"✓ RePrints CSV: {reprints_csv}")
    if signature_pdf:
        print(f"✓ Standard signature plots PDF: {signature_pdf}")
    print(f"✓ RePrint plots PDF: {reprint_pdf}")
    print(f"✓ RePrint plots PNG: {reprint_output_dir}/")
    print("\nDone!")

if __name__ == "__main__":
    main()

