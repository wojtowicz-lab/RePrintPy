import argparse
import os
import tempfile
import pandas as pd
import matplotlib.pyplot as plt
from reprint_tool.plot_static import generate_all_reprints
from reprint_tool.core import reprint
from reprint_tool.analyze import (
    create_heatmap_with_custom_sim,
    calculate_rmse,
    calculate_cosine,
)


def main():
    parser = argparse.ArgumentParser(description="Generate PDF plots for all signatures in a DataFrame.")
    parser.add_argument('--input', required=True, help='Path to input CSV/TSV file (with signatures as columns)')
    parser.add_argument('--output_dir', default='output', help='Directory to save PDF files (default: output)')
    parser.add_argument('--prefix', default='reprint_', help='Prefix for output PDF files (default: reprint_)')
    parser.add_argument('--sep', default='\t', help='Column separator in input file (default: tab)')
    parser.add_argument('--save_reprint', default=None, help='If set, save computed reprint DataFrame to this file (CSV/TSV)')
    parser.add_argument('--export_png', action='store_true', help='Also export per-signature plots to PNG files')
    parser.add_argument('--analyze_png', default=None, help='If set, generate similarity heatmap and save to this PNG path')
    parser.add_argument('--analyze_metric', default='rmse', choices=['rmse', 'cosine'], help='Similarity metric for heatmap')
    parser.add_argument('--analyze_colorscale', default='Blues', help='Colorscale for analyze heatmap')
    parser.add_argument('--analyze_hide_heatmap', action='store_true', help='Hide heatmap, show only dendrograms')
    parser.add_argument('--analyze_method', default='complete', help='Linkage method for clustering (e.g., complete, average)')
    args = parser.parse_args()
    
    # Convert escape sequences like '\t' to actual tab character
    sep = args.sep.encode().decode('unicode_escape')

    df = pd.read_csv(args.input, sep=sep, index_col=0)
    df_reprint = reprint(df)

    # Save reprint DataFrame to CSV (needed for plot_static functions)
    if args.save_reprint:
        os.makedirs(os.path.dirname(args.save_reprint) or '.', exist_ok=True)
        df_reprint.to_csv(args.save_reprint, sep=sep)
        print(f"Saved reprint DataFrame to {args.save_reprint}")
        reprint_csv_path = args.save_reprint
    else:
        # Create temporary CSV file for plot_static functions
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp_file:
            df_reprint.to_csv(tmp_file.name, sep=sep)
            reprint_csv_path = tmp_file.name
    
    # Generate reprint plots using plot_static
    if args.export_png:
        # Generate PNG files
        generate_all_reprints(
            csv_path=reprint_csv_path,
            output_pdf=None,
            output_dir=args.output_dir,
            sep=sep
        )
    else:
        # Generate PDF file with all plots
        pdf_path = os.path.join(args.output_dir, f"{args.prefix}all_reprints.pdf")
        os.makedirs(args.output_dir, exist_ok=True)
        generate_all_reprints(
            csv_path=reprint_csv_path,
            output_pdf=pdf_path,
            output_dir=None,
            sep=sep
        )
    
    # Clean up temporary file if we created one
    if not args.save_reprint:
        try:
            os.unlink(reprint_csv_path)
        except:
            pass

    if args.analyze_png:
        metric_func = calculate_rmse if args.analyze_metric == 'rmse' else calculate_cosine
        fig = create_heatmap_with_custom_sim(
            df_reprint,
            calc_func=metric_func,
            colorscale=args.analyze_colorscale,
            hide_heatmap=args.analyze_hide_heatmap,
            method=args.analyze_method,
        )
        out_path = args.analyze_png
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        fig.savefig(out_path, format='png', dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved analyze heatmap: {out_path}")


if __name__ == "__main__":
    main()