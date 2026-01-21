import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import os


def format_context_label(context):
    """
    Format context label - simple format without bold.
    Converts 'A[C>A]A' or 'ACA' to 'ACA' format.
    """
    if '[' in context and ']' in context:
        # Format: 'A[C>A]A' -> extract left, center, right
        parts = context.split('[')
        left = parts[0]
        mutation_and_right = parts[1].split(']')
        center = mutation_and_right[0].split('>')[0]  # Get ref base (C or T)
        right = mutation_and_right[1]
    else:
        # Format: 'ACA' -> extract left, center, right
        left = context[0]
        center = context[1]
        right = context[2]
    
    return f'{left}{center}{right}'


def create_main_dashboard(df, signature, title, yaxis_title, 
                          show_x_labels=True, show_y_labels=True,
                          figsize=None):
    """
    Original vertical bar chart implementation using matplotlib.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with signatures as columns
    signature : str
        Name of the signature column
    title : str
        Title for the plot
    yaxis_title : str
        Title for y-axis
    show_x_labels : bool, optional
        Whether to show X-axis labels (default: True)
    show_y_labels : bool, optional
        Whether to show Y-axis labels (default: True)
    figsize : tuple of float, optional
        Figure size (width, height) in inches. Default: (16, 5) for better X-axis visibility
    """
    frequencies = df[signature] * 1

    mutations = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    bases = ['A', 'C', 'G', 'T']
    contexts = [f'{x}[{m}]{y}' for m in mutations for x in bases for y in bases]

    colors = {
        'C>A': 'blue',
        'C>G': 'black',
        'C>T': 'red',
        'T>A': 'gray',
        'T>C': 'green',
        'T>G': 'pink'
    }

    if figsize is None:
        figsize = (16, 5)  # Wider aspect ratio for better X-axis visibility
    fig, ax = plt.subplots(figsize=figsize)
    
    # Calculate positions for bars
    n_contexts = len(contexts)
    x_pos = np.arange(n_contexts)
    width = 0.8  # Width of each bar
    
    # Plot bars for each mutation type
    for mutation in mutations:
        mutation_contexts = [c for c in contexts if f'[{mutation}]' in c]
        mutation_frequencies = [frequencies[mc] if mc in frequencies.index else 0 for mc in mutation_contexts]
        
        # Find positions for this mutation's contexts
        positions = [x_pos[contexts.index(mc)] for mc in mutation_contexts]
        
        ax.bar(positions, mutation_frequencies, 
               width=width, label=mutation, color=colors[mutation], alpha=0.8)

    y_max = frequencies.max()

    ax.set_title(title, fontsize=14, pad=20)
    
    if show_x_labels:
        # Remove X-axis label
        ax.set_xlabel('')
        # Set ticks but make them invisible
        ax.set_xticks(x_pos)
        # Format labels (simple format) - shift right by half bar width to center with bars
        formatted_labels = [format_context_label(ctx) for ctx in contexts]
        ax.set_xticklabels(formatted_labels, rotation=90, ha='center', fontsize=8, family='monospace')
        # Remove tick marks (keep labels)
        ax.tick_params(axis='x', which='both', length=0)
    else:
        ax.set_xlabel('')
        ax.set_xticks([])
        ax.set_xticklabels([])
    
    ax.set_xlim(-0.7, n_contexts - 0.3) 

    if show_y_labels:
        ax.set_ylabel(yaxis_title, fontsize=12, labelpad=2)
        ax.tick_params(axis='y', which='major', pad=2, labelsize=10)
    else:
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.set_yticklabels([])

    ax.set_ylim(0, y_max * 1.05)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Remove frame (all spines)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Legend for mutation types, placed outside but close to the plotting area
    legend_handles = []
    legend_labels = []
    for mut in mutations:
        legend_handles.append(mpatches.Patch(color=colors[mut], label=mut))
        legend_labels.append(mut)
    fig.legend(
        legend_handles,
        legend_labels,
        title=None,
        loc='upper right',
        bbox_to_anchor=(0.98, 0.98),
        frameon=False,
        fontsize=10,
        handlelength=1.0,
        handleheight=1.0,
    )

    # Adjust margins: no extra space on the left, reserve room on the right for legend
    plt.subplots_adjust(left=0.02, right=0.92, top=0.90, bottom=0.20)
    return fig


def create_main_dashboard_horizontal(df, signature, title, yaxis_title, figsize=None):
    """
    Horizontal bar chart with 3 rows, each showing a pair of complementary mutations.
    Each row has 2 groups of 16 bars (C-centered and T-centered contexts).
    Rows in order:
    - C>T (red, left 16) + T>G (light pink, right 16)
    - C>G (black, left 16) + T>C (light green, right 16)
    - C>A (light blue, left 16) + T>A (gray, right 16)
    
    X-axis shows 32 trinucleotide contexts: first 16 (ACA-TCT) then last 16 (ATA-TTT)
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with signatures as columns
    signature : str
        Name of the signature column
    title : str
        Title for the plot
    yaxis_title : str
        Title for y-axis
    figsize : tuple of float, optional
        Figure size (width, height) in inches. Default: (16, 6) for better X-axis visibility
    """
    frequencies = df[signature] * 1

    # Define the order of 32 trinucleotide contexts for X-axis
    # First 16: C-centered contexts
    c_contexts = [
        'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT',
        'GCA', 'GCC', 'GCG', 'GCT', 'TCA', 'TCC', 'TCG', 'TCT'
    ]
    # Last 16: T-centered contexts
    t_contexts = [
        'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG', 'CTT',
        'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC', 'TTG', 'TTT'
    ]
    # Combined order for X-axis
    context_order = c_contexts + t_contexts

    # Define mutation pairs: (C-mutation, T-mutation) for each row
    mutation_pairs = [
        ('C>T', 'T>G'),  # Top row
        ('C>G', 'T>C'),  # Middle row
        ('C>A', 'T>A')   # Bottom row
    ]
    
    def get_value_for_context(context_short, mutation_type):
        """
        Get value for a short context (e.g., 'ACA') and mutation type (e.g., 'C>A').
        Converts 'ACA' + 'C>A' -> 'A[C>A]A'
        """
        full_context = f'{context_short[0]}[{mutation_type}]{context_short[2]}'
        return frequencies[full_context] if full_context in frequencies.index else 0

    # Color scheme for each mutation pair
    color_schemes = {
        ('C>T', 'T>G'): {
            'left': 'red',        # C>T contexts (left 16)
            'right': 'lightpink'  # T>G contexts (right 16)
        },
        ('C>G', 'T>C'): {
            'left': 'black',      # C>G contexts (left 16)
            'right': 'lightgreen' # T>C contexts (right 16)
        },
        ('C>A', 'T>A'): {
            'left': 'lightblue',  # C>A contexts (left 16)
            'right': 'gray'       # T>A contexts (right 16)
        }
    }

    # Get overall min/max for consistent y-axis scaling
    all_values = []
    for c_mut, t_mut in mutation_pairs:
        for ctx in c_contexts:
            val = get_value_for_context(ctx, c_mut)
            all_values.append(val)
        for ctx in t_contexts:
            val = get_value_for_context(ctx, t_mut)
            all_values.append(val)
    
    y_min = min(all_values)
    y_max = max(all_values)
    y_range = max(abs(y_min), abs(y_max)) * 1.1  # Add 10% padding

    # Create subplots with 3 rows, one for each mutation pair
    if figsize is None:
        figsize = (16, 6)  # Wider aspect ratio for better X-axis visibility
    fig, axes = plt.subplots(
        len(mutation_pairs), 1,
        figsize=figsize,
        sharex=True,
        gridspec_kw={'hspace': 0.05}
    )
    
    # Ensure axes is a list even if there's only one subplot
    if len(mutation_pairs) == 1:
        axes = [axes]
    
    x_pos = np.arange(len(context_order))
    width = 0.4  # Width of bars

    for row_idx, (c_mut, t_mut) in enumerate(mutation_pairs):
        ax = axes[row_idx]
        
        # Get values for C-mutation contexts (first 16)
        c_values = [get_value_for_context(ctx, c_mut) for ctx in c_contexts]
        # Get values for T-mutation contexts (last 16)
        t_values = [get_value_for_context(ctx, t_mut) for ctx in t_contexts]
        
        # Rotate upper two rows (multiply by -1)
        if row_idx <= 2:
            c_values = [-v for v in c_values]
            t_values = [-v for v in t_values]
        
        colors = color_schemes[(c_mut, t_mut)]
        
        # Positions for C and T contexts
        c_positions = x_pos[:len(c_contexts)]
        t_positions = x_pos[len(c_contexts):]
        
        # Add bars for C-mutation contexts (left 16 bars)
        ax.bar(c_positions, c_values, width=width, color=colors['left'], 
               label=c_mut, alpha=0.8)
        
        # Add bars for T-mutation contexts (right 16 bars)
        ax.bar(t_positions, t_values, width=width, color=colors['right'], 
               label=t_mut, alpha=0.8)
        
        # Add horizontal line at y=0 for baseline
        ax.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.3)
        
        # Update y-axis for this subplot - hide axes, symmetric range
        ax.set_ylim(-y_range, y_range)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.grid(False)

    # Update x-axis (only for bottom subplot)
    axes[-1].set_xticks(x_pos)
    # Remove X-axis label
    axes[-1].set_xlabel('')
    # Format labels (simple format)
    formatted_labels = [format_context_label(ctx) for ctx in context_order]
    axes[-1].set_xticklabels(formatted_labels, rotation=90, ha='center', fontsize=8, family='monospace')
    # Remove tick marks (keep labels)
    axes[-1].tick_params(axis='x', which='both', length=0)
    axes[-1].spines['top'].set_visible(False)
    axes[-1].spines['right'].set_visible(False)
    axes[-1].spines['left'].set_visible(False)
    axes[-1].spines['bottom'].set_visible(False)
    
    # Hide x-axis labels for upper rows
    for row_idx in range(len(mutation_pairs) - 1):
        axes[row_idx].set_xticks([])
        axes[row_idx].set_xticklabels([])
    
    plt.suptitle(title, fontsize=14, y=0.98)
    
    # Adjust margins: reduce left, increase right
    plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.15, hspace=0.05)
    
    return fig


def save_all_signatures_to_pdf(df, output_dir=".", prefix="reprint_", yaxis_title="Probabilities", figsize=None):
    """
    Creates and saves PDF plots for all signatures in the DataFrame.
    Each signature is saved as a separate PDF file.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with signatures as columns
    output_dir : str, optional
        Directory to save PDF files (default: ".")
    prefix : str, optional
        Prefix for output PDF files (default: "reprint_")
    yaxis_title : str, optional
        Title for y-axis (default: "Probabilities")
    figsize : tuple of float, optional
        Figure size (width, height) in inches. Default: (16, 5) for better X-axis visibility
    """
    os.makedirs(output_dir, exist_ok=True)
    for signature in df.columns:
        # Format title as "SBS1 signature" style
        title = f"{signature} signature"
        
        fig = create_main_dashboard(
            df,
            signature=signature,
            title=title,
            yaxis_title=yaxis_title,
            figsize=figsize
        )
        pdf_path = os.path.join(output_dir, f"{prefix}{signature}.pdf")
        fig.savefig(pdf_path, format="pdf", bbox_inches='tight', dpi=300)
        plt.close(fig)
        print(f"Saved: {pdf_path}")


def save_all_signatures_to_single_pdf(df, output_pdf, prefix="", yaxis_title="Probabilities", figsize=None):
    """
    Creates and saves all signature plots in a single PDF file.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with signatures as columns
    output_pdf : str
        Path to output PDF file
    prefix : str
        Prefix for signature names in titles
    yaxis_title : str
        Title for y-axis
    figsize : tuple of float, optional
        Figure size (width, height) in inches. Default: (16, 5) for better X-axis visibility
    """
    print(f"\nCombining {len(df.columns)} plots into PDF...")
    with PdfPages(output_pdf) as pdf:
        for i, signature in enumerate(df.columns, 1):
            print(f"  Processing {i}/{len(df.columns)}: {signature}")
            try:
                # Format title as "SBS1 signature" style
                title = f"{signature} signature"
                
                fig = create_main_dashboard(
                    df,
                    signature=signature,
                    title=title,
                    yaxis_title=yaxis_title,
                    figsize=figsize
                )
                
                # Save to PDF with high DPI
                pdf.savefig(fig, bbox_inches='tight', dpi=300)
                plt.close(fig)
            except Exception as e:
                print(f"    Warning: Could not add {signature} to PDF: {e}")
                continue
        
    print(f"✓ Combined PDF saved: {output_pdf}")
