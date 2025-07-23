from src.project_config import get_paths, COLORS_MODELS, COLORS_SECONDARY
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import math
from matplotlib import rcParams




# A helper function for plot_rank_dotplot()
def round_up_nice(x):
    """Round up to a clean, compact tick value near x."""
    if x == 0:
        return 0
    exponent = math.floor(math.log10(x))
    fraction = x / 10**exponent

    # More fine-grained thresholds
    if fraction <= 1.5:
        nice = 1.5
    elif fraction <= 2:
        nice = 2
    elif fraction <= 3:
        nice = 3
    elif fraction <= 5:
        nice = 5
    else:
        nice = 10

    return nice * 10**exponent


# Find optimal number of bins for histogram
def optimal_bin_count(n):
    """Adaptive bin count: fewer bins for massive data"""
    if n > 2_000_000:
        return 60 # hard cap for huge datasets
    elif n > 500_000:
        return 30
    else:
        return max(10, math.ceil(2 * n ** (1/3)))  # Rice Rule






# Main Function

def plot_rank_dotplot(model, alpha=0.1, sample="human", sample_size=1000, save=False):
    """
    Plots a rank dotplot and histogram of raw and rank pathogenicity scores for protein variants.
    Parameters:
    model (str): The model used for scoring ('AlphaMissense' or 'ESM').
    alpha (float): The transparency level of the scatter plot points (default is 0.1).
    sample (str): The type of sample to use ('human', 'nout', or 'multispan').
    sample_size (int): The number of protein variants to sample (default is 1000).
    save (bool): If True, saves the plot as a PNG file (default is False).
    Raises:
    ValueError: If sample_size exceeds the maximum allowed for the specified sample type.
    The function prepares datasets, loads protein IDs, collects raw and ranked scores,
    and generates a scatter plot with a histogram to visualize the relationship between
    raw scores and rank scores. It also saves the plot to a specified directory if requested.
    """
    # Validate model input
    if sample not in ["human", "nout", "multispan"]:
        raise ValueError("Sample must be 'human', 'nout', or 'multispan'.")

    # Prepare datasets and folder pathways
    protein_ids_csv = get_paths()["processed"] / "Protein_IDs_Per_Experiment" / "intersection_protein_ids_to_be_ranked.csv"
    if model == "AlphaMissense":
        raw_folder = get_paths()['data'] / "raw" / "AlphaMissense_csv"
        ranked_folder = get_paths()['am_path'] 
    elif model == "ESM":
        raw_folder = get_paths()['processed'] / "ESM1b_no_isoforms"
        ranked_folder = get_paths()['esm_path'] 

    # Standard amino acids order for matrix data structure
    aa_list = [
        "A", "V", "L", "I", "M", "F", "W",  # Hydrophobic
        "S", "T", "N", "Q", "Y", "C",      # Polar uncharged
        "K", "R", "H",                     # Positively charged
        "D", "E",                          # Negatively charged
        "G", "P"                           # Special
    ]


    # Load protein IDs
    if sample == "human":
        protein_ids = pd.read_csv(protein_ids_csv)['Protein_ID'].sample(n=sample_size, random_state=42).tolist()

    elif sample == "nout":
        if sample_size > 1534: # N-Out proteome has 1534 proteins
            raise ValueError("Sample size for N-Out proteome cannot exceed 10,000.")
        protein_ids = pd.read_csv(get_paths()["n_out_id"])['entry'].sample(n=sample_size, random_state=42).tolist()

    elif sample == "multispan":
        if sample_size > 2822:
            raise ValueError("Sample size for Multispan proteome cannot exceed 2,822.")
        protein_ids = pd.read_csv(get_paths()["processed"] / "Multispan_Proteome_cleaned.csv")['Entry'].sample(n=sample_size, random_state=42).tolist()

   
    # Prepare lists to hold raw and rank scores
    raw_scores = []
    rank_scores = []

    # Modify CSV files for proteins for correct format
    if model == "ESM":
        for protein_id in tqdm(protein_ids, desc="Reading ESM1b scores"):
            # Load corresponding raw and ranked data for the protein
            # Reorder rows to standard AA list
            raw_file = os.path.join(raw_folder, f"{protein_id}_LLR.csv")
            ranked_file = os.path.join(ranked_folder, f"{protein_id}_LLR_rank.csv")

            if os.path.exists(raw_file) and os.path.exists(ranked_file):

                # Reorder rows to standard AA list and replace 0 for NaN for raw data like already in ranked data
                raw_df = pd.read_csv(raw_file, index_col = 0).reindex(aa_list).replace(0, np.nan)
                ranked_df = pd.read_csv(ranked_file, index_col = 0)

                # Collect scores, ignoring NaN values
                raw_scores.extend(raw_df.iloc[:, 1:].values.flatten()[~np.isnan(raw_df.iloc[:, 1:].values.flatten())])
                rank_scores.extend(ranked_df.iloc[:, 1:].values.flatten()[~np.isnan(ranked_df.iloc[:, 1:].values.flatten())])

    elif model == "AlphaMissense":
        for protein_id in tqdm(protein_ids, desc="Reading AlphaMissense scores"):
            # Load corresponding raw and ranked data for the protein
            # Reorder rows to standard AA list
            raw_file = os.path.join(raw_folder, f"{protein_id}.csv")
            ranked_file = os.path.join(ranked_folder, f"{protein_id}_rank.csv")

            if os.path.exists(raw_file) and os.path.exists(ranked_file):
                # Load the raw AlphaMissense CSV
                raw_df = pd.read_csv(raw_file)
                ranked_df = pd.read_csv(ranked_file, index_col=0)

                # Create column labels like "M 1", "A 2", etc.
                raw_df["column_label"] = raw_df["residue"] + " " + raw_df["residue_position"].astype(str)

                # Pivot the data: rows=variant AA (e.g., "A", "V"), cols=positions like "M 1", values=am_pathogenicity
                raw_df = raw_df.pivot(index="variation", columns="column_label", values="am_pathogenicity")
                raw_df = raw_df[sorted(raw_df.columns, key=lambda x: int(x.split()[1]))]
                
                # Clean up: set index to variation, drop unwanted row/column if it appears
                raw_df.index.name = None  # remove index name "variation" if shown
                raw_df.columns.name = None  # remove column level name "column_label" if shown
                ranked_df.index.name = None  # remove index name "variation" if shown
            
                # Reindex rows to standard amino acid list (some may be missing → become NaN)
                raw_df = raw_df.reindex(aa_list)
                
                # Collect flattened non-NaN scores
                raw_scores.extend(raw_df.iloc[:, 1:].values.flatten()[~np.isnan(raw_df.iloc[:, 1:].values.flatten())])
                rank_scores.extend(ranked_df.iloc[:, 1:].values.flatten()[~np.isnan(ranked_df.iloc[:, 1:].values.flatten())])

                

    ### --- Plotting --- ###

    # Create a DataFrame for plotting
    plot_df = pd.DataFrame({'Raw Score': raw_scores, 'Rank': rank_scores})
    # Set general style
    sns.set_style("whitegrid")
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.titlesize': 12,
        'axes.titleweight': 'bold',
        'axes.labelsize': 11,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10
    })


    # Define colors for AM and ESM models
    scatter_color = COLORS_MODELS["ESM"] if model == "ESM" else COLORS_MODELS["AM"]  # Colorblind-safe green = # ESM color '#29ab88'
    # Slightly darker version of the main color
    

    # Create figure with 2 subplots: scatter + histogram (right side)
    fig = plt.figure(figsize=(6.8, 6.8), dpi=300)
    gs = fig.add_gridspec(nrows=2, ncols=2, width_ratios=[6, 1], height_ratios=[1, 6], wspace=0.05, hspace=0.05)

    # Scatter plot on the left
    ax_scatter = fig.add_subplot(gs[1,0])
    scatter = ax_scatter.scatter(
        plot_df['Rank'], plot_df['Raw Score'],
        alpha=alpha, s=18, color=scatter_color)
    

    xlabel_text = (
        r"$\bf{Rank\ Score}$" + "\n" + "Pathogenicity Increases →"
        if model == "ESM"
        else r"$\bf{Rank\ Score}$" + "\n" + "Pathogenicity Increases →"
    )
    ax_scatter.set_xlabel(xlabel_text, fontsize=8, multialignment='center')
    # Set the ylabel with different styles for each line
    ylabel_text = (
        r"$\bf{Raw\ LLR\ Score}$" + "\n" + "Pathogenicity Increases →"
        if model == "ESM"
        else r"$\bf{Raw\ Score}$" + "\n" + "Pathogenicity Increases →"
    )

    ax_scatter.set_ylabel(ylabel_text, fontsize=8, multialignment='center')
    ax_text = fig.add_subplot(gs[0, 1])
    ax_text.axis('off')  # Hide all axes lines, ticks, etc.

    # Compose text
    title_text = (
        'ESM1b: Rank Dotplot\nand Histogram of\nPathogenicity Scores'
        if model == "ESM"
        else 'AlphaMissense:\nRank Dotplot\nand Histogram of\nPathogenicity Scores'
    )
    subtitle_text = f'Protein Sample Size: {sample_size:,}\nVariants: {len(raw_scores):,}'

    # Add text to the axes
    ax_text.text(
        0.0, 1, title_text,
        fontsize=6, fontweight='bold', ha='left', va='top'
    )
    ax_text.text(
        0.0, 0.45, subtitle_text,
        fontsize=6, ha='left', va='top'
    )
    ax_text.text(0.0, 0.2, "Rank\nTransformation", fontsize=6., fontweight='bold', ha='left', va='top', )
 
    
    ax_scatter.grid(True, linestyle='--', alpha=0.2)
    ax_scatter.tick_params(axis='both', labelsize=8)
    
  
    
    if model == "ESM":
        ax_scatter.axhline(y=-7.5, color=COLORS_SECONDARY["PURPLE"], linestyle='-', linewidth=1.5, xmin=-0.1, xmax=0.48)
        ax_scatter.axvline(x=0.48271, color=COLORS_SECONDARY["PURPLE"], linestyle='-', linewidth=1.5, ymin=0.0, ymax=0.480)
        ax_scatter.text(-0.0425, -8.2, 'Pathogenicity Threshold ↑', color=COLORS_SECONDARY["PURPLE"], fontsize=6, va='center')
        ax_scatter.text(-0.130, -7.5, '-7.5', color=COLORS_SECONDARY["PURPLE"], fontsize=7, va='center')
        ax_scatter.text(0.45, 16.5, '0.483', color=COLORS_SECONDARY["PURPLE"], fontsize=7, va='center')
        ax_scatter.text(0.49, 4, 'Rank Pathogenicity Threshold ↑', color=COLORS_SECONDARY["PURPLE"], fontsize=6, va='center', rotation = 270)

    elif model == "AlphaMissense":
        ax_scatter.axhline(y=0.564, color=COLORS_SECONDARY["PURPLE"], linestyle='-', linewidth=1.5, xmin=-0.1, xmax=0.559)
        ax_scatter.axvline(x=0.56548, color=COLORS_SECONDARY["PURPLE"], linestyle='-', linewidth=1.5, ymin=0.0, ymax=0.554)
        ax_scatter.text(-0.0425, 0.585 , 'Pathogenicity Threshold ↑', color=COLORS_SECONDARY["PURPLE"], fontsize=6, va='center')
        ax_scatter.text(-0.135 , 0.564, '0.56', color=COLORS_SECONDARY["PURPLE"], fontsize=7, va='center')
        ax_scatter.text(0.54 , -0.06, '0.56', color=COLORS_SECONDARY["PURPLE"], fontsize=7, va='center')
        ax_scatter.text(0.575, 0.25 , 'Rank Pathogenicity Threshold ↑', color=COLORS_SECONDARY["PURPLE"], fontsize=6, va='center', rotation = 270)
        
            
            
        

    # Flip y-axis: more negative = more pathogenic → should be higher
    if model == "ESM": ax_scatter.invert_yaxis()
        
    n_points = len(raw_scores)/5
    num_bins = optimal_bin_count(n_points)

    # Histogram on the right — aligned with y-axis of scatter
    ax_hist = fig.add_subplot(gs[1,1], sharey=ax_scatter)
    hist_values, bins, _ = ax_hist.hist(
        plot_df['Raw Score'], bins = 200, orientation='horizontal',
        color=scatter_color, alpha=0.6, edgecolor='black', linewidth=0.3)

    # Only show ticks at 0 and max
    # After histogram is created:
    max_freq = max(hist_values)
    #nice_max = round_up_nice(max_freq)

    ax_hist.set_xticks([0, round(max_freq)])
    ax_hist.set_xticklabels(['0', rf'$ {int(max_freq / 10 ** int(math.log10(max_freq)))} \times 10^{{{int(math.log10(max_freq))}}} $'], fontsize=8)

    # Clean up histogram axis
    ax_hist.tick_params(axis='y', labelleft=False)  # Hide y-axis ticks on histogram
    ax_hist.set_xlabel('Frequency', fontsize=8)
    ax_hist.grid(False)



     # Histogram on the top — aligned with x-axis of scatter
    ax_hist = fig.add_subplot(gs[0,0], sharex=ax_scatter)
    hist_values, bins, _ = ax_hist.hist(
        plot_df['Rank'], bins = 200, orientation='vertical',
        color=scatter_color, alpha=0.6, edgecolor='black', linewidth=0.3)

    # Only show ticks at 0 and max
    #max_freq = max(hist_values)
    #nice_max = round_up_nice(max_freq)
    ax_hist.set_yticks([0, round(max_freq)])
    ax_hist.set_yticklabels(['0', rf'$ {int(max_freq / 10 ** int(math.log10(max_freq)))} \times 10^{{{int(math.log10(max_freq))}}} $'], fontsize=8)

    # Clean up histogram axis
    ax_hist.tick_params(axis='x', labelbottom=False)  # Hide x-axis ticks on histogram
    ax_hist.set_ylabel('Frequency', fontsize=8, labelpad=-8)
    ax_hist.grid(False)






     # Save plot
    if save:
        images_path = get_paths()['images_path'] / '5.1.Ranking Dotplot'
        images_path.mkdir(parents=True, exist_ok=True)
        output_path = images_path / f"ESM_Rank_Dotplot_Sample_Size_{sample_size}.png" if model == "ESM" else images_path / f"AlphaMissense_Rank_Dotplot_Sample_Size_{sample_size}.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

        relative_output = os.path.relpath(output_path, start=get_paths()["project_root"])
        print(f"✅ Rank Dotplot saved to '{relative_output}'.")

    plt.tight_layout()
    plt.show()