import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from src.project_config import get_paths, get_paths_protein, COLORS_MODELS
from src.plotting.plot_heatmap import plot_heatmap
from src.plotting.plot_pathogenicity import plot_pathogenicity
from matplotlib.patches import Rectangle


def plot_heatmap_pathogenicity(protein_id: str, span: str = "singlespan", curves: str = "both", 
                               smoothing_window: int = 5, save: bool = False, title: bool = True,
                               rASA=False, secondary_str_version: str = "dssp"):
    """
    Plots a composite figure with heatmaps (AlphaMissense, ESM1b, Difference) and pathogenicity line plot.
    
    Parameters:
    - protein_id (str): UniProt protein ID.
    - smoothing_window (int): Smoothing window for pathogenicity curves.
    - save (bool): Whether to save the figure as PNG.
    """

    paths = get_paths_protein(protein_id)
    pathways = get_paths()
    images_path = pathways["images_path"] / "5.6.Combo_Heatmap_Average_Pathogenicity"

    # Read difference heatmap data (just for adjusting x-lenght of the figure)
    data_diff = pd.read_csv(paths["difference_path"], index_col=0)

    # Create figure and gridspec layout
    fig = plt.figure(figsize=(len(data_diff.columns) / 15, 18))

    if rASA:
        gs = gridspec.GridSpec(nrows=9, ncols=1, height_ratios=[1, 0, 1, 0.075, 1, 0.1, 2, 0.075, 0.5])
    else:
        gs = gridspec.GridSpec(nrows=7, ncols=1, height_ratios=[1, 0, 1, 0, 1, 0.1, 2])

    ax_am   = fig.add_subplot(gs[0])
    ax_esm  = fig.add_subplot(gs[2], sharex=ax_am)
    ax_clusters = fig.add_subplot(gs[3])
    ax_diff = fig.add_subplot(gs[4])
    ax_path = fig.add_subplot(gs[6])
    ax_rasa = fig.add_subplot(gs[8], sharex=ax_path) if rASA else None

    # Plot heatmaps
    plot_heatmap(
        protein_id=protein_id,
        axes=[ax_am, ax_esm, ax_diff],
        add_colorbars=True
    )
    
    # Load cluster annotation data
    try:
        cluster_csv_path = get_paths()["clusters_path"] / f"{protein_id}_clusters.csv" 
        cluster_df = pd.read_csv(cluster_csv_path)

        # Plot each cluster region as a rectangle
        for _, row in cluster_df.iterrows():
            try:
                start, end = row['range'].split('-')
                start, end = int(start), int(end)
            except:
                continue

            color = COLORS_MODELS['AM'] if row['more pathogenic'] == 'AlphaMissense' else COLORS_MODELS['ESM']
            # Temporary debug: draw background line
            ax_clusters.hlines(y=0.5, xmin=start, xmax=end, colors=color, linestyles="-", linewidth=10, rasterized=True)
    except:
        print(f"No cluster data found for {protein_id}. Skipping cluster annotation.")





    ax_clusters.set_ylim(0, 1)
    ax_clusters.set_xlim(0, len(data_diff.columns))
    ax_clusters.set_yticks([])
    ax_clusters.set_xticks([])
    ax_clusters.set_facecolor("white")
    ax_clusters.tick_params(axis="x", which="both", length=0)







    # Plot pathogenicity line plot
    plot_pathogenicity(
        protein_id=protein_id,
        curves=curves,
        rank=True,
        span=span,
        method="mean",
        highlight=["Helix", "Beta strand", "Transmembrane", "Juxtamembrane", "Coil; Bend; Disordered"],
        smoothing_window=smoothing_window,
        show_std=False,
        ax=ax_path,
        fig=fig,
        secondary_str_version=secondary_str_version
    )


    if rASA:
        # Get paths and select the path and rankscore
        paths = get_paths_protein(protein_id=protein_id)
        data = paths["dssp_protein_path"]
        df = pd.read_csv(data, index_col=0 if secondary_str_version == "uniprot" else None)
        # Define x-axis column
        x_col = 'residue_position'
        if x_col not in df.columns:
            raise KeyError(f"Column '{x_col}' missing in CSV.")


        ax_rasa.tick_params(axis='x', labelsize=8)
        # Plot rASA as a line (or points) on its own axis
        ax_rasa.plot(df[x_col], df['rASA'], color='black', linewidth=2, rasterized=True)
        ax_rasa.set_ylabel('rASA', fontsize=20, labelpad=2, fontweight='bold')
        ax_rasa.set_ylim(0, 1.01)
        ax_rasa.tick_params(axis='x', labelsize=30, rotation = 45)
        ax_rasa.tick_params(axis='y', labelsize=20)
        ax_rasa.set_yticks([0, 0.5, 1.0])
        ax_rasa.tick_params(labelbottom=True, labelsize=18)
        ax_rasa.grid(False)
        ax_rasa.set_xlabel("Residue Position", fontsize=36, fontweight='bold')
        
        # Hide redundant x-axis on top plot
        
        ax_path.tick_params(labelbottom=False)
        ax_path.set_xlabel("")




    # Y-axis labels
    ax_am.set_ylabel("AlphaMissense\nRank Scores", fontsize=20, fontweight='bold')
    ax_esm.set_ylabel("ESM-1b\nRank Scores", fontsize=20, fontweight='bold')
    ax_diff.set_ylabel("Difference\nRank Scores", fontsize=20, fontweight='bold')


    # X-axis label and x-ticks only on bottom plot
    if not rASA: 
        ax_path.set_xlabel("Residue Position",fontsize=36, fontweight='bold')


    # Title
    if title == True:
        ax_am.set_title(f'Rank Pathogenicity of AlphaMissense and ESM1b for {protein_id}', fontsize=48, fontweight='bold')

    
    # Remove x ticks from heatmaps of AlphaMissense and ESM1b
    ax_am.set_xticks([])
    ax_am.tick_params(axis="x", which="both", length=0)
    

    

    # Adjust layout
    plt.subplots_adjust(
        hspace=0.03,
        bottom=0.1,
        top=0.95,
        left=0.05,
        right=0.98
    )

    # Save plot
    if save:
        images_path.mkdir(parents=True, exist_ok=True)
        output_path = images_path / f"combo_heatmap_pathogenicity_{protein_id}_.png" 
        plt.savefig(output_path, dpi=100, bbox_inches='tight')
        print(f"Plot saved to {output_path}")

    plt.show()