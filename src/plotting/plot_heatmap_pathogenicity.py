import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from src.project_config import get_paths, get_paths_protein
from src.plotting.plot_heatmap import plot_heatmap
from src.plotting.plot_pathogenicity import plot_pathogenicity


def plot_heatmap_pathogenicity(protein_id: str, span: str = "singlespan", curves: str = "both", 
                               smoothing_window: int = 5, save: bool = True, title: bool = True):
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
    gs = gridspec.GridSpec(nrows=7, ncols=1, height_ratios=[1, 0, 1, 0, 1, 0.1, 2])

    ax_am   = fig.add_subplot(gs[0])
    ax_esm  = fig.add_subplot(gs[2], sharex=ax_am)
    ax_diff = fig.add_subplot(gs[4])
    ax_path = fig.add_subplot(gs[6])

    # Plot heatmaps
    plot_heatmap(
        protein_id=protein_id,
        axes=[ax_am, ax_esm, ax_diff],
        add_colorbars=True
    )
 

    # Plot pathogenicity line plot
    plot_pathogenicity(
        protein_id=protein_id,
        curves=curves,
        rank=True,
        span=span,
        method="mean",
        highlight=["Helix", "Beta strand", "Transmembrane", "Juxtamembrane", "Disordered"],
        smoothing_window=smoothing_window,
        show_std=False,
        ax=ax_path,
        fig=fig
    )

    # Y-axis labels
    ax_am.set_ylabel("AlphaMissense Rank", fontsize=20, fontweight='bold')
    ax_esm.set_ylabel("ESM-1b Rank", fontsize=20, fontweight='bold')
    ax_diff.set_ylabel("Difference Rank", fontsize=20, fontweight='bold')

    # X-axis label only on bottom plot
    ax_path.set_xlabel("Residue Position")

    # Title
    if title == True:
        ax_am.set_title(f'Rank Pathogenicity of AlphaMissense and ESM1b for {protein_id}', fontsize=48, fontweight='bold')

    
    # Remove x ticks from heatmaps of AlphaMissense and ESM1b
    ax_am.set_xticks([])
    ax_am.tick_params(axis="x", which="both", length=0)
    

    # Adjust layout
    plt.subplots_adjust(
        hspace=0.03,
        bottom=0.05,
        top=0.95,
        left=0.05,
        right=0.98
    )

    # Save plot
    if save:
        images_path.mkdir(parents=True, exist_ok=True)
        output_path = images_path / f"combo_heatmap_pathogenicity_{protein_id}.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_path}")

    plt.show()