import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from src.project_config import get_paths_protein, get_paths
import os



def plot_heatmap(
    protein_id="P05067",
    axes=None,
    figsize=None,
    add_colorbars=True,
    save=False
):
    """
    Unified heatmap plotting function that works for both standalone and embedded cases.
    
    Parameters:
    - protein_id: Identifier for the protein
    - axes: Optional list of three Matplotlib Axes objects. If None, a new figure is created.
    - figsize: Tuple for figure size if axes not provided.
    - add_colorbars: Whether to add colorbars (only shown if axes not provided or if desired for embedded use).
    
    Returns:
    - fig: Figure object (or None if no new figure created)
    - axes: List of three Axes objects
    """
    
    # Load data and pathways for the specified protein
    paths = get_paths_protein(protein_id)
    data_diff = pd.read_csv(paths['difference_path'], index_col=0)
    data_am = pd.read_csv(paths['am_path'], index_col=0)
    data_esm = pd.read_csv(paths['esm_path'], index_col=0)

    # Whether to plot the heatmaps standalone or embedded
    created_locally = axes is None
    if created_locally:
        fig, axes = plt.subplots(3, 1, figsize=figsize or (len(data_diff.columns) / 5, 18), sharex=True)
    else:
        fig = None
        if len(axes) != 3:
            raise ValueError("Three axes are required: one each for AlphaMissense, ESM, and Difference heatmaps.")
    
    # Unpack axes
    ax_am, ax_esm, ax_diff = axes

    # --- AlphaMissense Heatmap ---
    hm_am = sns.heatmap(data_am, ax=ax_am, cmap="RdBu_r", vmin=0, vmax=1, cbar=add_colorbars and created_locally,
        cbar_kws={'shrink': 1, 'aspect': 5, 'pad': 0.01} if add_colorbars else None
    )
    ax_am.set_ylabel("AlphaMissense Rank", fontsize=28, fontweight='bold')
    ax_am.set_yticks([i + 0.5 for i in range(len(data_diff.index))])  # Ensure all y-ticks are shown
    ax_am.set_yticklabels(data_am.index, fontsize=8)  # Set labels for y-ticks
    ax_am.tick_params(axis='x', labelbottom=False)

    # --- ESM Heatmap ---
    hm_esm = sns.heatmap(data_esm, ax=ax_esm, cmap="viridis", vmin=0, vmax=1, cbar=add_colorbars and created_locally,
        cbar_kws={'shrink': 1, 'aspect': 5, 'pad': 0.01} if add_colorbars else None
    )
    ax_esm.set_ylabel("ESM-1b Rank", fontsize=28, fontweight='bold')
    ax_esm.set_yticks([i + 0.5 for i in range(len(data_diff.index))])  # Ensure all y-ticks are shown
    ax_esm.set_yticklabels(data_esm.index, fontsize=8)  # Set labels for y-ticks
    ax_esm.tick_params(axis='x', labelbottom=False)

    # --- Difference Heatmap ---
    max_diff = max(abs(data_diff.max().max()), abs(data_diff.min().min()))
    hm_diff = sns.heatmap(data_diff, ax=ax_diff, cmap="BrBG", vmin=-max_diff, vmax=max_diff, 
                          center=0, cbar=add_colorbars and created_locally, 
                          cbar_kws={'shrink': 1, 'aspect': 5, 'pad': 0.01} if add_colorbars else None
    )
    ax_diff.set_ylabel("Difference Rank", fontsize=28, fontweight='bold')
    ax_diff.set_yticks([i + 0.5 for i in range(len(data_diff.index))])  # Ensure all y-ticks are shown
    ax_diff.set_yticklabels(data_diff.index, fontsize=8)  # Set labels for y-ticks
    ax_diff.set_xticks([i + 0.5 for i in range(len(data_diff.columns))])
    ax_diff.set_xticklabels(data_diff.columns, rotation=90, fontsize=8 if created_locally else 5)  # Set labels for x-ticks

    # --- Custom Colorbar Labels ---
    if add_colorbars and created_locally:
        cbar_am = hm_am.collections[0].colorbar
        cbar_am.ax.text(1.5, 0.5, "AlphaMissense Rank", rotation=90, va="center", ha="left",
                        transform=cbar_am.ax.transAxes, color='black', fontsize=30, fontweight='bold')

        cbar_esm = hm_esm.collections[0].colorbar
        cbar_esm.ax.text(1.5, 0.5, "ESM1b Rank", rotation=90, va="center", ha="left",
                         transform=cbar_esm.ax.transAxes, color='black', fontsize=30, fontweight='bold')

        cbar_diff = hm_diff.collections[0].colorbar
        cbar_diff.ax.text(1.5, 0.5, "Difference", rotation=90, va="center", ha="left",
                          transform=cbar_diff.ax.transAxes, color='black', fontsize=30, fontweight='bold')

    # --- Layout ---
    if created_locally:
        plt.tight_layout()
        if save:
            images_path = get_paths()['images_path'] / '5.6.Combo_Heatmap_Average_Pathogenicity'
            output_path = images_path / "heatmap.png"
            plt.savefig(output_path, dpi=300, bbox_inches='tight')

            relative_output = os.path.relpath(output_path, start=get_paths()["project_root"])
            print(f"✅ Heatmap saved to '{relative_output}'.")
    
        plt.show()

    print(f"Maximum difference value: {max_diff}")
    return fig, axes
