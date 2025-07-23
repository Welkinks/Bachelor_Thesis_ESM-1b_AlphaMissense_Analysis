import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pickle
from matplotlib.lines import Line2D
from src.project_config import get_paths_protein, RAW_DIR, IMAGES_DIR, COLORS_CPP, COLORS_MODELS, COLORS_SECONDARY

# helper function for plot_pathogenicity()
def _add_structured_legend(ax, highlight, model_labels=None, models="both",
                           span="singlespan", feature_labels=None, region_labels=None,
                           legend_loc='lower left', anchor=(0.0, -0.01), created_locally=False):
    """
    Add a structured legend to an axis, filtering by highlight regions.

    Parameters:
        ax             : matplotlib axis to apply the legend to
        highlight      : list of region/feature names to include
        model_labels   : list of model names (default: AlphaMissense, ESM-1b)
        feature_labels : list of all feature labels (e.g. Transmembrane)
        region_labels  : list of all region labels (e.g. Helix, Beta strand)
        legend_loc     : location of legend (default: 'lower left')
        anchor         : bbox_to_anchor (default: (0.0, 0.0))
    """

    handles, labels = ax.get_legend_handles_labels()
    label_to_handle = dict(zip(labels, handles))

    if models == "both":
        model_labels = model_labels or ['AlphaMissense', 'ESM-1b']
    elif models == "AlphaMissense":
        model_labels = model_labels or ['AlphaMissense', ""]
    elif models == "ESM":
        model_labels = model_labels or ['ESM-1b', ""]


    if span == "singlespan":
        feature_labels = [l for l in (feature_labels or ['Transmembrane', 'Juxtamembrane']) if l in highlight]
        region_labels  = [l for l in (region_labels  or ['Helix', 'Beta strand', 'Disordered']) if l in highlight]
    else:  # multispan
        
        region_labels  = [l for l in (region_labels  or ['Helix', 'Beta strand', 'Disordered']) if l in highlight]

    ordered_labels = model_labels[:]
    if feature_labels:
        ordered_labels += [' '] + feature_labels
    if region_labels:
        ordered_labels += [' '] + region_labels

    ordered_handles = [
        label_to_handle.get(l, Line2D([], [], color='none')) for l in ordered_labels
    ]

    ax.legend(
        ordered_handles, ordered_labels,
        loc=legend_loc,
        bbox_to_anchor=anchor if span == "singlespan" else (0.0, -0.05),
        frameon=False,
        fontsize=7 if created_locally else 18,
        ncol=3 if span == "singlespan" else 2,
        columnspacing=1.2,
    )






# Main function to plot pathogenicity
def plot_pathogenicity(
    protein_id,
    curves="both",
    rank=True,
    span="singlespan",    #"multispan", "singlespan"
    method="mean",
    highlight=None,
    smoothing_window=5,
    show_std=True,
    rASA=False,
    secondary_str_version="dssp",    # "uniprot" or "dssp"
    figsize=(12, 5),
    dpi=300,
    title=None,
    save_path=None,
    fig=None,
    ax=None
):
    """
    Generate a pathogenicity plot for a given protein, optionally highlighting structural regions.

    Parameters
    ----------
    protein_id : str
        UniProt or PDB identifier corresponding to a CSV file named
        '{protein_id}_statistics.csv' in `data_dir`.
    curves : {'AlphaMissense', 'ESM', 'both'}, default 'both'
        Which model(s) to plot:
        - 'AlphaMissense' : plot only AM_smoothed_pathogenicity
        - 'ESM'          : plot only ESM_smoothed_pathogenicity
        - 'both'         : plot both curves
    highlight : list of str or None, default None
        List of regions to highlight as bars. Valid keys:
        ['Helix', 'Beta strand', 'Disordered', 'Transmembrane', 'Juxtamembrane'].
        If None, no highlights are drawn.
    data_dir : str, default '/Users/doma/.../Protein_Statistics'
        Directory containing the input CSV files.
    smoothing_window : int, default 5
        Rolling window size for smoothing the pathogenicity scores.
    figsize : tuple, default (12, 5)
        Figure size in inches (width, height).
    dpi : int, default 300
        Resolution of the figure.
    title : str or None
        Plot title. If None, generated as 'Pathogenicity for {protein_id}'.
    save_path : str or None
        If provided, the figure is saved to this file path (PNG format).
        If None, the plot is not saved.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object.
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes object.

    Raises
    ------
    FileNotFoundError
        If the input CSV file does not exist in `data_dir`.
    ValueError
        If `curves` argument is not one of 'AlphaMissense', 'ESM', or 'both'.
    """

    # define Variables
    TRANSMEMBRANE = 'Transmembrane_mask'
    JUXTAMEMBRANE = 'Juxtamembrane_mask'
    HELIX_MASK = 'Helix_mask'

    """
    # Get paths to correct CSV files
    paths = get_paths_protein(protein_id=protein_id)
    if rank == True:
        rankscore = "rank"
        if span == "singlespan":
            data = paths["n_out_rank_path"]
        else:
            data = paths["multispan_rank_path"]

    else:
        rankscore = "normalized"
        if span == "singlespan":
            data = paths["protein_statistics_path"]
        else:
            data = paths["multispan_rank_path"]
    """


    # Choose correct key for data pathway based on span and secondary structure version
    data_key_map = {
        ("singlespan", "uniprot"): "n_out_rank_path" if rank else "protein_statistics_path",
        ("singlespan", "dssp"): "dssp_protein_path",
        ("multispan", "uniprot"): "multispan_rank_path",
        ("multispan", "dssp"): "dssp_protein_path"
    }
    
    data_key = data_key_map.get((span, secondary_str_version))
    if data_key is None:
        raise ValueError("Invalid combination of span and secondary structure version.")




    # check for correct input for secondary structure version
    if secondary_str_version not in ["dssp", "uniprot"]:
        raise ValueError("secondary_str_version must be 'dssp' or 'uniprot'")
    
    # Get paths and select the path and rankscore
    paths = get_paths_protein(protein_id=protein_id)
    data = paths[data_key]


    # Load data
    if not os.path.exists(data):
        raise FileNotFoundError(f"CSV not found: {data}")
    
    df = pd.read_csv(data, index_col=0 if secondary_str_version == "uniprot" else None)



    # Define x-axis column
    x_col = 'residue_position'
    if x_col not in df.columns:
        raise KeyError(f"Column '{x_col}' missing in CSV.")

    if secondary_str_version == "uniprot":
        # Preprocess unified Helix mask
        if span == "singlespan":
            df[HELIX_MASK] = df['Helix'].fillna(0).eq(1)
            df[JUXTAMEMBRANE] = pd.to_numeric(df['Juxtamembrane'], errors='coerce').fillna(0).eq(1)
        else:
            # For multispan, helices are indicated in Transmembrane as Helical_1, Helical_2, etc.
            df[TRANSMEMBRANE] = df['Transmembrane'].astype(str).str.lower().str.startswith('helical')
            df[HELIX_MASK] = df[TRANSMEMBRANE]  # In multispan, helices are in TM
            df['Transmembrane_Helix_mask'] = df[TRANSMEMBRANE]
            df[JUXTAMEMBRANE] = False  # or drop if unused

        df['Beta strand'] = df['Beta strand'].fillna(0).eq(1)
        df[TRANSMEMBRANE] = df['Transmembrane'].astype(str).str.lower().str.startswith('helical')
        df['Region_mask'] = df['Region'].astype(str).str.lower().eq('disordered')

    elif secondary_str_version == "dssp":
        # Normalize DSSP codes "-" for missing values
        dssp = df["2D_Structures"].fillna("-").astype(str)

        # DSSP to masks
        df[HELIX_MASK] = dssp.isin(["H", "G", "I"])       # Helix-related
        df["Beta strand"] = dssp.isin(["E", "B"])         # Beta-sheet/strand
        df["Region_mask"] = dssp.isin(["-", "C", "S", "T"])  # Coil/loop/disordered

        # For compatibility with highlight logic
        df[TRANSMEMBRANE] = df['Transmembrane'].astype(str).str.lower().str.startswith('helical')  # DSSP doesn't provide this
        df[JUXTAMEMBRANE] = df[JUXTAMEMBRANE] = pd.to_numeric(df['Juxtamembrane'], errors='coerce').fillna(0).eq(1)

        """
        structure_codes = {"H": "Alpha helix",
                           "B": "Beta bridge",
                           "E": "Strand",
                           "G": "Helix-3",
                           "I": "Helix-5",
                           "T": "Turn",
                           "S": "Bend"}

            df["Turn"] = dssp == "T"
            df["Bend"] = dssp == "S"
            df["310 Helix"] = dssp == "G"
        """




    # Smooth scores
    if method == "mean":
        df['AM_smoothed_pathogenicity'] = df['AM_mean'].rolling(
            smoothing_window, center=True, min_periods=1).mean()
        df['ESM_smoothed_pathogenicity'] = df['ESM_mean'].rolling(
            smoothing_window, center=True, min_periods=1).mean()
    elif method == "median":
        df['AM_smoothed_pathogenicity'] = df['AM_median'].rolling(
            smoothing_window, center=True, min_periods=1).mean()
        df['ESM_smoothed_pathogenicity'] = df['ESM_median'].rolling(
            smoothing_window, center=True, min_periods=1).mean()

    # Define curve specs
    all_curves = {
        'AlphaMissense': {'y_col':'AM_smoothed_pathogenicity','label':'AlphaMissense',
                          'color':COLORS_MODELS["AM"],'std_col':'AM_std'},
        'ESM':           {'y_col':'ESM_smoothed_pathogenicity','label':'ESM-1b',
                          'color':COLORS_MODELS["ESM"],'std_col':'ESM_std'}}
    
    

    if curves=='both':
        curve_list = [all_curves['AlphaMissense'], all_curves['ESM']]
        

    elif curves in all_curves:
        curve_list = [all_curves[curves]]
        
    else:
        raise ValueError("curves must be 'AlphaMissense','ESM', or 'both'")


    # Prepare highlights based on span type
    if span == "singlespan":
        all_highlights = {
            'Helix':{'column':HELIX_MASK,'color':COLORS_SECONDARY["BLUE"],'alpha':1,'shade_alpha':0.4,'position':'bottom'},
            'Beta strand':{'column':'Beta strand','color':COLORS_SECONDARY["PURPLE"],'alpha':1,'shade_alpha':0.4,'position':'bottom'},
            'Disordered':{'column':'Region_mask','color':COLORS_SECONDARY["PINK"],'alpha':1,'shade_alpha':0.4,'position':'bottom'},
            'Transmembrane':{'column':TRANSMEMBRANE,'color':COLORS_CPP["TRANS"],'alpha':1,'shade_alpha':0.0,'position':'bottom'},
            'Juxtamembrane':{'column':JUXTAMEMBRANE,'color':COLORS_CPP["JUXTA"],'alpha':1,'shade_alpha':0.0,'position':'bottom'}
        }
        hl_regions = {k:all_highlights[k] for k in (highlight or []) if k in all_highlights}

    else:
        # In multispan avoid Juxtamembrane 
        all_highlights = {
            'Helix':{'column':HELIX_MASK,'color':COLORS_SECONDARY["BLUE"],'alpha':1,'shade_alpha':0.4,'position':'bottom'},
            'Beta strand':{'column':'Beta strand','color':COLORS_SECONDARY["PURPLE"],'alpha':1,'shade_alpha':0.4,'position':'bottom'},
            'Disordered':{'column':'Region_mask','color':COLORS_SECONDARY["PINK"],'alpha':1,'shade_alpha':0.4,'position':'bottom'},
            'Transmembrane':{'column':TRANSMEMBRANE,'color':COLORS_CPP["TRANS"],'alpha':1,'shade_alpha':0.0,'position':'bottom'}
        }
        hl_regions = {k:all_highlights[k] for k in (highlight or []) if k in all_highlights}

    # Initialize plot
    if fig is None or ax is None:
        if not rASA:
            fig, ax = plt.subplots(figsize=figsize, dpi=dpi)    # otherwise used for function combining plot_pathogenicity and heatmap
            ax_rasa = None
        elif rASA:
            fig = plt.figure(figsize=figsize, dpi=dpi, constrained_layout=False)  # Slightly taller for rASA
            gs = gridspec.GridSpec(2, 1, height_ratios=[5, 0.6], hspace=0.05)
            ax = fig.add_subplot(gs[0])
            ax_rasa = fig.add_subplot(gs[1], sharex=ax)

        created_locally = True
    else:
        created_locally = False


    if rASA and created_locally:
        # Plot rASA as a line (or points) on its own axis
        ax_rasa.plot(df[x_col], df['rASA'], color='black', linewidth=0.6, rasterized=True)
        ax_rasa.set_ylabel('rASA', fontsize=8, labelpad=2, fontweight='bold')
        ax_rasa.set_ylim(0, 1.01)
        ax_rasa.tick_params(axis='x', labelsize=8)
        ax_rasa.tick_params(axis='y', labelsize=8)
        ax_rasa.set_yticks([0, 0.5, 1.0])
        ax_rasa.grid(False)
        ax_rasa.set_xlabel("Residue Position", fontsize=20 if created_locally==False else 9, fontweight='bold')
        
        # Hide redundant x-axis on top plot
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(axis='x', which='both', bottom=False, top=False)







    # Set x-axis and y-axis labels and title
    if rASA == False:
        ax.set_xlabel("Residue Position", fontsize=20 if created_locally==False else 9, fontweight='bold')

    # set y-axis label 
    ylabel_base = "Mean Predicted Pathogenicity Rank" if rank else "Mean Predicted Pathogenicity"
    smoothing_info = f"\n(Smoothing: {smoothing_window})" if smoothing_window > 1 else ""
    ax.set_ylabel(f"{ylabel_base}{smoothing_info}",
                  fontsize=20 if not created_locally else 9, fontweight='bold')

    # y-ticks
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])  
    ax.tick_params(axis='y', labelsize=8 if created_locally else 20)  # Set y-tick label size

    # Set title if created locally and title allowed (not None)
    if created_locally and title is not None:
        ax.set_title(f'Pathogenicity Landscape of {protein_id} Predicted by AlphaMissense and ESM-1b', fontsize=10, fontweight='bold')
    
    # Plot curves and optional std
    for spec in curve_list:
        ycol, lbl, clr, sdcol = spec['y_col'], spec['label'], spec['color'], spec['std_col']
        ax.plot(df[x_col], df[ycol], label=lbl, color=clr, linewidth=1.25 if created_locally else 3, rasterized=True)
        if show_std and sdcol in df:
            y, err = df[ycol], df[sdcol]
            ax.fill_between(df[x_col], y-err, y+err, color=clr, alpha=0.2, zorder=0)

    # Highlight regions
    y_min, y_max = 0, 1.005
    bar_h = 0.015
    used = set()
    # Put Juxtamembrane last
    sorted_regions = sorted(hl_regions.items(), key=lambda x: x[0] == 'Juxtamembrane')
    for i, (name, cfg) in enumerate(sorted_regions):
        col = cfg['column']; pos = cfg['position']; clr = cfg['color']
        ab, asha = cfg['alpha'], cfg['shade_alpha']
        mask = df[col].astype(bool).values
        padded = np.r_[False, mask, False]
        idx = np.where(padded[:-1] != padded[1:])[0]
        for s, e in zip(idx[0::2], idx[1::2]):
            xs, xe = df[x_col].iloc[s], df[x_col].iloc[e-1]
            if asha > 0:
                ax.axvspan(xs, xe, color=clr, alpha=asha, zorder=0)  # Still background
            base = y_max - 1.5 * bar_h if pos == 'top' else y_min
            ax.axhspan(
                base, base + bar_h,
                xmin=(xs - df[x_col].min()) / (df[x_col].max() - df[x_col].min()),
                xmax=(xe - df[x_col].min()) / (df[x_col].max() - df[x_col].min()),
                color=clr, alpha=ab,
                label=name if name not in used else None,
                zorder=5 + i  # Ensure later-drawn regions layer on top
            )
        used.add(name)
    # Final styling
    ax.set_xlim(df[x_col].min(), df[x_col].max()); ax.set_ylim(y_min, 1)

    ### Add legend ###
    if span == "singlespan":
        _add_structured_legend(ax, highlight, models=curves, span=span, 
                               created_locally=created_locally if created_locally else False)
    else: 
        _add_structured_legend(ax, highlight, models=curves, span=span,
                               created_locally=created_locally if created_locally else False)

    if created_locally == False:
        # Restore x-ticks explicitly (important if ax is shared)
        ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True) 

    # Optional: show specific x-ticks and y-ticks
    ax.set_xticks(df[x_col][::10 if created_locally==False else 50])  # or ::20 to show fewer labels
    ax.tick_params(axis='x', rotation=45)  # Rotate x-tick labels for better readability

    # gridlines for reference
    ax.grid(True, which='major', axis='y', linestyle='--', alpha=0.25, zorder=0)
    ax.grid(True, which='major', axis='x', linestyle='--', alpha=0.25, zorder=0) 
    ax_rasa.grid(True, which='major', axis='y', linestyle='--', alpha=0.25, zorder=0) if rASA else None
    ax_rasa.grid(True, which='major', axis='x', linestyle='--', alpha=0.25, zorder=0) if rASA else None

    # Adjust the text size of the x-tick labels
    ax.set_xticklabels(ax.get_xticks(), fontsize=18 if created_locally==False else 8)  # Set fontsize to desired value

    # Save settings
    rankscore = "rank" if rank else "normalized"    # for file naming purposes
    if save_path:
        if show_std:        
            save_path_1 = os.path.join(IMAGES_DIR / "5.3.Visualization",
                                        f"{protein_id}_{curves}_{rankscore}_{method}std_plot.png")
        else:
            save_path_1 = os.path.join(IMAGES_DIR / "5.3.Visualization",
                                        f"{protein_id}_{curves}_{rankscore}_{method}plot.png")

    # --- Layout --- only if created here
    if created_locally:
        #plt.tight_layout()
        fig.subplots_adjust(hspace=0.05, top=0.95, bottom=0.12, left=0.08, right=0.98)
        fig.savefig(save_path_1, dpi=dpi, bbox_inches='tight')  # High quality
        plt.show()

    return fig, ax
