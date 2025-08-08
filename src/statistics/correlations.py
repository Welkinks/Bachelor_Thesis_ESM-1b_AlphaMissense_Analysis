# --- Imports ---
import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from src.project_config import AM_PATH, ESM_PATH, DSSP_PROTEINS_PATH, N_OUT_PROTEIN_ID, PROTEIN_IDS_CSV, MULTISPAN_PROTEIN_ID, PROCESSED_DIR
from tqdm.notebook import tqdm




def load_variant_matrices(am_path, esm_path):
    """
    Load AlphaMissense and ESM matrices and ensure they have the same shape.
    Both matrices must be 20 rows (variants) x protein length columns.
    Returns:
        (am_matrix, esm_matrix) if shapes match
        None if shapes do not match
    """

    try:
        am = pd.read_csv(am_path, index_col=0)  # AM scores are tab-separated
        esm = pd.read_csv(esm_path, index_col=0)  # ESM scores are comma-separated

    
        if am.shape != esm.shape:
            print(f"Shape mismatch: {am.shape} != {esm.shape}")
            return None
        
        length = esm.shape[1]  

        am_flat = am.values.flatten()
        esm_flat = esm.values.flatten()

        # Create a mask where BOTH AM and ESM are not NaN
        valid_mask = ~np.isnan(am_flat) & ~np.isnan(esm_flat)

        # Return aligned flattened arrays
        return am_flat[valid_mask], esm_flat[valid_mask], length  # return as numpy arrays for correlation calculation


    except Exception as e:
        print(f"Error loading files {am_path} and {esm_path}: {e}")
        return None

def load_residue_stats(path, model, corr_with, struct_codes=None, wanted_structs=None):
    """
    Load per-residue stats CSV and return arrays for correlation:
    - model_mean: "AM_mean" or "ESM_mean"
    - rASA: float values
    - structure: one of the raw DSSP codes
    """
    df = pd.read_csv(path)
    length_protein = df["residue_position"].max()

    if model == "AM_vs_ESM_Per_Residue":
        x = df["AM_mean"].values
        y = df["ESM_mean"].values
    else:
        x = df[f"{model}_mean"].values


    if corr_with == "rASA":
        y = df["rASA"].values
    elif corr_with == "secondary_structure":
        # map DSSP codes to full names, then filter to wanted_structs
        df["structure_full"] = df["2D_Structures"].map(struct_codes)
        mask = df["structure_full"].isin(wanted_structs)
        
        
        if model == "AM_vs_ESM_Per_Residue":
            x = x[mask]
            y = y[mask]
        else:
            x = x[mask]
            y = df["rASA"][mask]  # use rASA as dummy, or could map another property

    else:
        raise ValueError("corr_with must be 'rASA' or 'secondary_structure'")
    valid = ~np.isnan(x) & ~np.isnan(y)
    return x[valid], y[valid], length_protein

def correlate_vectors(x, y, method="pearson"):
    """
    Compute correlation between two 1D arrays.
    """
    if len(x) < 2:
        return np.nan
    if method == "pearson":
        return pearsonr(x, y)[0]
    elif method == "spearman":
        return spearmanr(x, y)[0]
    else:
        raise ValueError("method must be 'pearson' or 'spearman'")

def compute_correlations(
    proteome, 
    model="both", 
    corr_range="global", 
    correlation_with=None,
    structures=None,
    output_folder=None,
    method="pearson"
):
    """
    Main function to compute correlations.
    
    Arguments:
    - proteome: list of proteome names ["human_proteome", "N_out", "multispan"]
    - model: "both", "AM", or "ESM"
    - corr_range: "global" or "per_protein"
    - correlation_with: None if model="both"; else "rASA" or "secondary_structure"
    - structures: list of structure names to include if correlation_with="secondary_structure"
    - id_csv: dict mapping proteome -> path to CSV of protein IDs
    - id_col: column name in ID CSV
    - folders: dict with keys:
        - "AM_scores", "ESM_scores" for variant matrices
        - "stats" for per-residue stats
    - output_folder: path to save per-protein CSVs
    
    Outputs:
    - Prints global correlations
    - Saves per-protein CSVs if corr_range="per_protein"
    """
    struct_codes = {
        "H": "Alpha helix", "G": "Helix-3", "I": "Helix-5", "P": "κ‐helix",
        "E": "Strand", "B": "Beta bridge", "T": "Turn", "S": "Bend", "-": "Disordered"
    }
    
    global count_missmatch
    count_missmatch = 0  # global counter for shape mismatches
    

    # paths for protein IDs
    paths = {"N_out": (N_OUT_PROTEIN_ID, "entry"),
             "human_proteome": (PROTEIN_IDS_CSV,"Protein_ID"),
             "multispan": (MULTISPAN_PROTEIN_ID, "Entry")}


    for p in proteome:
        id_csv = paths[p][0] if p in paths else None
        id_col = paths[p][1] if p in paths else None
        ids = pd.read_csv(id_csv, delimiter="\t" if p == "multispan" else ',')[id_col].tolist()
        flat_x, flat_y = [], []
        per_prot_results = []
        
        for prot in tqdm(ids, desc=f"Processing"):
   
            try:
                if model == "both":
                    paths = (os.path.join(AM_PATH, f"{prot}_rank.csv"),
                             os.path.join(ESM_PATH, f"{prot}_LLR_rank.csv"))
                    matrices = load_variant_matrices(*paths)
                    
                    
                    if matrices is None:
                        continue  # Skip proteins with shape mismatch or load errors
                    
                    am, esm, length_prot = matrices
                    x, y, length_protein = am, esm, length_prot
               
                elif model in ["AM", "ESM", "AM_vs_ESM_Per_Residue"]:
                    # load stats for one model vs. physico-chemical
                    stats_path = os.path.join(DSSP_PROTEINS_PATH, f"{prot}_statistics.csv")
                    x, y, length_protein = load_residue_stats(stats_path, model, correlation_with, struct_codes, structures)
            
                if len(x) != len(y):
                    print("AM and ESM-1b have missmatch.")
                    continue

                # append for global
                flat_x.append(x)
                flat_y.append(y)
               
                
                if corr_range == "per_protein":
                    rho = correlate_vectors(x, y, method)
                    per_prot_results.append((prot, rho, length_protein))
                    
            except Exception as e:
                print(f"Error: {e}")
                # skip proteins with issues
                count_missmatch += 1
                print(f"Protein {prot} ({count_missmatch}) has shape mismatch.")

                continue
        
        # global correlation
        gx = np.concatenate(flat_x)
        gy = np.concatenate(flat_y)
        global_rho = correlate_vectors(gx, gy, method)
        print(f"{p} - {model} global {correlation_with or 'AM_vs_ESM'} rho = {global_rho:.4f}")
        

        output_folder = os.path.join(PROCESSED_DIR, "1.3.Correlations")
        # per-protein CSV
        if corr_range == "per_protein" and output_folder:
            df_out = pd.DataFrame(per_prot_results, columns=["protein_id", "correlation", "length"])
            df_out.to_csv(os.path.join(output_folder, f"{p}_{model}_{correlation_with or 'both'}_{structures}per_protein.csv"), index=False)


