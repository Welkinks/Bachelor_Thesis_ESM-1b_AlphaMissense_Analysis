from src.project_config import get_paths, get_paths_protein, get_aa_list
import pandas as pd
from tqdm import tqdm
import numpy as np
from os import mkdir
import pandas as pd
import numpy as np
from tqdm import tqdm

def compute_variant_bias(dataset):
    """
    Compute ESM and AM 20x20 substitution matrices for a given dataset and their average pathogenicity 
    scores for both models (ESM and AlphaMissense).

    Args:
        dataset (str): The dataset to process. Options are "human_proteome", "N_out_proteome", "Multispan_proteome".

    Returns:
        None: Saves the matrices as CSV files.
    """
    # Load protein IDs based on the dataset
    if dataset == "human_proteome":
        protein_ids = pd.read_csv(get_paths()["protein_ids_intersection"])
        output_suffix = "human_proteome"
        correct_entry = "Protein_ID"
    elif dataset == "N_out_proteome":
        protein_ids = pd.read_csv(get_paths()["processed"] / "N_out_proteome_cleaned.csv")
        output_suffix = "N_out_proteome"
        correct_entry = "entry"
    elif dataset == "Multispan_proteome":
        protein_ids = pd.read_csv(get_paths()["processed"] / "Multispan_proteome_cleaned.csv")
        output_suffix = "multispan_proteome"
        correct_entry = "Entry"
    else:
        raise ValueError("Invalid dataset. Choose from 'human_proteome', 'N_out_proteome', or 'Multispan_proteome'.")

    # Get amino acid list
    aa_list = get_aa_list()

    # Initialize matrices
    esm_matrix = pd.DataFrame(0.0, index=aa_list, columns=aa_list)
    am_matrix = pd.DataFrame(0.0, index=aa_list, columns=aa_list)

    # Loop over each variant amino acid
    for variant in tqdm(aa_list, desc="Outer loop over variant AAs"):
        esm_all_variant_scores = []
        am_all_variant_scores = []

        esm_ref_scores_by_aa = {ref_aa: [] for ref_aa in aa_list}
        am_ref_scores_by_aa = {ref_aa: [] for ref_aa in aa_list}

        # Inner loop over proteins
        for protein in tqdm(protein_ids[correct_entry], desc=f"Inner loop over proteins. Variant: {variant}", leave=False):
            # Load data
            esm_file = pd.read_csv(get_paths_protein(protein)["esm_path"], index_col=0)
            am_file = pd.read_csv(get_paths_protein(protein)["am_path"], index_col=0)

            # Extract the row for the variant
            esm_variant_row = esm_file.loc[variant]
            am_variant_row = am_file.loc[variant]

            # Flatten and clean scores
            esm_values = esm_variant_row.dropna().values
            am_values = am_variant_row.dropna().values
            esm_all_variant_scores.extend(esm_values)
            am_all_variant_scores.extend(am_values)

            # Iterate over each position (column) to extract the reference AA
            for col in esm_variant_row.index:
                if pd.isna(esm_variant_row[col]):
                    continue
                try:
                    ref_aa = col.split()[0]  # e.g., "M 1" -> "M"
                    if ref_aa in aa_list:
                        esm_ref_scores_by_aa[ref_aa].append(esm_variant_row[col])
                except Exception:
                    continue

            for col in am_variant_row.index:
                if pd.isna(am_variant_row[col]):
                    continue
                try:
                    ref_aa = col.split()[0]
                    if ref_aa in aa_list:
                        am_ref_scores_by_aa[ref_aa].append(am_variant_row[col])
                except Exception:
                    continue

        # Compute average variant score and fill diagonal with 0
        esm_avg_scores = {
            ref: np.mean(scores) if ref != variant and scores else 0.0
            for ref, scores in esm_ref_scores_by_aa.items()
        }

        am_avg_scores = {
            ref: np.mean(scores) if ref != variant and scores else 0.0
            for ref, scores in am_ref_scores_by_aa.items()
        }

        # Fill matrices
        for ref_aa in aa_list:
            esm_matrix.loc[variant, ref_aa] = esm_avg_scores.get(ref_aa, 0.0)
            am_matrix.loc[variant, ref_aa] = am_avg_scores.get(ref_aa, 0.0)

    # Save the matrices
    output_dir = get_paths()["processed"] / "5.6.Bias"
    output_dir.mkdir(exist_ok=True)
    esm_matrix.to_csv(output_dir / f"esm_20x20_matrix_{output_suffix}.csv")
    am_matrix.to_csv(output_dir / f"am_20x20_matrix_{output_suffix}.csv")
    print(f"Saved ESM and AM 20x20 matrices for {dataset}.")