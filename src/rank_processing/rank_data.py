import os
import pandas as pd
import numpy as np
from scipy.stats import rankdata
from glob import glob
from tqdm import tqdm
from src.project_config import get_paths

def rank_data(input_dir: str, output_dir: str, model: str, intersection_csv_path: str):
    """
    Ranks AlphaMissense pathogenicity scores and outputs per-protein ranked pivot tables.
    
    Parameters:
        input_dir (str): Directory containing CSV files with raw AlphaMissense data.
        output_dir (str): Directory where ranked output CSVs will be saved.
        model (str): Model name used for the ranking process
    """
    # Standard amino acids order
    aa_list = [
        "A", "V", "L", "I", "M", "F", "W",  # Hydrophobic
        "S", "T", "N", "Q", "Y", "C",      # Polar uncharged
        "K", "R", "H",                     # Positively charged
        "D", "E",                          # Negatively charged
        "G", "P"                           # Special
    ]

    # Intersection proteins available in both datasets
    csv_files = pd.read_csv(intersection_csv_path, usecols=["Protein_ID"])

    if model == "AlphaMissense":
        # First pass: collect all scores
        print("Collecting AM scores...")
        all_scores = []

        for protein_id in tqdm(csv_files["Protein_ID"], desc="Reading scores"):
            csv_file_input_path = os.path.join(input_dir, f"{protein_id}.csv")
            df = pd.read_csv(csv_file_input_path, usecols=["am_pathogenicity"])
            all_scores.extend(df["am_pathogenicity"].astype("float32").values)
        
        # Second pass: compute global ranks and write per-protein outputs
        print("Computing global ranks...")
        all_scores = np.array(all_scores, dtype="float32")
        ranks = rankdata(all_scores, method="average") 
        normalized_ranks = np.round(ranks / len(all_scores), 5)

        # Prepare output directory
        print("Processing per-file and writing outputs...")
        score_index = 0
        os.makedirs(output_dir, exist_ok=True)

        # Columns to keep in AlphaMissense raw CSV files for ranking
        columns_to_keep = [
        "uniprot_id", "variation", "am_pathogenicity",
        "residue", "residue_position"]

        for protein_id in tqdm(csv_files["Protein_ID"], desc="Processing files"):
            csv_file_input_path = os.path.join(input_dir, f"{protein_id}.csv")
            
            df = pd.read_csv(csv_file_input_path, usecols=columns_to_keep)

            n_rows = len(df)
            df["rank_score"] = normalized_ranks[score_index:score_index + n_rows]
            score_index += n_rows

            df["pos_label"] = df["residue"] + " " + df["residue_position"].astype(str)

            for uniprot_id, group in df.groupby("uniprot_id"):
                pivot = group.pivot(index="variation", columns="pos_label", values="rank_score")
                pivot = pivot.reindex(aa_list)

                # Sort columns like "M 1", "R 2", ...
                try:
                    pivot = pivot[sorted(pivot.columns, key=lambda x: int(x.split()[1]))]
                except Exception:
                    pass  # In case columns are missing or misformatted

                output_file = os.path.join(output_dir, f"{uniprot_id}_rank.csv")
                pivot.to_csv(output_file)

        print("AlphaMissense Scores were successfully ranked.")



    elif model == "ESM":
        # First pass: collect all LLR values
        print("Collecting LLR scores...")
        all_scores = []
        file_map = []
        
        for protein_id in tqdm(csv_files["Protein_ID"], desc="Reading LLR scores"):
            file_path = os.path.join(input_dir, f"{protein_id}_LLR.csv")

            df = pd.read_csv(file_path, index_col=0)
            flat_scores = df.values.flatten()
            valid_mask = (~np.isnan(flat_scores)) & (flat_scores != 0)  # Exclude NaNs and zeros
            all_scores.extend(flat_scores[valid_mask])
            file_map.append((file_path, df))


        # Global ranking
        print("Computing global ranks...")
        all_scores = np.array(all_scores, dtype="float32")
        ranks = rankdata(-all_scores, method="average")
        normalized_ranks = np.round(ranks / len(all_scores), 5)

        # Second pass: rebuild each matrix with rank scores
        print("Writing ranked matrices...")
        score_index = 0
        
        for file, df in tqdm(file_map):
            flat_llrs = df.values.flatten()
            flat_ranks = np.full_like(flat_llrs, np.nan, dtype="float32")

            valid_mask = (~np.isnan(flat_llrs)) & (flat_llrs != 0)
            flat_ranks[valid_mask] = normalized_ranks[score_index:score_index + valid_mask.sum()]
            score_index += valid_mask.sum()

            rank_matrix = pd.DataFrame(
                flat_ranks.reshape(df.shape),
                index=df.index,
                columns=df.columns
            )

            # Optional: reorder rows to standard AA list
            rank_matrix = rank_matrix.reindex(aa_list)

            outname = os.path.basename(file).replace(".csv", "_rank.csv")
            outpath = os.path.join(output_dir, outname)
            rank_matrix.to_csv(outpath, float_format="%.5f")

        print("ESM1b Scores were successfully ranked.")