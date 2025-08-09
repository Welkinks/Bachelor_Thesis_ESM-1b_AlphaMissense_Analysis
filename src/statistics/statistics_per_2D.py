# --- Imports ---
from src.project_config import get_paths_protein, get_paths
import pandas as pd










def compute_2D_statistics(protein_id, feature):
    """
    Computes statistics for 2D structures of a given protein ID.
    
    Args:
        protein_id (str): The protein ID to compute statistics for.
        feature (str): Whether Topology or secondary structure (2D) statistics 
        
    Returns:
        pd.DataFrame: A DataFrame containing the average rank score per 2D structure and count occurrences.
    """

    # mapping of 2D structure codes to their descriptions from the newest DSSP
    structure_codes = {"H": "Alpha helix",
                "B": "Beta bridge",
                "E": "Strand",
                "G": "Helix-3",
                "I": "Helix-5",
                "T": "Turn",
                "S": "Bend",
                "P": "κ‐helix",
                "-": "Disorded"}

    # Check if correct arguments
    if feature not in ["Topology", "2D"]:
        raise ValueError("Feature must be either 'Topology' or '2D'.")

    # Load Rank Score statistics for a protein
    protein_statistics = pd.read_csv(get_paths_protein(protein_id)["dssp_protein_path"], index_col=0)

    # Compute average rank score per 2D_Structure and count occurrences
    statistics = (
        protein_statistics.groupby("2D_Structures" if feature == "2D" else "Topological domain")
        .agg(
            AM_mean=("AM_mean", "mean"),
            AM_std=("AM_mean", "std"),
            ESM_mean=("ESM_mean", "mean"),
            ESM_std=("ESM_mean", "std"),
            count=("AM_mean", "count")  
        )
        .reset_index()
        .sort_values("count", ascending=True)
    )




    if feature == "2D":
        statistics["Structure"] = statistics["2D_Structures"].map(structure_codes)

    else: # feature == "Topology" - to include Transmembrane statistics
        # Preprocess the "Transmembrane" column to group all "Helical" variants as "Helical"
        protein_statistics["Transmembrane"] = protein_statistics["Transmembrane"].astype(str)
        protein_statistics["Transmembrane"] = protein_statistics["Transmembrane"].apply(
            lambda x: "Helical" if x.startswith("Helical") else x
        )

        statistics_transmembrane = (
                protein_statistics.groupby("Transmembrane")
                .agg(
                    AM_mean=("AM_mean", "mean"),
                    AM_std=("AM_mean", "std"),
                    ESM_mean=("ESM_mean", "mean"),
                    ESM_std=("ESM_mean", "std"),
                    count=("AM_mean", "count")  
                )
                .reset_index()
                .sort_values("count", ascending=True)
                )

        # Locate the row where "Transmembrane" is "Helical"
        helical_row = statistics_transmembrane.loc[statistics_transmembrane["Transmembrane"] == "Helical"]
        
        if not helical_row.empty:
            # Access specific values from the row
            am_mean = helical_row["AM_mean"].values[0]
            am_std = helical_row["AM_std"].values[0]
            esm_mean = helical_row["ESM_mean"].values[0]
            esm_std = helical_row["ESM_std"].values[0]
            count = helical_row["count"].values[0]

            # Assign these values to the statistics DataFrame
            statistics = pd.concat([statistics,
                pd.DataFrame([{
                    "Topological domain": "Transmembrane",
                    "AM_mean": am_mean,
                    "AM_std": am_std,
                    "ESM_mean": esm_mean,
                    "ESM_std": esm_std,
                    "count": count,
                }])],
                ignore_index=True,
            )



    return statistics


