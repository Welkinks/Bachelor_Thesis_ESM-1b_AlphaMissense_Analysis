# --- Imports ---
from src.project_config import get_paths_protein
import pandas as pd

# --- Function to identify clusters of high differences in a protein's difference matrix ---
def identify_clusters(id, matrix_diff, 
                      threshold=0.074,
                      search_window_size = 5):
    

    cluster_dict = {}

    avg_diff = matrix_diff.mean(axis=0)  # shape: (N,)
    protein_length = avg_diff.shape[0]
    
    # Step 1: find all high-scoring windows
    high_score_windows = []
    for i in range(protein_length - search_window_size + 1):
        window = avg_diff[i:i + search_window_size]
        window_avg = window.mean()

        if abs(window_avg) >= threshold:
            high_score_windows.append((i + 1, i + search_window_size))

    
    # Step 2: merge overlapping/adjacent windows
    merged = []
    if high_score_windows:
        current_start, current_end = high_score_windows[0]

        for start, end in high_score_windows[1:]:
            if start <= current_end + 1:  # overlapping or adjacent
                current_end = max(current_end, end)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = start, end
        merged.append((current_start, current_end))  # Add last region


    # Step 3: compute average difference for full merged region
    clusters = []
    for start, end in merged:
        region_avg = round(avg_diff[start - 1 : end].mean(), 4) # back to 0-indexed
        method = "ESM1b" if region_avg > 0 else "AlphaMissense"
        clusters.append({
            "range": f"{start}-{end}",
            "avg_diff": region_avg,
            "more pathogenic": method
        })


    if clusters:
        cluster_dict[id] = clusters
            

    return cluster_dict