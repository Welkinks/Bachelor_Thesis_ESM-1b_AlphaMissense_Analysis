# --- Project Setup ---
from setup_notebook import setup_project_root
setup_project_root()


# --- Project Setup ---
import os 
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import requests
from matplotlib.colors import LinearSegmentedColormap
from src.project_config import get_paths



# For every file in the esm_path folder, compute the difference with the corresponding file in am_path
# Every CSV file is a matrix and both models have exactly aligned cells

def compute_difference(esm_path, am_path, difference_path):
    esm_files = {file.replace("_LLR_rank.csv", "") for file in os.listdir(esm_path)}
    am_files = {file.replace("_rank.csv", "") for file in os.listdir(am_path)}

    count = 0
    mismatch = []

    # Find common files in both directories
    common_files = sorted(esm_files.intersection(am_files))
    print(f"Common files: {len(common_files)}")
    
    for file in tqdm(common_files, total=len(common_files)):
        esm_file = f"{file}_LLR_rank.csv"
        am_file = f"{file}_rank.csv"

        esm_df = pd.read_csv(esm_path / esm_file, index_col=0)
        am_df = pd.read_csv(am_path / am_file, index_col=0)

        # Ensure both dataframes have the same shape
        if esm_df.shape != am_df.shape:
            count += 1
            tqdm.write(f"Shape mismatch for {file} - {count} Skipping.")
            mismatch.append(file)
            continue

        # Compute the difference and round to 5 decimal places
        difference_df = (esm_df - am_df).round(5)

        # Save the difference dataframe to a new CSV file
        difference_filename = f"{file}_rank_difference.csv"
        difference_df.to_csv(difference_path / difference_filename)


    mismatch = pd.DataFrame(mismatch)
    mismatch.to_csv(difference_path / "000mismatch_files.csv", index=False, header=False)
    print("Difference computation completed for all common files.") 




# Ensure the script runs only when executed directly
if __name__ == "__main__":
    # Run the compute_difference function
    paths = get_paths()
    esm_path = paths["esm_path"]
    am_path = paths["am_path"]
    difference_path = paths["difference_path"]

    compute_difference(esm_path, am_path, difference_path)