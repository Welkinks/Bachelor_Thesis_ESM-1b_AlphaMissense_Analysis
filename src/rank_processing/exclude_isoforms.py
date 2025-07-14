import os
import shutil
from glob import glob
from tqdm import tqdm
from src.project_config import get_paths

def exclude_isoforms(raw, output):
    """
    Copies all CSV files from `src_dir` to `dest_dir` that do not contain a hyphen in their filename.
    
    Parameters:
        raw (str or Path): Directory containing the original CSV files.
        output (str or Path): Target directory where non-isoform CSVs will be copied.
    """
    os.makedirs(output, exist_ok=True)

    # Get all CSV files from source directory
    all_csv_files = glob(os.path.join(raw, "*.csv"))

    # Filter: exclude files with '-' in the filename
    filtered_files = [f for f in all_csv_files if "-" not in os.path.basename(f)]

    # Copy the filtered files
    for file_path in tqdm(filtered_files, desc="Copying non-isoform CSVs"):
        filename = os.path.basename(file_path)
        dest_path = os.path.join(output, filename)
        shutil.copy(file_path, dest_path)

    relative_output = os.path.relpath(output, start=get_paths()["project_root"])
    print(f"✅ Done: copied {len(filtered_files)} files to '{relative_output}'.")
