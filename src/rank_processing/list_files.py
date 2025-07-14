import os
from src.project_config import get_paths

def list_csv_files(folder_path, suffix, output_file, column_name):
    """
    Lists all .csv files in the specified folder with the given suffix and writes their base names to an output file.

    Args:
        folder_path (str): Path to the folder to search for .csv files.
        suffix (str): Suffix to filter .csv files (e.g., "_statistics.csv").
        output_file (str): Path to the output .csv file to save the base names.
    """
    base_names = []

    # Iterate through files in the folder
    for file_name in os.listdir(folder_path):
        if file_name.endswith(suffix) and file_name.endswith(".csv"):
            # Extract base name by removing the suffix
            base_name = file_name.replace(suffix, "")
            base_names.append(base_name)

    # Write base names to the output file
    with open(output_file, 'w') as f:
        f.write(column_name + '\n')  # Write the column name as the header
        for name in base_names:
            f.write(name + '\n')

    relative_output = os.path.relpath(output_file, start=get_paths()["project_root"])
    print(f"✅ Done: {len(base_names)} files listed in '{relative_output}'.")

# Example usage
# list_csv_files_with_suffix('/path/to/folder', '_statistics.csv', '/path/to/output.csv')