# --- Project Setup ---
from setup_notebook import setup_project_root
setup_project_root()
import os 
import pandas as pd
from src.project_config import DATA_DIR
from tqdm import tqdm


input_file = DATA_DIR / "raw/AlphaMissense_aa_substitutions.tsv"
output_dir = DATA_DIR / "raw/AlphaMissense_csv"
os.makedirs(output_dir, exist_ok=True)

chunk_size = 100000
reader = pd.read_table(input_file, sep='\t', skiprows=3, header=0, chunksize=chunk_size)

chunk_num = 0
for chunk in tqdm(reader, desc="Processing chunks"):
    chunk_num += 1
    print(f"Processing chunk {chunk_num}...")

    variants = chunk['protein_variant'].str.extract(r'(?P<residue>[A-Z])(?P<residue_position>\d+)(?P<variation>[A-Z])')
    chunk = chunk.join(variants)

    for uniprot_id, group in chunk.groupby('uniprot_id'):
        filename = os.path.join(output_dir, f"{uniprot_id}.csv")
        write_header = not os.path.exists(filename)
        group.to_csv(filename, index=False, mode='a', header=write_header)


print(f"All files saved successfully! Folder: {output_dir}")




