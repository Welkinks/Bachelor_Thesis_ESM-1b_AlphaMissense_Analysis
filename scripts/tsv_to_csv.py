import os
import pandas as pd

input_file = '/Users/doma/Documents/Bachelor_Arbeit/Code/data/raw/AlphaMissense_aa_substitutions.tsv'
output_dir = '/Users/doma/Documents/Bachelor_Arbeit/Code/data/raw/AlphaMissense2_csv'

os.makedirs(output_dir, exist_ok=True)

chunk_size = 100000
reader = pd.read_table(input_file, sep='\t', skiprows=3, header=0, chunksize=chunk_size)

chunk_num = 0
for chunk in reader:
    chunk_num += 1
    print(f"Processing chunk {chunk_num}...")

    variants = chunk['protein_variant'].str.extract(r'(?P<residue>[A-Z])(?P<residue_position>\d+)(?P<variation>[A-Z])')
    chunk = chunk.join(variants)

    for uniprot_id, group in chunk.groupby('uniprot_id'):
        filename = os.path.join(output_dir, f"{uniprot_id}.csv")
        write_header = not os.path.exists(filename)
        group.to_csv(filename, index=False, mode='a', header=write_header)

print("All files saved successfully!")
