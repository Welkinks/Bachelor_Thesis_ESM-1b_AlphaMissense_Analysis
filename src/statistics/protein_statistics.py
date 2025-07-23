# --- Imports ---
import os 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm
import requests
import re
import time
import requests
import itertools
import pickle
import traceback
from src.project_config import PROJECT_ROOT, RAW_DIR, PROCESSED_DIR, ESM_PATH, AM_PATH, PROTEIN_IDS_CSV


# --- Helper Functions --- #

def safe_get_uniprot_json(uniprot_id, retries=5, delay=3):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, timeout=10)
            r.raise_for_status()
            return r.json()
        except (requests.exceptions.SSLError, requests.exceptions.ConnectionError) as e:
            print(f"[Retry {attempt}/{retries}] SSL error for {uniprot_id}: {e}")
            time.sleep(delay)
        except requests.exceptions.RequestException as e:
            print(f"[Request error] {uniprot_id}: {e}")
            break
    print(f"[FAIL] Could not fetch UniProt data for {uniprot_id}")
    return None



def get_structural_features(uniprot_id, feature_types):
    
    if feature_types is None:
        return ("Feature Types: 'Turn', 'Glycosylation', 'Binding site', 'Motif', "
                "'Sequence conflict', 'Signal', 'Site', 'Modified residue', 'Helix', "
                "'Natural variant', 'Domain', 'Mutagenesis', 'Cross-link', 'Beta strand', "
                "'Disulfide bond', 'Alternative sequence', 'Region', 'Compositional bias', "
                "'Transmembrane', 'Topological domain', 'Peptide', 'Chain'")


    


    data = safe_get_uniprot_json(uniprot_id)
    if data is None:
        return [], None  # fail safely

    features = data.get('features', [])
    feature_list = []

    for feature in features:
        if feature['type'] in feature_types:
            if 'end' in feature['location']:
                start = int(feature['location']['start']['value'])
                end = int(feature['location']['end']['value'])
            else:
                start = end = int(feature['location']['start']['value'])
            feature_list.append({
                'type': feature['type'],
                'description': feature.get('description', ''),
                'start': start,
                'end': end
            })

    # sequence length with fallback
    seq_length = None
    seq_info = data.get('sequence')
    if seq_info:
        if 'length' in seq_info:
            try:
                seq_length = int(seq_info['length'])
            except (ValueError, TypeError):
                pass
        elif 'value' in seq_info:
            try:
                seq_length = len(seq_info['value'])
            except Exception:
                pass

    return feature_list, seq_length




# Add structural features 
# trying add juxtamembrane positions 
def map_features_to_residues_multi(
        seq_length,
        feature_list,
        feature_types,
        jmd_pad=10        # ±-window around each TM span
    ):

    """
    Build a per-residue annotation frame.
    • data accessed from UniProt

    • feature_types are copied 1-to-1 into their own columns  
    • a new 'Juxtamembrane' column is filled with '1' for
      residues lying jmd_pad residues just outside any TM span
    """
    # Check if feature_types is not None
    if feature_types is None:
        print('Feature Types: "Region", "Site", "Helix", "Beta strand", "Transmembrane", "Modified residue", "Glycosylation"')


    # ➊ allocate empty strings for each feature column
    annot_dict = {ft: [''] * seq_length for ft in feature_types}
    jmd        = [''] * seq_length      # new column

    # ➋ first pass – fill normal feature columns and collect TM ranges
    tm_ranges = []                      # (start,end) tuples, 1-based
    for feat in feature_list:
        ftype = feat['type']
        if ftype not in feature_types:
            continue
        start, end = feat['start'], feat['end']
        label      = feat.get('description', '') or '1'

        for pos in range(start, end + 1):
            if 1 <= pos <= seq_length:
                idx = pos - 1
                annot_dict[ftype][idx] = (
                    f"{annot_dict[ftype][idx]},{label}" if annot_dict[ftype][idx] else label
                )

        if ftype == "Transmembrane":
            tm_ranges.append((start, end))

    # ➌ second pass – flag Juxtamembrane (±10 outside each TM core)
    for s, e in tm_ranges:
        for pos in range(max(1, s - jmd_pad), s):           # before TM
            jmd[pos - 1] = '1'
        for pos in range(e + 1, min(seq_length, e + jmd_pad) + 1):  # after TM
            jmd[pos - 1] = '1'

    # ➍ build the output DataFrame
    df = pd.DataFrame({'residue_position': range(1, seq_length + 1)})
    for ftype in feature_types:
        df[ftype] = annot_dict[ftype]
    df['Juxtamembrane'] = jmd          # <- new column

    return df



def annotate_top_variants(seq_len, top_df):
    annot = [''] * seq_len
    for rank, row in enumerate(top_df.itertuples(), 1):
        annot[row.residue_position - 1] += f"{rank}{row.variation},"
    annot = [a.rstrip(',') for a in annot]
    return annot



def matrix_transformation_df(esm1b_path, model="am_score"):
    df = pd.read_csv(esm1b_path, index_col=0)

    # Clip all values in the DataFrame to [0, 1] before reshaping
    df = df.clip(lower=0, upper=1)

    # Columns: e.g. 'M 1', 'L 2', ...
    columns = [re.match(r'([A-Z])\s+(\d+)', c) for c in df.columns]
    wt_res = [m.group(1) for m in columns]
    positions = [int(m.group(2)) for m in columns]
    tidy = []
    for col, wt, pos in zip(df.columns, wt_res, positions):
        for mut_aa in df.index:
            value = df.loc[mut_aa, col]
            if pd.isnull(value): continue
            tidy.append({
                'residue_position': pos,
                'residue': wt,
                'variation': mut_aa,
                model: value
            })
    tidy_df = pd.DataFrame(tidy)
    return tidy_df


# --- Main Function --- #


def protein_statistics(store_pathway,  
                       dataset,
                       topology = None
                       ):
    """
    Compute protein-level statistics from pathogenicity score matrices and add structural/topological annotations.
    """

    if topology not in ["singlespan", "multispan", None]:
        raise ValueError("Invalid topology specified. Use 'singlespan', 'multispan', or None.")

    if dataset not in ["N_out", "Multispan", "DSSP", "all_human"]:
        raise ValueError("Invalid dataset specified. Use 'N_out', 'Multispan', 'DSSP' or 'all_human'")

    # Load DSSP dataset
    with open(RAW_DIR / '0000DSSP_result_complete_dict.pkl', 'rb') as f:
        data = pickle.load(f)

    # Define protein IDs based on the dataset
    if dataset == "N_out":
        protein_ids = pd.read_csv(PROCESSED_DIR / "Human_N_Out_Proteome.csv", usecols=["entry"])
    elif dataset == "Multispan":
        protein_ids =  pd.read_csv(RAW_DIR  / "multipass.tsv", usecols=["Entry"], sep='\t')
    elif dataset == "DSSP":
        protein_ids = data.keys()
        dssp_human_ids = PROCESSED_DIR / "Protein_IDs_Per_Experiment" / "DSSP_human_protein_ids.csv"
        protein_ids = set(protein_ids).intersection(set(pd.read_csv(dssp_human_ids, usecols=["entry"])['entry']))
    elif dataset == "all_human":
        protein_ids = pd.read_csv(PROTEIN_IDS_CSV, usecols=["Protein_ID"])['Protein_ID'].tolist()

        

    # Special cases
    if dataset == "DSSP" and topology is not None:
        if topology == "singlespan":
            n_out_ids = pd.read_csv(PROCESSED_DIR / "Human_N_Out_Proteome.csv", usecols=["entry"])
            protein_ids = set(protein_ids).intersection(set(n_out_ids['entry']))
        elif topology == "multispan":
            multispan_ids = pd.read_csv(RAW_DIR / "multipass.tsv", usecols=["Entry"], sep='\t')
            protein_ids = set(protein_ids).intersection(set(multispan_ids['Entry']))

    # create store pathway and .csv file for storing misaligned or missing proteins
    os.makedirs(store_pathway, exist_ok=True)
    misaligned = store_pathway / "misaligned_proteins_rank.csv" 
    skipped_proteins = []


    # Main loop over all protein IDs
    for protein in tqdm(protein_ids, desc="Processing proteins", unit="protein"):
        try:
            # Load the .csv file from AM and ESM and check whether both exist
            protein_csv_am = protein + "_rank.csv"
            protein_csv_esm = protein + "_LLR_rank.csv"

            if protein_csv_am not in os.listdir(AM_PATH) or protein_csv_esm not in os.listdir(ESM_PATH):
                print(f"[WARNING] Skipping {protein}: missing rank files.")
                skipped_proteins.append((protein, None, None))
                continue
            
            protein_pathway_am = os.path.join(AM_PATH, protein_csv_am)
            protein_pathway_esm = os.path.join(ESM_PATH, protein_csv_esm)

            df_protein_AM = matrix_transformation_df(protein_pathway_am, model="am_rank_score")
            df_protein_ESM = matrix_transformation_df(protein_pathway_esm, model="esm_rank_score")

            grouped_AM = df_protein_AM.groupby(['residue_position', 'residue'])
            grouped_ESM = df_protein_ESM.groupby(['residue_position', 'residue'])

            stats_AM = grouped_AM['am_rank_score'].agg(['mean', 'median', 'std', 'max', 'min']).round(4).add_prefix('AM_')
            stats_ESM = grouped_ESM['esm_rank_score'].agg(['mean', 'median', 'std', 'max', 'min']).round(4).add_prefix('ESM_')
            stats = stats_AM.merge(stats_ESM, left_index=True, right_index=True).reset_index()

            threshold_AM = 0.7
            threshold_ESM = 0.7
            common_keys = grouped_AM.groups.keys() & grouped_ESM.groups.keys()

            top_mutants_dict_AM = {}
            top_mutants_dict_ESM = {}

            for key in common_keys:
                group_AM = grouped_AM.get_group(key)
                group_ESM = grouped_ESM.get_group(key)

                high_am = group_AM[group_AM['am_rank_score'] > threshold_AM]
                top_mutants_dict_AM[key] = ','.join(high_am.sort_values('am_rank_score', ascending=False)['variation'].astype(str).str.upper())

                high_esm = group_ESM[group_ESM['esm_rank_score'] > threshold_ESM]
                top_mutants_dict_ESM[key] = ','.join(high_esm.sort_values('esm_rank_score', ascending=False)['variation'].astype(str).str.upper())

            stats['key'] = list(zip(stats['residue_position'], stats['residue']))
            stats['AM_top_var'] = stats['key'].map(top_mutants_dict_AM).fillna('')
            stats['ESM_top_var'] = stats['key'].map(top_mutants_dict_ESM).fillna('')
            stats.drop(columns='key', inplace=True) 
            stats = stats.reset_index()

            # UniProt structural features

            feature_types=['Region', 'Topological domain',
                           'Transmembrane', 'Domain', 'Signal']

            seq_length_dataset = int(stats["residue_position"].max())
            features, seq_length_uniprot = get_structural_features(protein,
                                                    feature_types=feature_types)



            if seq_length_uniprot is not None and seq_length_dataset != seq_length_uniprot:
                print(f"[WARNING] Skipping {protein}: seq length mismatch (AM: {seq_length_dataset}, UniProt: {seq_length_uniprot})")
                skipped_proteins.append((protein, seq_length_dataset, seq_length_uniprot))
                continue

            df_annot = map_features_to_residues_multi(seq_length_dataset, 
                                                      features, 
                                                      feature_types=feature_types)
            
            merged = stats.merge(df_annot, on='residue_position', how='left')

            # 2D Structures from DSSP
            merged['2D_Structures'] = merged['residue_position'].apply(
                lambda x: data[protein][x][2] if x in data[protein] else None)
            merged['rASA'] = merged['residue_position'].apply(
                lambda x: data[protein][x][3] if x in data[protein] else None)

            # Top 10 mutations in TMD regions
            tmd_mask = merged['Transmembrane'].astype(bool)
            tmd_positions = set(merged.loc[tmd_mask, 'residue_position'])

            df_tmd_AM = df_protein_AM[df_protein_AM['residue_position'].isin(tmd_positions)]
            df_tmd_ESM = df_protein_ESM[df_protein_ESM['residue_position'].isin(tmd_positions)]

            N = 10
            top_tmd_AM = df_tmd_AM.nlargest(N, 'am_rank_score')
            top_tmd_ESM = df_tmd_ESM.nlargest(N, 'esm_rank_score')

            top_annotation_am_tmd = annotate_top_variants(seq_length_dataset, top_tmd_AM)
            annot_dict_am_tmd = {i + 1: val for i, val in enumerate(top_annotation_am_tmd)}
            merged['AM_TMD_top10'] = merged['residue_position'].map(annot_dict_am_tmd).fillna('')

            top_annotation_esm_tmd = annotate_top_variants(seq_length_dataset, top_tmd_ESM)
            annot_dict_esm_tmd = {i + 1: val for i, val in enumerate(top_annotation_esm_tmd)}
            merged['ESM_TMD_top10'] = merged['residue_position'].map(annot_dict_esm_tmd).fillna('')

            if topology == "singlespan":
                jmd_mask = merged['Juxtamembrane'].astype(bool)
                tmdjmd_mask = tmd_mask | jmd_mask
                tmdjmd_positions = set(merged.loc[tmdjmd_mask, 'residue_position'])

                df_tmdjmd_AM = df_protein_AM[df_protein_AM['residue_position'].isin(tmdjmd_positions)]
                df_tmdjmd_ESM = df_protein_ESM[df_protein_ESM['residue_position'].isin(tmdjmd_positions)]

                top_tmdjmd_AM = df_tmdjmd_AM.nlargest(N, 'am_rank_score')
                top_tmdjmd_ESM = df_tmdjmd_ESM.nlargest(N, 'esm_rank_score')

                top_annotation_am_tmdjmd = annotate_top_variants(seq_length_dataset, top_tmdjmd_AM)
                annot_dict_am_tmdjmd = {i + 1: val for i, val in enumerate(top_annotation_am_tmdjmd)}
                merged['AM_TMDJMD_top10'] = merged['residue_position'].map(annot_dict_am_tmdjmd).fillna('')

                top_annotation_esm_tmdjmd = annotate_top_variants(seq_length_dataset, top_tmdjmd_ESM)
                annot_dict_esm_tmdjmd = {i + 1: val for i, val in enumerate(top_annotation_esm_tmdjmd)}
                merged['ESM_TMDJMD_top10'] = merged['residue_position'].map(annot_dict_esm_tmdjmd).fillna('')

            elif topology == "multispan":
                merged.drop(columns=['Juxtamembrane', 'AM_TMDJMD_top10', 'ESM_TMDJMD_top10',
                                     'AM_TMD_top10', 'ESM_TMD_top10'], inplace=True, errors='ignore')

                merged['Transmembrane'] = merged['Transmembrane'].replace(
                    to_replace=r'^Helical(;.*)?$', value='Helical', regex=True)

                helix_id = 1
                new_labels = []
                in_helix = False
                for val in merged['Transmembrane']:
                    if val == 'Helical':
                        if not in_helix:
                            label = f'Helical_{helix_id}'
                            helix_id += 1
                            in_helix = True
                        new_labels.append(label)
                    else:
                        in_helix = False
                        new_labels.append(val)
                merged['Transmembrane'] = new_labels

                am_annot = [''] * seq_length_dataset
                esm_annot = [''] * seq_length_dataset

                helical_labels = sorted(set(lab for lab in new_labels if lab and lab.startswith("Helical_")))

                for helix_label in helical_labels:
                    helix_mask = merged['Transmembrane'] == helix_label
                    helix_positions = set(merged.loc[helix_mask, 'residue_position'])

                    df_helix_AM = df_protein_AM[df_protein_AM['residue_position'].isin(helix_positions)]
                    df_helix_ESM = df_protein_ESM[df_protein_ESM['residue_position'].isin(helix_positions)]

                    top_helix_AM = df_helix_AM.nlargest(N, 'am_rank_score')
                    top_helix_ESM = df_helix_ESM.nlargest(N, 'esm_rank_score')

                    top_am = annotate_top_variants(seq_length_dataset, top_helix_AM)
                    top_esm = annotate_top_variants(seq_length_dataset, top_helix_ESM)

                    for i in range(seq_length_dataset):
                        if top_am[i]:
                            am_annot[i] += top_am[i] + ','
                        if top_esm[i]:
                            esm_annot[i] += top_esm[i] + ','

                am_annot = [s.rstrip(',') for s in am_annot]
                esm_annot = [s.rstrip(',') for s in esm_annot]

                merged['AM_Top10_Helix'] = merged['residue_position'].map(lambda x: am_annot[x - 1] if 0 < x <= seq_length_dataset else '')
                merged['ESM_Top10_Helix'] = merged['residue_position'].map(lambda x: esm_annot[x - 1] if 0 < x <= seq_length_dataset else '')


            # Adjust order of columns in CSV files
            cols = [
                'residue_position', 'residue', 'AM_mean', 'ESM_mean',
                'AM_median', 'ESM_median', 
                'AM_std', 'ESM_std', 
                'AM_max', 'ESM_max',
                'AM_min', 'ESM_min',
                'AM_top_var', 'ESM_top_var',
                '2D_Structures', 'rASA',
                'Topological domain', 'Transmembrane', 'Juxtamembrane',
                'Domain', 'Signal', 'Region'
            ]

            # Add the last two columns only if topology is not None
            if topology == "singlespan":
                cols.extend(['AM_TMD_top10', 'ESM_TMD_top10'])
            elif topology is not None:
                cols.extend(['AM_Top10_Helix', 'ESM_Top10_Helix'])

            # Reindex the DataFrame to ensure the columns are in the desired order
            merged = merged.reindex(columns=cols)  
                   
            # Save results
            output_path = os.path.join(store_pathway, f"{protein}_statistics.csv")
            merged.to_csv(output_path, index=False)

        except Exception as e:
             # Get the traceback details
            tb = traceback.format_exc()
            print(f"[ERROR] Failed processing {protein}: {e}")
            print(f"Traceback details:\n{tb}")
            continue

    skipped_proteins_df = pd.DataFrame(skipped_proteins, columns=["protein_id", "seq_length_AM", "seq_length_UniProt"])
    skipped_proteins_df.to_csv(misaligned, index=False)