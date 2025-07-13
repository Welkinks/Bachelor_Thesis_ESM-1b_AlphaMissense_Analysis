# Utility functions for loading data, saving results, classes,
                    # helper functions, common imports, etc.
# add __init__.py to the utils folder to make it a package
# In other python files use: from src import utils

import pandas as pd
import os
from pathlib import Path
from dotenv import load_dotenv

# Function to save a DataFrame to a CSV file without overwriting existing files
def save_csv_no_overwrite(df, filepath):
    base, ext = os.path.splitext(filepath)  
    counter = 1
    new_filepath = filepath
    while os.path.exists(new_filepath):
        new_filepath = f"{base}_{counter}{ext}"
        counter += 1
    df.to_csv(new_filepath, index=False)
    print(f"Saved to {new_filepath}")

# Function to compute highlight regions for a given protein
def compute_highlight_regions(model, protein_id, pathway, df_predicted_cleaned, TMD_region, rolling_average):

    # model: 'AlphaMissense', 'ESM1b'
    # protein_id = 'P05067'
    # mean_pathway = '/Users/doma/Documents/Bachelor_Arbeit/Code/data/processed'
    # df_predicted_cleaned 


    # Choose correct pathway to different models
    if model == 'AlphaMissense':
        mean_pathway = os.path.join(pathway, 'AlphaMissense_mean_csv')
    elif model == 'ESM1b':
        mean_pathway = os.path.join(pathway,'ESM_mean_csv')
    
    # Read the protein-specific CSV file
    file_suffixes = {
        'AlphaMissense': '_avg.csv',
        'ESM1b': '_mean_LLR.csv'}

    if model in file_suffixes:
        entry_protein_pathway = os.path.join(mean_pathway, protein_id + file_suffixes[model])
        
        if not os.path.exists(entry_protein_pathway):
            print(f"File {entry_protein_pathway} does not exist.")
            return {
            "highlight_start": None,
            "highlight_end": None,
            "protein_data": None,
            "average_cleavage_region": None,
            "TMD_start": None,
            "TMD_end": None
            }
        
        # Read the CSV file into a DataFrame
        protein_data = pd.read_csv(entry_protein_pathway)
    else:
        raise ValueError("Model must be either 'AlphaMissense' or 'ESM1b'")
    

    # Get the row corresponding to the protein from df_predicted_cleaned
    row = df_predicted_cleaned.loc[df_predicted_cleaned['entry'] == protein_id].squeeze()

    # Compute TMD region
    TMD_start = row['len_signal_pep'] + row['len_top_n'] + 1
    TMD_end = TMD_start + row['len_tmd'] - 1


    if TMD_region == '11AA':
        # Compute tmd cleavage region = 24. to 34. AA in JMD-TMD-JMD-c region
        highlight_start = TMD_end - 6
        highlight_end = TMD_end + 4

    elif TMD_region == '40AA':
        # Compute JMD-N TMD JMD-C region
        highlight_start = TMD_start - 10
        highlight_end = TMD_end + 10

    else:
        raise ValueError("TMD_region must be either '11AA' or '40AA'")
    

    # Compute rolling average for plotting and unify column names
    if model == 'AlphaMissense':
        protein_data['smoothed_pathogenicity'] = protein_data['avg_pathogenicity'].rolling(window=rolling_average, center=True).mean()
        protein_data = protein_data.rename(columns={'avg_pathogenicity': 'residue_pathogenicity'})

    elif model == 'ESM1b':
        protein_data['smoothed_pathogenicity'] = protein_data['LLR_score'].rolling(window=rolling_average, center=True).mean()
        protein_data = protein_data.rename(columns={'position': 'residue_position', 'LLR_score': 'residue_pathogenicity', 'amino_acid': 'residue'})
        



    # Boolean mask for highlight range - only for plotting
    #if model == 'AlphaMissense':
        # For AlphaMissense, the highlight region is based on the rolling average
    #    highlight_mask = (protein_data['residue_position'] >= highlight_start) & (protein_data['residue_position'] <= highlight_end)
    #elif model == 'ESM1b':
    #    highlight_mask = (protein_data['position'] >= highlight_start) & (protein_data['position'] <= highlight_end)
    #else:
    #    raise ValueError("Model must be either 'AlphaMissense' or 'ESM1b'")
    


    # Set "residue_position" as the index for indexing purposes
    protein_data = protein_data.set_index('residue_position')
    
    # Compute average pathogenicity score of the highlighted region
    average_cleavage_region = protein_data['residue_pathogenicity'].loc[highlight_start:highlight_end].mean()
    
    # Reset the index back to default if needed later
    protein_data = protein_data.reset_index()

    # Convert all values to integers except for "average_cleavage_region"
    result = {
        "highlight_start": int(highlight_start),
        "highlight_end": int(highlight_end),
        #"highlight_mask": highlight_mask,
        "protein_data": protein_data,
        "average_cleavage_region": average_cleavage_region,  # Keep original datatype
        "TMD_start": int(TMD_start),
        "TMD_end": int(TMD_end)
    }

    # Return the result dictionary
    return result





def group_by_pred_class(df):
    """
    This function takes a dataframe and groups it by the 'pred_class' column.
    It returns a dictionary where each key is a unique value from 'pred_class'
    and the value is the corresponding dataframe.
    """
    # Group the dataframe by "pred_class" 
    grouped_df = df.groupby("pred_class")

    # Create a dictionary of dataframes for each pred_class category
    grouped_dfs = {group: data for group, data in grouped_df}

    return grouped_dfs






 # Load environment variables from .env (in current working dir or parent)
load_dotenv(dotenv_path=Path(__file__).resolve().parent.parent / ".env")

def get_project_paths():
    """
    Returns a dictionary of relevant project paths relative to the notebook location.
    """
    
    # Get the project root directory
    project_root = Path(os.getenv("PROJECT_ROOT"))
    
    paths = {
        "project_root": project_root,
        "esm_path": project_root / "data" / "processed" / "ESM1b_rank_csv",
        "am_path": project_root / "data" / "processed" / "AlphaMissense_rank_csv",
        "difference_path": project_root / "data" / "processed" / "difference_rank_csv",
        "images_path": project_root / "results" / "images"
    }

    return paths
