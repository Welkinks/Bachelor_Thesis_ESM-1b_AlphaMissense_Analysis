from pathlib import Path
from dotenv import load_dotenv
import os
import sys

# === Load environment variables ===
# Looks for .env in the project root (assumed to be 2 levels up)
dotenv_path = Path(__file__).resolve().parents[1] / ".env"
load_dotenv(dotenv_path=dotenv_path)

# === Get project root from .env OR fallback ===
project_root_env = os.getenv("PROJECT_ROOT")

if project_root_env:
    PROJECT_ROOT = Path(project_root_env)
else:
    # Fallback to this file's grandparent 
    PROJECT_ROOT = Path(__file__).resolve().parents[1]
    print(f"⚠️  PROJECT_ROOT not found in .env — using fallback: {PROJECT_ROOT}")

# === Add project root to Python path ===
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

# === Define standard project paths ===
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
IMAGES_DIR = RESULTS_DIR / "images"

ESM_PATH = PROCESSED_DIR / "ESM1b_rank_csv"
ESM_PATH_CSV = "_LLR_rank.csv"

AM_PATH = PROCESSED_DIR / "AlphaMissense_rank_csv"
AM_PATH_CSV = "_rank.csv"

DIFF_PATH = PROCESSED_DIR / "difference_rank_csv"
DIFF_PATH_CSV = "_rank_difference.csv"

N_OUT_RANK_PATH = PROCESSED_DIR / "N_Out_Statistics_Rank"
N_OUT_RANK_PATH_CSV = "_statistics.csv"

MULTISPAN_PATH = PROCESSED_DIR / "Multispan_Statistics_Rank"
MULTISPAN_PATH_CSV = "_statistics.csv"

N_OUT_PROTEIN_ID = PROCESSED_DIR / "Human_N_Out_Proteome.csv"
PROTEIN_IDS_CSV = PROCESSED_DIR / "Protein_IDs_Per_Experiment" / "intersection_protein_ids_to_be_ranked.csv"
DSSP_PROTEINS_PATH = PROCESSED_DIR / "5.2.Protein_Statistics" / "DSSP_5800"


# === Expose paths via function ===
def get_paths():
    return {
        "project_root": PROJECT_ROOT,
        "data": DATA_DIR,
        "processed": PROCESSED_DIR,
        "esm_path": ESM_PATH,
        "am_path": AM_PATH,
        "difference_path": DIFF_PATH,
        "results_path": RESULTS_DIR,
        "images_path": IMAGES_DIR,
        "n_out_rank": N_OUT_RANK_PATH,
        "multispan_rank": MULTISPAN_PATH,
        "esm_csv_suffix": ESM_PATH_CSV,
        "am_csv_suffix": AM_PATH_CSV,
        "diff_csv_suffix": DIFF_PATH_CSV,
        "n_out_rank_csv_suffix": N_OUT_RANK_PATH_CSV,
        "multispan_rank_csv_suffix": MULTISPAN_PATH_CSV,
        "protein_ids_intersection": PROTEIN_IDS_CSV,
        "n_out_protein_id": N_OUT_PROTEIN_ID,
        "dssp_proteins_path": DSSP_PROTEINS_PATH,
    }


# === Get paths for a specific protein ID ===
def get_paths_protein(protein_id="P05067"):
    """
    Returns a dictionary of relevant project paths for a specific protein ID.
    """
    paths = get_paths()
    return {
        "esm_path": paths["esm_path"] / f"{protein_id}{paths['esm_csv_suffix']}",
        "am_path": paths["am_path"] / f"{protein_id}{paths['am_csv_suffix']}",
        "difference_path": paths["difference_path"] / f"{protein_id}{paths['diff_csv_suffix']}",
        "n_out_rank_path": paths["n_out_rank"] / f"{protein_id}{paths['n_out_rank_csv_suffix']}",
        "multispan_rank_path": paths["multispan_rank"] / f"{protein_id}{paths['multispan_rank_csv_suffix']}",
        "protein_statistics_path": paths["processed"] / "Protein_Statistics" / f"{protein_id}_statistics.csv",
        "dssp_protein_path": paths["dssp_proteins_path"] / f"{protein_id}_statistics.csv",
    }



def get_aa_list():
    return [
        "A", "V", "L", "I", "M", "F", "W",  # Hydrophobic
        "S", "T", "N", "Q", "Y", "C",      # Polar uncharged
        "K", "R", "H",                     # Positively charged
        "D", "E",                          # Negatively charged
        "G", "P"                           # Special
    ]



# --- Define Colors for the entire Project --- #
# Models
COLORS_MODELS = {
    'AM': '#ab294c',
    'ESM': '#29ab88',
}
# Secondary Colors
COLORS_SECONDARY = {
    'BLUE': '#c7eaf7',
    'PURPLE': '#d5c7f7',
    'PINK': '#f7c7e4',
}
# CCP Colors
COLORS_CPP = {
    'TRANS': '#00F89D',
    'JUXTA': '#0000FF',
}