# Bioinformatics Thesis Project

This repository contains the Python code, data processing steps, and analysis notebooks related to my bachelor thesis in bioinformatics. The project focuses on comparing and visualizing pathogenicity predictions from **AlphaMissense** and **ESM-1b**, and integrating these predictions with structural protein features.

---

## Installation

```bash
git clone <https://github.com/Welkinks/Bachelor_Thesis_PLM_VEP_Analysis.git>
cd <Bachelor_Thesis_PLM_VEP_Analysis>

# Create and activate a virtual environment (recommended)
python -m venv venv
source venv/bin/activate   # or venv\Scripts\activate on Windows

# Install dependencies
pip install -r requirements.txt
```

---

## Project Setup

To use the notebooks, always run the following snippet at the beginning to set the project root:

```python
# --- Project Setup ---
from setup_notebook import setup_project_root
setup_project_root()
```

---

## Raw Data

Large raw data files required for the analysis are **not included** in this repository due to size and version control limitations. Please download them manually:

1. **AlphaMissense**  
   - File: `AlphaMissense_aa_substitutions.tsv.gz`  
   - Source: [Zenodo Repository](https://zenodo.org/records/8208688)

2. **ALL_hum_isoforms_ESM1b_LLR**  
   - File: `ALL_hum_isoforms_ESM1b_LLR.zip`  
   - Source: [Hugging Face - ntranoslab](https://huggingface.co/spaces/ntranoslab/esm_variants/tree/main)

After downloading:
- Unzip files (if needed).
- Place the extracted files/folders into: `data/raw/`

---

## Data Preprocessing

1. **AlphaMissense TSV to CSV**  
   Run `scripts/1.AM_TSV_to_CSV.py` to split the TSV file into per-protein CSV files.

2. **ESM-1b**  
   Already in appropriate per-protein CSV format.

3. **Rank Transformation**  
   - Run `notebooks/0.1.Rank_scores.ipynb` to rank-transform AlphaMissense and ESM-1b scores.
   - Output stored in:  
     - `data/processed/AlphaMissense_rank_csv`  
     - `data/processed/ESM1b_rank_csv`
   - Visualizations: `notebooks/0.2.Ranking_dotplot.ipynb`
   - Precomputed datasets (~1.5 GB x 2 = 3 GB) available on request.

4. **DSSP Annotations**  
   - Run `notebooks/0.3.Per_Residue_Annotations_DSSP.ipynb` to add rASA and secondary structure info.  
   - Runtime: ~12 hours for full human proteome.  
   - Precomputed dataset (~1.5 GB) available on request.

---

## Experiments and Visualizations

1. **TMD/JMD Experiment**  
   Compare rank scores between Gamma-secretase substrates and non-substrates.  
   File: `notebooks/1.1.TMD_JMD_Experiment.ipynb`

2. **Visualization Function Example**  
   Plot per-residue pathogenicity scores with optional rASA and secondary structure overlays.  
   File: `notebooks/1.2.Pathogenicity_Plots.ipynb`

3. **rASA Correlation Analysis**  
   Correlate pathogenicity scores with rASA, and compute per-variant correlations.  
   File: `notebooks/1.3.rASA_Corr_Pathog.ipynb`

4. **Dendrogram Heatmap**  
   Visualize similarity trends and hydrophobicity clusters.  
   Files:  
   - `notebooks/1.4.1.Dendrogram_Heatmap.ipynb`  
   - `notebooks/1.4.2.Heatmaps_bias_computation.ipynb`

5. **Clustering by Secondary Structure**  
   Explore average pathogenicity scores per secondary structure.  
   File: `notebooks/1.5.Clustering_2D.ipynb`

6. **Gene-Level Pathogenicity**  
   Compute mean pathogenicity per protein to assess evolutionary relevance.  
   File: `notebooks/1.6.Gene_Level_agg.ipynb`

---

## Acknowledgments

This project makes use of the following resources:

- **AlphaMissense**: Cheng J. et al. (2023), [Zenodo Repository](https://zenodo.org/records/8208688), Paper: (https://www.science.org/doi/10.1126/science.adg7492)
  
- **ESM-1b**: Rives A. et al. (2021), [Hugging Face - ntranoslab](https://huggingface.co/spaces/ntranoslab/esm_variants/tree/main), Paper: (https://pmc.ncbi.nlm.nih.gov/articles/PMC10484790/)
  
- **AlphaFold2** Protein Structure Database: Varadi M. et al. (2022), "AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models", Nucleic Acids Research, 50(D1): D439–D444 (https://alphafold.ebi.ac.uk/download)
  
- **DSSP**: The newest paper and algorithm used: (https://onlinelibrary.wiley.com/doi/10.1002/pro.70208), the original paper for DSSP: (https://pubmed.ncbi.nlm.nih.gov/6667333/)
