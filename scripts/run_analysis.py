# Main script to run analysis using utility functions

from src import utils

data = utils.load_data('data/processed/sample_data.csv')
print(data.head())
