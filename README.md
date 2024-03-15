Tool for visualization of molecular docking results
# Installation
run 
```
conda env create -f environment.yml
conda activate dock_visualizer
```
# Usage
```streamlit run visualizer.py -- --csv DATA_PATHS.csv```

Note that the double dash is required. 
```DATA_PATHS.csv``` is a csv file with columns [_complex_name, mol_pred, mol_cond, mol_true_]

## Example usage
```streamlit run visualizer.py -- --csv data/example.csv```

