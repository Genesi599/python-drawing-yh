import pandas as pd
from pathlib import Path

base_dir = Path(r"D:\Projects\B_Cell_Aging\B_cell_ratio")
csv_files = base_dir.glob("*.csv")

for file_path in csv_files:
    df = pd.read_csv(file_path)
    df['b_cell_ratio'] = (df['B cells'] + df['Plasma cells']) / df['total_cells']
    df.to_csv(file_path, index=False)
    print(file_path)