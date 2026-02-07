#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : test.py
@Date    : 2026/2/7 17:28
@Author  : yh109
"""
import os
import glob
import pandas as pd

data_dir = r"D:\Projects\B_Cell_Aging\B_cell_ratio"
csv_files = glob.glob(os.path.join(data_dir, "**/*.csv"), recursive=True)

for fpath in csv_files:
    df = pd.read_csv(fpath)
    if 'Tissue' in df.columns:
        df['Tissue'] = df['Tissue'].replace('AdiposeTissue', 'BAT')
        df.to_csv(fpath, index=False)
        print(f"Updated: {fpath}")

print(f"Total files processed: {len(csv_files)}")