import os
from read_file import *
from predicted_viability import predect_viability
from heatmap_with_colorbar import get_heapmap

cell_name = 'Karpas299'
drug_names = ['JQ-1', 'CC-90009']
date = '20250120'

merged_list = []
folder_path = "data/raw"
entries = os.listdir(folder_path)
for i in range(len(entries)):
    para_df = create_para_df(cell_name, drug_names, i, date)
    file_path = os.path.join(folder_path, entries[i])
    if os.path.isfile(file_path):
        merged_list += read_file(file_path, para_df, drug_names)

get_viablity(merged_list, cell_name, drug_names, date)

file_path = f'data/viability/{cell_name}-{drug_names[0]}+{drug_names[1]}-{date}.csv'
predect_viability(file_path, drug_names)

get_heapmap(drug_names, 'viability')