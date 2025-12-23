import json, csv, pathlib

in_file  = pathlib.Path(r'C:\Users\yh109\OneDrive\桌面\讨论\gene_annotation_simple.json')
out_file = in_file.with_name('gene_annotation_simple_wide.csv')

data = json.loads(in_file.read_text(encoding='utf-8'))

# 先算最大文件数，决定列宽
max_files = max(len(info['files']) for info in data.values())
header    = ['gene'] + [f'{tag}{i}' for i in range(1, max_files+1) for tag in ('file','tissue','contrast','direction')] + ['geneSets']

rows = []
for gene, info in data.items():
    row = {'gene': gene, 'geneSets': '|'.join(info['geneSets'])}
    for idx, (tissue, contrast, direction) in enumerate(info['files'], 1):
        row[f'file{idx}']      = f'{tissue}|{contrast}|{direction}'
        row[f'tissue{idx}']    = tissue
        row[f'contrast{idx}']  = contrast
        row[f'direction{idx}'] = direction
    rows.append(row)

with out_file.open('w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=header)
    writer.writeheader()
    writer.writerows(rows)

print('✅ 宽表已生成：', out_file)