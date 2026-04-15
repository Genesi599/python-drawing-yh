import pandas as pd
from pathlib import Path
import scanpy as sc
path = Path(r"D:\leukocyte_single_cell\Monkey\monkey-B_cell\subluster\data")
name = "after_Annotation.h5ad"
adata = sc.read_h5ad(path / name)

# 1. 查看adata中有多少基因
print(f"总基因数: {len(adata.var_names)}")

# 2. 打印前20个基因名，检查格式
print("前20个基因名:")
print(adata.var_names[:20].tolist())

# 3. 搜索您要找的基因
search_genes = ['HOPX', 'LITAF', 'PLEK', 'ZBTB32']
for gene in search_genes:
    if gene in adata.var_names:
        print(f"✓ 找到: {gene}")
    else:
        print(f"✗ 未找到: {gene}")

# 4. 检查是否是大小写问题
print("\n检查大小写问题:")
gene_names_lower = adata.var_names.str.lower()
for gene in search_genes:
    matches = adata.var_names[gene_names_lower == gene.lower()].tolist()
    if matches:
        print(f"  {gene} -> 实际名称: {matches}")
    else:
        print(f"  {gene} -> 未找到匹配")