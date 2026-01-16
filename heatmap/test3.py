import pandas as pd
import sys

EXPR  = r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\abundance_sample_x_protein.csv"
PAT   = r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\Pattern_Analysis\all_proteins_pattern.csv"

e = pd.read_csv(EXPR, index_col=0)
p = pd.read_csv(PAT,  index_col=0)

print("=== 表达矩阵（列名即蛋白名）前 5 个示例 ===")
print(e.columns[:5].tolist())
print("\n=== pattern 文件 index（蛋白名）前 5 个示例 ===")
print(p.index[:5].tolist())

overlap = e.columns.intersection(p.index)
print(f"\n=== 两者交集数量：{len(overlap)} ===")
if len(overlap) < 10:
    print("交集太少，请检查命名规则！")
else:
    print("交集看起来正常，继续下一步。")