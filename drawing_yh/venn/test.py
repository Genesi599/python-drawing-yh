# -*- coding: utf-8 -*-
"""
合并 monkey*up.csv 并按 gene 去重
"""
import pandas as pd
from pathlib import Path

folder = Path(r"D:\Projects\Thymus_Aging\overlap\thymus_data")
out_file = folder / "merged_monkey_down.csv"

# 找到所有 monkey*up.csv
files = list(folder.glob("monkey*down.csv"))
if not files:
    raise FileNotFoundError("目录下没有 monkey*down.csv")

# 逐个读入并拼接
dfs = [pd.read_csv(f) for f in files]
merged = pd.concat(dfs, ignore_index=True)

# 按 gene 去重（保留第一次出现）
if 'gene' not in merged.columns:
    # 如果列名不是 gene，尝试第 1 列
    merged.rename(columns={merged.columns[0]: 'gene'}, inplace=True)
dedup = merged.drop_duplicates(subset='gene')

# 保存
dedup.to_csv(out_file, index=False)
print(f"✅ 已合并 {len(files)} 个文件，共 {len(dedup)} 条基因\n-> {out_file}")