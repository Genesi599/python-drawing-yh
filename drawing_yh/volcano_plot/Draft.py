# batch_volcano.py
from pathlib import Path
from math import log10
from volcano import create_volcano_plot   # 确保 volcano.py 在同一目录或在 PYTHONPATH

# ========== 唯一需要手动改的地方 ==========
root_dir = Path(r"D:\Projects\Neutrophil_Aging\NET_bulk\Ctrl_vs_NETs_vs_NETs_DNase_I")   # 根目录
# =========================================

def main():
    # 递归找所有 *all.csv
    csv_files = root_dir.rglob('*all.csv')
    for csv_path in csv_files:
        print(f'Processing: {csv_path}')
        out_name = csv_path.with_suffix('').name + '_volcano.png'  # xxx_volcano.png
        out_path = csv_path.parent / out_name                     # 同目录

        create_volcano_plot(
            input_file=str(csv_path),
            output_file=str(out_path),
            x_threshold=0.5,
            y_threshold=-log10(0.05),
            lfc_col='log2FoldChange',
            p_col='padj',
            id_col='GeneName'
        )
        print(f'  -> saved: {out_path}\n')

if __name__ == '__main__':
    main()