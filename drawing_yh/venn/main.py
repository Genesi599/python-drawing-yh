import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib_venn import venn2, venn3
import numpy as np
from pathlib import Path
from matplotlib.patches import FancyArrowPatch


def plot_venn2_or_3(gene_sets: dict, labels: list, colors: list, out_png: str, title: str):
    n = len(labels)
    assert n in {2, 3}

    sets = [gene_sets[lbl] for lbl in labels]
    inter = set.intersection(*sets)
    if not inter:
        print("❌ 交集为空，跳过")
        return False

    fig = plt.figure(figsize=(12, 12))
    ax_main = fig.add_axes([0.1, 0.35, 0.8, 0.6])
    ax_text = fig.add_axes([0.1, 0.05, 0.8, 0.28])
    ax_text.axis('off')
    ax_main.set_facecolor('none')

    sizes, weighted = [len(s) for s in sets], [np.log10(len(s) + 1) for s in sets]
    alpha = 0.3

    if n == 2:
        v = venn2(subsets=(weighted[0], weighted[1], min(weighted) * 0.5),
                 set_labels=labels, ax=ax_main, set_colors=colors, alpha=alpha)
        patch_ids = ['10', '01', '11']
    else:
        w0, w1, w2 = weighted
        overlap = min(weighted) * 0.3
        v = venn3(subsets=(w0, w1, overlap, w2, overlap, overlap, overlap),
                 set_labels=labels, ax=ax_main, set_colors=colors, alpha=alpha)
        patch_ids = ['100', '010', '001', '110', '101', '011', '111']

    for pid in patch_ids:
        p = v.get_patch_by_id(pid)
        if p:
            p.set_edgecolor('none')
            p.set_zorder(1)

    for t in v.set_labels:
        if t:
            t.set_fontsize(30)
            t.set_fontweight('bold')

    for pid in patch_ids:
        lbl = v.get_label_by_id(pid)
        if lbl:
            real_sz = len(set.intersection(*[sets[i] for i, c in enumerate(pid) if c == '1']))
            lbl.set_text(str(real_sz))
            lbl.set_fontsize(25)
            lbl.set_fontweight('bold')

    center_patch = v.get_patch_by_id('11' if n == 2 else '111')
    if center_patch:
        verts = center_patch.get_path().vertices
        cx, cy = np.mean(verts[:, 0]), np.mean(verts[:, 1])
        y_min, y_max = ax_main.get_ylim()
        arrow = FancyArrowPatch((cx, cy - 0.05), (0, y_min - (y_max - y_min) * 0.05),
                               arrowstyle='->', mutation_scale=30, linewidth=4,
                               color='grey', connectionstyle="arc3,rad=0.1", zorder=100, clip_on=False)
        ax_main.add_patch(arrow)

    gene_list = sorted(inter)
    lines = [",  ".join(gene_list[i:i + 3]) for i in range(0, len(gene_list), 3)]
    gene_txt = "\n".join(lines)
    inter_title = f"[{' ∩ '.join(labels)}]  ({len(gene_list)} genes)"

    title_text = ax_text.text(0.5, 0.85, inter_title, fontsize=30, fontweight='bold',
                             ha='center', va='top', color='#1A5490', transform=ax_text.transAxes)
    title_text.set_path_effects([path_effects.withStroke(linewidth=4, foreground='white', alpha=0.9)])
    genes_text = ax_text.text(0.5, 0.7, gene_txt, fontsize=30, fontweight='bold',
                             ha='center', va='top', family='monospace', transform=ax_text.transAxes)
    genes_text.set_path_effects([path_effects.withStroke(linewidth=4, foreground='white', alpha=0.9)])

    ax_main.set_title(title, fontsize=30, fontweight='bold', pad=15)
    ax_main.axis('off')
    plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✅ 已保存 -> {out_png}")
    return True
