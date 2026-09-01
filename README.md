# python-drawing-yh (`drawing_yh`)

Publication-quality scientific plotting toolkit for Python — the plotting layer behind multi-module slide-style research reports.

`import drawing_yh` applies the house style automatically (Arial, 8 pt base, editable text in SVG/PDF, three-format export), and every chart family ships as a reusable template instead of ad-hoc scripts.

## Highlights

- **Chart families** — bar / box / heatmap / scatter / violin / volcano / PCA / venn / pie / network / chord / dotplot / dumbbell / lollipop / embedding / dose-response / Michaelis-Menten … each with a documented template module.
- **Report suite (`drawing_yh.report`)** — Markdown → multi-module slide-style HTML websites (three-level nav: module / chapter / section) and same-source native PPTX inheriting a master theme. See the companion skill [local-report-builder](https://github.com/Genesi599/local-report-builder) for the full playbook.
- **Text audit helpers** — `get_text_bboxes` / `find_overlaps` to programmatically verify no text overlaps or canvas overflow (the two hard rules every figure must pass).
- **Style constants** — Okabe-Ito / DEEP-20 / age / sex / species palettes, `compute_figsize` layout helpers, `save_fig` one-call pdf+png+svg export.
- **Example gallery** — `example_gallery/` renders one preview per template; rebuilt automatically by a pre-commit hook when templates change.

## Install

```bash
pip install -e .          # from a checkout
```

Requires Python ≥ 3.11, matplotlib / numpy / pandas / scipy, plus `matplotlib-venn` and `pycirclize` for those chart families. The report suite additionally needs system `pandoc`, `python-pptx` and `Pillow`.

## Usage sketch

```python
import drawing_yh                    # applies rcParams (Arial 8pt, editable fonts)
from drawing_yh import save_fig      # pdf + png + svg in one call

fig, ax = plt.subplots()
...
save_fig(fig, "my_figure.pdf", also=(".png", ".svg"))
```

## Conventions (STANDARDS.md)

- Base font 8 pt locked; figure size adapts to content, never the other way.
- All text in English inside figures; significance markers as plain black `*`.
- Every figure ships with a plotted-data CSV sidecar where applicable.

See `STANDARDS.md` for the full checklist.

## License

MIT
