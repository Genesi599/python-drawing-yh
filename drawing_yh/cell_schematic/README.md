# Flat eukaryotic cell schematic

Reusable editable template for subcellular-localization figures. The template
contains a phospholipid bilayer, nucleus connected to rough ER, Golgi,
mitochondria, vesicles, secreted proteins, ECM, a membrane projection, and an
apposed-cell junction. Additional organelle SVGs are inserted from a placement
list, so their number and positions are controlled in code.

```python
from drawing_yh.cell_schematic import copy_flat_cell_template

copy_flat_cell_template("figure/cell_localization.svg")
```

Customize quantity, position, size, or rotation:

```python
from drawing_yh.cell_schematic import IconPlacement, copy_flat_cell_template

placements = [
    IconPlacement("lysosome_1", (430, 430), 72, asset="lysosome"),
    IconPlacement("lysosome_2", (520, 470), 64, asset="lysosome"),
    IconPlacement("peroxisome", (1100, 600), 64),
]
copy_flat_cell_template("figure/cell_localization.svg", placements=placements)
```

Coordinates use the bundled 1400 x 1182 reference preview. Use
`default_anchor_map()` for the validated localization and branch endpoints.
The output includes a sibling asset directory and an attribution file; keep
them beside the SVG when moving the figure.
