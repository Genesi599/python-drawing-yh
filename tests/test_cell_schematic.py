from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from drawing_yh.cell_schematic import (
    DEFAULT_ICON_PLACEMENTS,
    copy_flat_cell_template,
    default_anchor_map,
    flat_cell_preview_path,
    flat_cell_template_path,
)


SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"


def test_flat_cell_template_resources_exist() -> None:
    assert flat_cell_template_path().is_file()
    assert flat_cell_preview_path().is_file()
    anchors = default_anchor_map()
    assert anchors["reference_size"] == [1400, 1182]
    assert anchors["summary_anchors_px"]["Secreted / extracellular"] == [92, 235]


def test_copy_flat_cell_template_is_portable(tmp_path: Path) -> None:
    output = copy_flat_cell_template(tmp_path / "cell.svg")
    assert output.is_file()
    assert (tmp_path / "cell_ATTRIBUTION.md").is_file()

    root = ET.parse(output).getroot()
    detail_group = root.find(f".//{{{SVG_NS}}}g[@id='detail_organelle_icons']")
    assert detail_group is not None
    detail_images = list(detail_group.findall(f"{{{SVG_NS}}}image"))
    assert len(detail_images) == len(DEFAULT_ICON_PLACEMENTS)

    for image in root.iter(f"{{{SVG_NS}}}image"):
        href = image.get("href") or image.get(f"{{{XLINK_NS}}}href")
        assert href
        assert (output.parent / href).is_file(), href
