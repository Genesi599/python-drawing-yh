"""Reusable flat eukaryotic-cell schematic template."""

from __future__ import annotations

import json
import shutil
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple, Union


SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
REFERENCE_SIZE = (1400.0, 1182.0)
VIEWBOX = (180.0, 45.0, 1220.0, 1030.0)

PACKAGE_DIR = Path(__file__).resolve().parent
TEMPLATE_DIR = PACKAGE_DIR / "templates"
BASE_TEMPLATE = TEMPLATE_DIR / "flat_eukaryotic_cell_base.svg"
PREVIEW = TEMPLATE_DIR / "flat_eukaryotic_cell_preview.png"
ANCHORS = TEMPLATE_DIR / "flat_eukaryotic_cell_anchors.json"
ATTRIBUTION = TEMPLATE_DIR / "ATTRIBUTION.md"
ASSET_DIR = TEMPLATE_DIR / "assets"


@dataclass(frozen=True)
class IconPlacement:
    """One editable SVG icon placed in reference-image pixel coordinates."""

    name: str
    center_px: Tuple[float, float]
    width_px: float
    rotation_deg: float = 0.0
    asset: Optional[str] = None


DEFAULT_ICON_PLACEMENTS = (
    IconPlacement("collagen", (700, 65), 300),
    IconPlacement("endosome", (340, 430), 76),
    IconPlacement("lysosome", (430, 430), 72),
    IconPlacement("vacuole", (385, 520), 76),
    IconPlacement("peroxisome", (1100, 600), 64),
    IconPlacement("proteasome", (1015, 810), 42, asset="proteasome_simplified"),
    IconPlacement("centrosome", (1020, 650), 52, asset="centrosome_simplified"),
)


def flat_cell_template_path() -> Path:
    """Return the packaged editable base SVG path."""

    return BASE_TEMPLATE


def flat_cell_preview_path() -> Path:
    """Return the packaged PNG preview path."""

    return PREVIEW


def flat_cell_attribution_path() -> Path:
    """Return source and license notes for every bundled asset."""

    return ATTRIBUTION


def default_anchor_map() -> Dict[str, object]:
    """Load summary and branch anchors in 1400 x 1182 reference pixels."""

    return json.loads(ANCHORS.read_text(encoding="utf-8"))


def _asset_aspect_ratio(path: Path) -> float:
    root = ET.parse(path).getroot()
    viewbox = root.get("viewBox")
    if not viewbox:
        raise ValueError(f"SVG asset has no viewBox: {path}")
    _, _, width, height = (float(value) for value in viewbox.replace(",", " ").split())
    if width <= 0 or height <= 0:
        raise ValueError(f"Invalid SVG viewBox: {path}")
    return width / height


def _reference_to_svg(
    point: Tuple[float, float],
) -> Tuple[float, float]:
    x, y = point
    ref_width, ref_height = REFERENCE_SIZE
    vb_x, vb_y, vb_width, vb_height = VIEWBOX
    return (
        vb_x + x * vb_width / ref_width,
        vb_y + y * vb_height / ref_height,
    )


def copy_flat_cell_template(
    destination: Union[str, Path],
    *,
    placements: Iterable[IconPlacement] = DEFAULT_ICON_PLACEMENTS,
) -> Path:
    """Create a portable editable SVG plus its local asset directory.

    Add, remove, or replace ``IconPlacement`` entries to control organelle
    quantity, position, size, and rotation. Coordinates use the 1400 x 1182
    reference preview, matching ``default_anchor_map()``.
    """

    destination = Path(destination)
    if destination.suffix.lower() != ".svg":
        destination = destination.with_suffix(".svg")
    destination.parent.mkdir(parents=True, exist_ok=True)

    output_assets = destination.parent / f"{destination.stem}_assets"
    shutil.copytree(ASSET_DIR, output_assets, dirs_exist_ok=True)

    ET.register_namespace("", SVG_NS)
    ET.register_namespace("xlink", XLINK_NS)
    tree = ET.parse(BASE_TEMPLATE)
    root = tree.getroot()

    for image in root.iter(f"{{{SVG_NS}}}image"):
        href = image.get("href") or image.get(f"{{{XLINK_NS}}}href")
        if href and href.startswith("assets/"):
            rewritten = f"{output_assets.name}/{href[len('assets/'):]}"
            image.set("href", rewritten)
            image.set(f"{{{XLINK_NS}}}href", rewritten)

    detail_group = ET.SubElement(root, f"{{{SVG_NS}}}g", {"id": "detail_organelle_icons"})
    vb_scale_x = VIEWBOX[2] / REFERENCE_SIZE[0]
    for placement in placements:
        asset_name = placement.asset or placement.name
        source_asset = ASSET_DIR / "detail" / f"{asset_name}.svg"
        if not source_asset.exists():
            raise FileNotFoundError(f"Unknown cell-schematic asset: {source_asset}")

        aspect = _asset_aspect_ratio(source_asset)
        width = placement.width_px * vb_scale_x
        height = width / aspect
        center_x, center_y = _reference_to_svg(placement.center_px)
        href = f"{output_assets.name}/detail/{asset_name}.svg"
        attrs = {
            "id": f"detail_{placement.name}",
            "data-location": placement.name,
            "x": f"{center_x - width / 2:.3f}",
            "y": f"{center_y - height / 2:.3f}",
            "width": f"{width:.3f}",
            "height": f"{height:.3f}",
            "preserveAspectRatio": "xMidYMid meet",
            "href": href,
            f"{{{XLINK_NS}}}href": href,
        }
        if placement.rotation_deg:
            attrs["transform"] = (
                f"rotate({placement.rotation_deg:.3f} {center_x:.3f} {center_y:.3f})"
            )
        ET.SubElement(detail_group, f"{{{SVG_NS}}}image", attrs)

    tree.write(destination, encoding="utf-8", xml_declaration=True)
    shutil.copy2(ATTRIBUTION, destination.with_name(f"{destination.stem}_ATTRIBUTION.md"))
    return destination


__all__ = [
    "IconPlacement",
    "DEFAULT_ICON_PLACEMENTS",
    "REFERENCE_SIZE",
    "VIEWBOX",
    "copy_flat_cell_template",
    "default_anchor_map",
    "flat_cell_attribution_path",
    "flat_cell_preview_path",
    "flat_cell_template_path",
]
