"""Overlay heatmap templates."""
from .heatmap_fan_species import plot_heatmap_fan_species
from .heatmap_with_icon import (
    load_species_icons,
    merge_species_data,
    normalize_species_name,
    plot_heatmap_with_icons,
)

__all__ = [
    "plot_heatmap_fan_species",
    "load_species_icons",
    "merge_species_data",
    "normalize_species_name",
    "plot_heatmap_with_icons",
]
