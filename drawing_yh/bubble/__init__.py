"""drawing_yh.bubble — generic bubble (dot) plot for categorical axes.

Exports
-------
bubble_plot
    Draw a grid of bubbles: x/y = category axes, size = numeric, color = category
    or diverging direction.  Not bound to any specific dataset.

two_tone_ticklabels
    Replace axis tick labels with multi-coloured AnnotationBbox segments.
    Returns the artist list needed for ``bbox_extra_artists`` in ``save_fig``.
"""
from .bubble_plot import bubble_plot, two_tone_ticklabels

__all__ = ["bubble_plot", "two_tone_ticklabels"]
