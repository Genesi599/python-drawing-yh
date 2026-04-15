"""
drawing_yh — Python 科研作图工具包
"""
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("drawing-yh")
except PackageNotFoundError:
    __version__ = "0.0.0"
