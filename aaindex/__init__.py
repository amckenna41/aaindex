from importlib.metadata import version, PackageNotFoundError

from .aaindex1 import *
from .aaindex2 import *
from .aaindex3 import *

# Single-source version from installed package metadata
try:
    __version__ = version("aaindex")
except PackageNotFoundError:
    __version__ = "0.0.0"

__author__ = "AJ McKenna: https://github.com/amckenna41"
__license__ = "MIT"

__all__ = ["AAIndex1", "aaindex1", "AAIndex2", "aaindex2", "AAIndex3", "aaindex3"]