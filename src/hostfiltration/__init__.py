"""
Host Filtration Package
ML-based host genome filtering for metagenomic sequencing data.
"""

__version__ = "0.1.0"
__author__ = "Tristan Ferdinand"
__email__ = "tferdinand@ucsd.edu"

from . import cli
from . import data
from . import models

__all__ = ["cli", "data", "models"]