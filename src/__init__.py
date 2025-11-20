"""
Host Filtration Package
ML-based host genome filtering for metagenomic sequencing data.
"""

__version__ = "0.1.0"
__author__ = "Tristan Ferdinand"
__email__ = "tferdinand@ucsd.edu"

from .hostfiltration import cli
from .hostfiltration import data
from .hostfiltration import models

__all__ = ["cli", "data", "models"]