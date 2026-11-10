"""Input and output locations for every script in this repository.

Defaults (relative to the repository root):
  data/Data_cleaned.csv       the 310-study dataset (main article figures and values)
  data/abstract_dataset.csv   the global abstract dataset (Supplementary File figures)
  output/                     generated figures and tables

Each can be overridden with an environment variable:
  META_ANALYSIS_DATA_FILE, META_ANALYSIS_GLOBAL_DATA_FILE, META_ANALYSIS_OUTPUT_DIR
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_FILE = Path(os.environ.get("META_ANALYSIS_DATA_FILE", ROOT / "data" / "Data_cleaned.csv"))
GLOBAL_DATA_FILE = Path(os.environ.get("META_ANALYSIS_GLOBAL_DATA_FILE", ROOT / "data" / "abstract_dataset.csv"))
OUTPUT_DIR = Path(os.environ.get("META_ANALYSIS_OUTPUT_DIR", ROOT / "output"))
