"""
Example execution script for the automated DFT pipeline.
"""

import os
import sys

# 1. Change the Current Working Directory to the example folder
example_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(example_dir)

# 2. Add the root project directory to the Python path to import main.py
project_root = os.path.dirname(os.path.dirname(example_dir))
sys.path.append(project_root)

from main import run_pipeline

if __name__ == "__main__":
    # The pipeline will now generate inputs/, outputs/, tmp/ and search for pseudo/ locally
    run_pipeline("NiTe2")