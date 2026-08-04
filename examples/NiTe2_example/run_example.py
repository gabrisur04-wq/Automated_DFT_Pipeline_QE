"""
Example execution script for the automated DFT pipeline.
"""

import sys
import os
import shutil

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, "../../"))

sys.path.append(project_root)

from main import run_pipeline

if __name__ == "__main__":
    prefix = "NiTe2"

    source_json = os.path.join(current_dir, f"{prefix}.json")
    target_json = os.path.join(project_root, f"{prefix}.json")

    try:
        shutil.copy(source_json, target_json)
        print(f"Configuration {prefix}.json copied to project root.")
    except FileNotFoundError:
        sys.exit(f"Error: {prefix}.json not found in the example directory.")

    # Switch to project root to generate standard output directories (inputs/, outputs/, tmp/)
    os.chdir(project_root)

    run_pipeline(prefix)