# NiTe2 Electronic Structure Example

This directory contains a complete configuration to calculate the electronic band structure of the Dirac semimetal NiTe2 using the Automated DFT Pipeline.

## Prerequisites

1. **Quantum ESPRESSO**: Ensure the executables (e.g., `pw.x`) are installed and available in your system's PATH.
2. **Pseudopotentials**: Ensure the required `.UPF` files for Nickel and Tellurium are located in the `pseudo/` directory.

## How to Run

You can launch the example directly using the provided script. From your terminal, navigate to this directory and run:

python run_example.py

## Expected Output

The script will dynamically create the following directories in the project root (if they do not already exist):
* `inputs/`: Contains the generated Quantum ESPRESSO input files.
* `tmp/`: Stores temporary files during the calculation.
* `outputs/`: Contains the final data, including the parsed band structure files and plots.

Once the pipeline finishes, check the `outputs/` folder for the graphical representation of the NiTe2 bands.