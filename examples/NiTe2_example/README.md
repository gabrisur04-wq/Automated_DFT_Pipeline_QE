# NiTe2 Electronic Structure Example

This directory contains a complete configuration to calculate the electronic band structure of the Dirac semimetal NiTe2 using the Automated DFT Pipeline.

## Prerequisites

1. **Quantum ESPRESSO**: Ensure the executables (e.g., `pw.x`, `dos.x`, `bands.x`, `fs.x`) are installed and available in your system's PATH.
2. **Pseudopotentials**: Ensure the required `.UPF` files for Nickel and Tellurium are located in the `pseudo/` directory.

## How to Run

Ensure your terminal is located within this specific directory (`examples/NiTe2_example/`) before executing any commands.

### 1. Automated Execution
You can launch the default pipeline sequentially using the provided wrapper script (Scalar-Relativistic calculation):

```bash
python run_example.py
```

### 2. Modular Execution via CLI
For precise control over the pipeline, error recovery, or to execute specific physical models (such as Spin-Orbit coupling), you can invoke the main infrastructure engine directly using the `../../main.py` relative path.

*Note: The commands below use `--num_cores 4` and `--npool 1`, which is the recommended parallelization configuration for multi-core desktop processors to maximize stability and memory efficiency.*

*   **Phase 1: Structural Setup & Convergence**
    ```bash
    python ../../main.py --prefix NiTe2 setup --num_cores 4 --npool 1
    ```
*   **Phase 2: Scalar-Relativistic (SR) Calculation**
    ```bash
    python ../../main.py --prefix NiTe2 run --mode SR --num_cores 4 --npool 1
    ```
    *Optional:* You can restrict the execution to a specific step by appending `--step bands`, `--step dos`, or `--step fs`.
*   **Phase 3: Fully Relativistic (FR) Calculation**
    ```bash
    python ../../main.py --prefix NiTe2 run --mode FR --num_cores 4 --npool 1
    ```
*   **Phase 4: Plotting and Analysis**
    ```bash
    # Plot SR results
    python main.py --prefix NiTe2 plot --mode SR

    # Overlay SR and FR results for comparison
    python main.py --prefix NiTe2 plot --mode compare

    # Visualize the Fermi Surface (requires FermiSurfer)
    python main.py --prefix NiTe2 plot --mode fs
    ```
    *Note: Plotting axes limits can be adjusted using the optional flags `--emin`, `--emax`, and `--dos_max`.*

## Expected Output

The scripts will dynamically create the following directories within this example directory (if they do not already exist):
* `inputs/`: Contains the generated Quantum ESPRESSO input files.
* `tmp/`: Stores temporary data and charge density files during the calculation.
* `outputs/`: Contains the final data, including the parsed band structure files and plots.

Once the pipeline finishes, check the `outputs/` folder for the graphical representations of the NiTe2 bands and Density of States.

> **Note:** The complete execution of the pipeline may take several minutes depending on your hardware capabilities, the chosen calculation parameters, and the number of parallel processes utilized.

## Pre-computed Results

If you wish to examine the final output without running the calculations, the `outputs/` directory contains the pre-computed electronic structure of NiTe2.

### Band Structure and Density of States (Comparison)
Below is the final comparative plot showing both the Scalar-Relativistic (SR) and Fully Relativistic (FR) calculations, highlighting the effects of spin-orbit coupling on the Dirac nodes.

![NiTe2 Bands and DOS Comparison](reference_outputs/NiTe2_bands_dos_compare.png)