import argparse
import json
import os
import sys
from input_builder import InputBuilder
from simulation_runner import SimulationRunner
from simulation_analyzer import SimulationAnalyzer


def parse_k_path(k_path_list):
    """Parses the k-point path from the JSON configuration and returns the formatted string and labels."""
    num_points = len(k_path_list)
    k_lines = [str(num_points)]
    labels = []

    for pt in k_path_list:
        k_lines.append(f"{pt['coord']} {pt['points']}")
        labels.append(pt['label'])

    return "\n".join(k_lines), labels

def handle_setup(args):
    """
    Handles the structural optimization phase.
    Runs convergence tests for cutoff energy and k-points, followed by a vc-relax.
    Saves the optimized geometry to a state file.
    """
    print(f"Starting SETUP phase for {args.prefix}...")

    # Ensure working directories exist
    os.makedirs("./tmp", exist_ok=True)
    os.makedirs("./inputs", exist_ok=True)
    os.makedirs("./outputs", exist_ok=True)

    # Load infrastructure and material configurations
    with open('config.json', 'r', encoding='utf-8') as f:
        config = json.load(f)

    try:
        with open(f'{args.prefix}.json', 'r', encoding='utf-8') as f:
            mat_config = json.load(f)
    except FileNotFoundError:
        sys.exit(f"Critical Error: Material file '{args.prefix}.json' not found. Setup aborted.")

    # Initialize core operational classes
    builder = InputBuilder(
        prefix=args.prefix,
        output_folder="./inputs",
        global_params=config['global'],
        step_overrides=config.get('step_overrides', {}),
        atomic_species=mat_config['structure']['atomic_species_SR'],
        atomic_positions=mat_config['structure']['atomic_positions']
    )

    runner = SimulationRunner(qe_path=args.qe_path, output_folder="./outputs", num_cores=args.npool)
    analyzer = SimulationAnalyzer(prefix=args.prefix, dos_path="", bands_gnu_path="", bands_out_path="")

    # ---------------------------------------------------------
    # CONVERGENCE TESTS
    # ---------------------------------------------------------
    print("\n[1/3] Running cutoff energy convergence test...")
    ecut_list = mat_config['convergence']['ecut_list']
    nat = float(mat_config['structure']['nat'])

    ecut_in_files = builder.build_ecut_convergence(ecut_list=ecut_list,
                                                   k_points=mat_config['convergence']['kpoints_list'][1])
    ecut_out_files = runner.run_batch(ecut_in_files)

    if len(ecut_out_files) == len(ecut_list):
        best_ecut = analyzer.plot_ecut_convergence(ecut_out_files, ecut_list, nat=nat, threshold=0.001)
        if best_ecut:
            print(f"-> Ecut convergence reached. Setting ecutwfc = {best_ecut} Ry and ecutrho = {best_ecut * 10} Ry.")
            builder.system_params['ecutwfc'] = best_ecut
            builder.system_params['ecutrho'] = best_ecut * 10
    else:
        sys.exit("Critical Error: Ecut convergence batch failed. Accuracy cannot be guaranteed. Setup aborted.")

    print("\n[2/3] Running K-points convergence test...")
    kpoints_list = mat_config['convergence']['kpoints_list']
    kpoints_in_files = builder.build_kpoints_convergence(kpoints_list)
    kpoints_out_files = runner.run_batch(kpoints_in_files)

    if len(kpoints_out_files) == len(kpoints_list):
        best_k = analyzer.plot_kpoints_convergence(kpoints_out_files, kpoints_list, nat=nat, threshold=0.0001)
        print(f"-> K-points convergence reached for grid: {best_k}.")
    else:
        sys.exit("Critical Error: K-points convergence batch failed. Setup aborted.")

    # ---------------------------------------------------------
    # STRUCTURAL RELAXATION (vc-relax)
    # ---------------------------------------------------------
    print("\n[3/3] Running variable-cell structural relaxation...")
    vcrelax_in = builder.build_vcrelax_input(k_points=mat_config['convergence']['kpoints_scf'])
    vcrelax_out = runner.run_pw(vcrelax_in)

    if vcrelax_out is None:
        sys.exit("Critical Error: vc-relax calculation failed. Setup aborted.")

    # ---------------------------------------------------------
    # DATA EXTRACTION & PERSISTENCE
    # ---------------------------------------------------------
    try:
        relaxed_data = SimulationAnalyzer.extract_relaxed_structure(vcrelax_out)
    except Exception as e:
        sys.exit(f"Error parsing relaxed structure: {e}")

    state_file = f"{args.prefix}_state.json"
    state_data = {
        "cell_parameters": relaxed_data['cell_parameters'],
        "atomic_positions": relaxed_data['atomic_positions'],
        "optimal_ecut": builder.system_params.get('ecutwfc', 80)
    }

    with open(state_file, 'w', encoding='utf-8') as f:
        json.dump(state_data, f, indent=4)

    print(f"\nSetup successfully completed. Relaxed geometry saved to '{state_file}'.")

def handle_run(args):
    """
    Executes the main computational pipeline (SCF, DOS, BANDS).
    Requires a valid state file generated by the 'setup' command.
    """
    print(f"Starting RUN phase (Mode: {args.mode}) for system {args.prefix}...")

    os.makedirs("./tmp", exist_ok=True)
    os.makedirs("./inputs", exist_ok=True)
    os.makedirs("./outputs", exist_ok=True)

    # Load persistent state and configurations
    state_file = f"{args.prefix}_state.json"
    if not os.path.exists(state_file):
        sys.exit(f"Error: State file '{state_file}' not found. Run the 'setup' command first.")

    with open(state_file, 'r', encoding='utf-8') as f:
        state_data = json.load(f)

    with open('config.json', 'r', encoding='utf-8') as f:
        config = json.load(f)

    with open(f'{args.prefix}.json', 'r', encoding='utf-8') as f:
        mat_config = json.load(f)

    # Apply global parameters dynamically based on the relaxed state
    global_params = config['global']
    global_params['system_params']['ecutwfc'] = state_data.get('optimal_ecut', 80)
    global_params['system_params']['ecutrho'] = state_data.get('optimal_ecut', 80) * 10
    global_params['system_params']['ibrav'] = 0

    # Configure Scalar Relativistic (SR) or Fully Relativistic (FR) settings
    if args.mode == 'FR':
        global_params['system_params']['noncolin'] = True
        global_params['system_params']['lspinorb'] = True
        global_params['system_params']['nbnd'] = mat_config['structure'].get('nbnd_FR', 60)
        atomic_species = mat_config['structure']['atomic_species_FR']
    else:
        global_params['system_params'].pop('noncolin', None)
        global_params['system_params'].pop('lspinorb', None)
        global_params['system_params']['nbnd'] = mat_config['structure'].get('nbnd_SR', 30)
        atomic_species = mat_config['structure']['atomic_species_SR']

    dummy_positions = "Dummy 0.0 0.0 0.0"
    k_path_str, _ = parse_k_path(mat_config['k_path'])

    # Initialize the input builder
    builder = InputBuilder(
        prefix=args.prefix,
        output_folder="./inputs",
        global_params=global_params,
        step_overrides=config.get('step_overrides', {}),
        atomic_species=atomic_species,
        atomic_positions=dummy_positions
    )

    runner = SimulationRunner(qe_path=args.qe_path, output_folder="./outputs", num_cores=args.npool)

    # Inject the optimized structure
    builder.update_structure(
        new_cell_parameters=state_data['cell_parameters'],
        new_atomic_positions=state_data['atomic_positions']
    )

    builder.prefix = f"{args.prefix}_{args.mode}"
    current_prefix = builder.prefix

    # ---------------------------------------------------------
    # PIPELINE EXECUTION
    # ---------------------------------------------------------
    if args.step in ['all', 'scf']:
        print(f"\n[1/7] Executing base SCF calculation...")
        scf_in = builder.build_scf_input(k_points=mat_config['convergence']['kpoints_scf'], suffix=f"scf")
        if runner.run_pw(scf_in) is None:
            sys.exit("Critical Error: SCF calculation failed. Pipeline aborted.")

    if args.step in ['all', 'dos']:
        print(f"\n[2/7] Executing NSCF calculation for DOS (npool={args.npool})...")
        runner.prepare_branch_dir(base_tmp="./tmp", target_tmp="./tmp_dos", prefix=current_prefix)
        nscf_dos_in = builder.build_nscf_input(k_points=mat_config['convergence']['kpoints_dos'], nosym=False,
                                               outdir="./tmp_dos", suffix=f"nscf_dos")
        if runner.run_pw(nscf_dos_in, npool=args.npool) is None:
            sys.exit("Error: NSCF calculation for DOS failed. Pipeline aborted.")

        print(f"\n[3/7] Extracting DOS...")
        dos_in = builder.build_dos_input(deltae=0.01, outdir="./tmp_dos")
        if runner.run_dos(dos_in) is None:
            sys.exit("Error: DOS extraction (dos.x) failed. Pipeline aborted.")

    if args.step in ['all', 'fs']:
        print(f"\n[4/7] Executing NSCF calculation for Fermi Surface (npool={args.npool})...")
        runner.prepare_branch_dir(base_tmp="./tmp", target_tmp="./tmp_fs", prefix=current_prefix)
        nscf_fs_in = builder.build_nscf_input(k_points=mat_config['convergence']['kpoints_scf'], nosym=True,
                                              outdir="./tmp_fs", suffix=f"nscf_fs")
        if runner.run_pw(nscf_fs_in, npool=args.npool) is None:
            sys.exit("Error: NSCF calculation for Fermi Surface failed. Pipeline aborted.")

        print(f"\n[5/7] Extracting Fermi Surface...")
        fs_in = builder.build_fs_input(outdir="./tmp_fs")
        if runner.run_fs(fs_in) is None:
            sys.exit("Error: Fermi Surface extraction (fs.x) failed. Pipeline aborted.")

    if args.step in ['all', 'bands']:
        print(f"\n[6/7] Executing BANDS calculation (npool={args.npool})...")
        runner.prepare_branch_dir(base_tmp="./tmp", target_tmp="./tmp_bands", prefix=current_prefix)
        bands_in = builder.build_bands_input(k_points=k_path_str, outdir="./tmp_bands")
        if runner.run_pw(bands_in, npool=args.npool) is None:
            sys.exit("Error: Bands calculation failed. Pipeline aborted.")

        print(f"\n[7/7] Post-processing Bands...")
        bands_pp_in = builder.build_bands_pp_input(outdir="./tmp_bands")
        if runner.run_bands_pp(bands_pp_in) is None:
            sys.exit("Error: Bands post-processing (bands.x) failed. Pipeline aborted.")

    print(f"\nExecution ({args.step}) successfully completed for {args.mode} mode.")

def get_analyzer(prefix_base, mode, k_labels, output_folder="./outputs"):
    """Helper function to initialize the SimulationAnalyzer with correct paths."""
    prefix_mode = f"{prefix_base}_{mode}"
    return SimulationAnalyzer(
        prefix=prefix_mode,
        dos_path=os.path.join(output_folder, f"{prefix_mode}.dos"),
        bands_gnu_path=os.path.join(output_folder, f"{prefix_mode}_bands.dat.gnu"),
        bands_out_path=os.path.join(output_folder, f"{prefix_mode}_bands_pp.out"),
        fs_path=os.path.join(output_folder, f"{prefix_mode}.bxsf"),
        k_labels=k_labels
    )

def handle_plot(args):
    """
    Handles data visualization and plots generation.
    """
    print(f"Starting PLOT phase for mode: {args.mode} (System: {args.prefix})")

    try:
        with open(f'{args.prefix}.json', 'r', encoding='utf-8') as f:
            mat_config = json.load(f)
        _, k_labels = parse_k_path(mat_config['k_path'])

        energy_win = (args.emin, args.emax)

        if args.mode == 'compare':
            analyzer_sr = get_analyzer(args.prefix, 'SR', k_labels)
            analyzer_fr = get_analyzer(args.prefix, 'FR', k_labels)

            out_plot_name = f"{args.prefix}_bands_dos_compare.png"
            analyzer_fr.plot_bands_dos_compare(
                analyzer_sr, out_plot=out_plot_name, energy_window=energy_win, dos_max=args.dos_max
            )
            print(f"Comparative plot generated successfully: {out_plot_name}")

        elif args.mode == 'fs':
            analyzer = get_analyzer(args.prefix, 'SR', k_labels)
            print("Opening FermiSurfer to visualize the Fermi Surface...")
            analyzer.show_fermi_surface()

        else:
            analyzer = get_analyzer(args.prefix, args.mode, k_labels)
            fig = analyzer.plot_bands_dos(energy_window=energy_win, dos_max=args.dos_max)
            out_plot_name = f"{args.prefix}_{args.mode}_bands_dos.png"
            fig.savefig(out_plot_name, dpi=300, bbox_inches='tight')
            print(f"Single plot ({args.mode}) generated successfully: {out_plot_name}")

    except Exception as e:
        sys.exit(f"\n[!] ERROR during plotting/visualization: {e}")

def main(args_list=None):
    """
    Main entry point for the pipeline.
    Parses command line arguments or accepts a manual list of arguments.
    """
    parser = argparse.ArgumentParser(
        description="Automated Generic Pipeline for DFT Simulations and Post-Processing (Quantum ESPRESSO).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--prefix', type=str, required=True,
                        help="Material prefix (must match the <prefix>.json configuration file).")
    parser.add_argument('--qe_path', type=str, default='',
                        help="Absolute path to QE binaries (e.g., /usr/local/bin). Leave empty to use system PATH.")

    subparsers = parser.add_subparsers(dest="command", required=True, help="Simulation phase to execute")

    parser_setup = subparsers.add_parser('setup', help="Executes convergence tests and structural relaxation.")
    parser_setup.add_argument('--npool', type=int, default=4, help="Number of MPI processes.")

    parser_run = subparsers.add_parser('run', help="Executes the main computational pipeline (SCF, NSCF, BANDS).")
    parser_run.add_argument('--mode', type=str, choices=['SR', 'FR'], required=True,
                            help="Physical mode: SR (Scalar Relativistic) or FR (Fully Relativistic).")
    parser_run.add_argument('--npool', type=int, default=4, help="Number of k-point pools for parallelization.")
    parser_run.add_argument('--step', type=str, choices=['all', 'scf', 'dos', 'fs', 'bands'], default='all',
                            help="Execute the entire pipeline ('all') or a specific logical block.")

    parser_plot = subparsers.add_parser('plot', help="Generates plots and visualizations from output data.")
    parser_plot.add_argument('--mode', type=str, choices=['SR', 'FR', 'compare', 'fs'], required=True,
                             help="Plotting mode: SR, FR, compare (overlay), or fs (FermiSurfer).")
    parser_plot.add_argument('--dos_max', type=float, default=None, help="Custom upper limit for the DOS X-axis.")
    parser_plot.add_argument('--emin', type=float, default=-2.0,
                             help="Minimum energy (eV) relative to the Fermi Level.")
    parser_plot.add_argument('--emax', type=float, default=2.0, help="Maximum energy (eV) relative to the Fermi Level.")

    # Parse either the provided list (from external scripts) or sys.argv (from terminal)
    args = parser.parse_args(args_list)

    if args.command == 'setup':
        handle_setup(args)
    elif args.command == 'run':
        handle_run(args)
    elif args.command == 'plot':
        handle_plot(args)

def run_pipeline(prefix):
    """
    Programmatic entry point for external wrapper scripts.
    Executes sequentially: setup -> run (SR) -> plot (SR).
    """
    print(f"\n{'=' * 50}\nStarting Automated Pipeline for: {prefix}\n{'=' * 50}")

    main(["--prefix", prefix, "setup"])
    main(["--prefix", prefix, "run", "--mode", "SR"])
    main(["--prefix", prefix, "plot", "--mode", "SR"])

if __name__ == "__main__":
    main()