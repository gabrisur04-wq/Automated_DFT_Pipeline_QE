import os


class InputBuilder:
    """
    Constructs and formats input files for Quantum ESPRESSO (pw.x, dos.x, bands.x, fs.x).
    It manages the dynamic injection of parameters, structural data, and specific overrides.
    """

    def __init__(self, prefix, output_folder, global_params, step_overrides,
                 atomic_species, atomic_positions, data_folder="./outputs"):

        self.prefix = prefix
        self.output_folder = output_folder
        self.data_folder = data_folder

        self.global_params = global_params
        self.step_overrides = step_overrides

        # Isolate base dictionaries to prevent modifying the original configuration
        self.control_params = global_params.get('control_params', {}).copy()
        self.system_params = global_params.get('system_params', {}).copy()
        self.electrons_params = global_params.get('electrons_params', {}).copy()

        self.atomic_species = atomic_species
        self.atomic_positions = atomic_positions
        self.cell_parameters = None

        os.makedirs(self.output_folder, exist_ok=True)
        os.makedirs(self.data_folder, exist_ok=True)

    def _apply_overrides(self, step_name):
        """Merges global parameters with step-specific overrides defined in the JSON configuration."""
        ctrl = self.control_params.copy()
        sys = self.system_params.copy()
        elec = self.electrons_params.copy()

        overrides = self.step_overrides.get(step_name, {})

        ctrl.update(overrides.get('control', {}))
        sys.update(overrides.get('system', {}))
        elec.update(overrides.get('electrons', {}))

        return ctrl, sys, elec

    def _format_namelist(self, name, params_dict):
        """Translates a Python dictionary into a formatted Fortran Namelist for Quantum ESPRESSO."""
        if params_dict is None:
            params_dict = {}

        if name.upper() == "CONTROL":
            params_dict['prefix'] = self.prefix

        lines = [f"&{name}"]

        for key, value in params_dict.items():
            # Translate Python native booleans to Fortran syntax
            if isinstance(value, bool):
                value = '.true.' if value else '.false.'
            # Format strings safely, avoiding double quotes for Fortran booleans
            elif isinstance(value, str):
                is_already_quoted = value.startswith("'") or value.startswith('"')
                is_fortran_bool = value.lower() in ['.true.', '.false.', 'true', 'false']
                if not (is_already_quoted or is_fortran_bool):
                    value = f"'{value}'"

            lines.append(f"    {key} = {value}")

        lines.append("/")
        return "\n".join(lines)

    def update_structure(self, new_cell_parameters, new_atomic_positions):
        """
        Updates the crystal structure for subsequent calculations.
        Forces ibrav=0, clears legacy celldm parameters, and injects explicit cell vectors.
        """
        self.system_params['ibrav'] = 0
        for i in range(1, 7):
            self.system_params.pop(f'celldm({i})', None)

        self.cell_parameters = new_cell_parameters
        self.atomic_positions = new_atomic_positions

    def _build_structure_blocks(self):
        """Constructs the ATOMIC_SPECIES, CELL_PARAMETERS, and ATOMIC_POSITIONS blocks."""
        blocks = ["ATOMIC_SPECIES", self.atomic_species.strip()]

        if self.system_params.get('ibrav') == 0:
            if self.cell_parameters and self.cell_parameters.strip():
                if not self.cell_parameters.lstrip().startswith("CELL_PARAMETERS"):
                    blocks.append("CELL_PARAMETERS (alat)")
                blocks.append(self.cell_parameters.strip())
            else:
                raise ValueError("Critical Error: ibrav=0 requires 'cell_parameters', but the block is empty.")
        else:
            if self.cell_parameters and self.cell_parameters.strip():
                blocks.append(self.cell_parameters.strip())

        if self.atomic_positions and self.atomic_positions.lstrip().startswith("ATOMIC_POSITIONS"):
            blocks.append(self.atomic_positions.strip())
        else:
            blocks.append("ATOMIC_POSITIONS (crystal)")
            if self.atomic_positions:
                blocks.append(self.atomic_positions.strip())
            else:
                raise ValueError("Critical Error: 'atomic_positions' is empty in the Builder.")

        return blocks

    def build_vcrelax_input(self, k_points, suffix="vcrelax"):
        ctrl, sys, elec = self._apply_overrides('vcrelax')
        ctrl['prefix'] = self.prefix

        control_str = self._format_namelist("CONTROL", ctrl)
        system_str = self._format_namelist("SYSTEM", sys)
        electrons_str = self._format_namelist("ELECTRONS", elec)

        # Retrieve specific parameters for ions and cell dynamics if defined in JSON
        ions_overrides = self.step_overrides.get('vcrelax', {}).get('ions', {'ion_dynamics': 'bfgs'})
        cell_overrides = self.step_overrides.get('vcrelax', {}).get('cell',
                                                                    {'cell_dynamics': 'bfgs', 'cell_dofree': 'all'})

        ions_str = self._format_namelist("IONS", ions_overrides)
        cell_str = self._format_namelist("CELL", cell_overrides)

        file_content = [control_str, system_str, electrons_str, ions_str, cell_str]
        file_content.extend(self._build_structure_blocks())
        file_content.extend(["K_POINTS (automatic)", k_points.strip()])

        final_string = "\n".join(file_content)
        file_path = os.path.join(self.output_folder, f"{self.prefix}_{suffix}.in")

        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(final_string + "\n")

        return file_path

    def build_scf_input(self, k_points, suffix="scf"):
        self.control_params['calculation'] = 'scf'
        self.control_params['prefix'] = self.prefix

        control_str = self._format_namelist("CONTROL", self.control_params)
        system_str = self._format_namelist("SYSTEM", self.system_params)
        electrons_str = self._format_namelist("ELECTRONS", self.electrons_params)

        file_content = [
            control_str,
            system_str,
            electrons_str,
            *self._build_structure_blocks(),
            "K_POINTS (automatic)",
            k_points.strip()
        ]
        final_string = "\n".join(file_content)

        file_path = os.path.join(self.output_folder, f"{self.prefix}_{suffix}.in")

        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(final_string + "\n")

        return file_path

    def build_ecut_convergence(self, ecut_list, k_points):
        """Generates a batch of SCF input files to test wavefunction cutoff energy convergence."""
        generated_files = []
        original_ecut = self.system_params.get('ecutwfc')
        original_ecutrho = self.system_params.get('ecutrho')

        for ecut in ecut_list:
            # Dynamically scale charge density cutoff alongside wavefunction cutoff
            self.system_params['ecutwfc'] = ecut
            self.system_params['ecutrho'] = ecut * 10

            file_path = self.build_scf_input(k_points, suffix=f"scf_ecut_{ecut}")
            generated_files.append(file_path)

        # Restore original configuration
        if original_ecut is not None:
            self.system_params['ecutwfc'] = original_ecut
        else:
            self.system_params.pop('ecutwfc', None)

        if original_ecutrho is not None:
            self.system_params['ecutrho'] = original_ecutrho
        else:
            self.system_params.pop('ecutrho', None)

        return generated_files

    def build_kpoints_convergence(self, kpoints_list):
        """Generates a batch of SCF input files to test k-points grid convergence."""
        generated_files = []
        for k_mesh in kpoints_list:
            parts = k_mesh.split()
            k_name = "x".join(parts[:3])  # Format output filename (e.g., '6x6x4')
            file_path = self.build_scf_input(k_points=k_mesh, suffix=f"scf_k_{k_name}")
            generated_files.append(file_path)

        return generated_files

    def build_nscf_input(self, k_points, nosym=False, outdir=None, suffix="nscf"):
        ctrl, sys, elec = self._apply_overrides('nscf')

        ctrl['prefix'] = self.prefix
        ctrl['outdir'] = outdir if outdir else self.control_params.get('outdir', './tmp/')

        # Programmatically remove smearing and degauss to prevent conflicts with 'tetrahedra' occupations
        sys.pop('smearing', None)
        sys.pop('degauss', None)

        # Disable symmetry reduction to ensure a uniform k-grid for specific post-processing (e.g., fs.x)
        if nosym:
            sys['nosym'] = '.true.'
        else:
            sys.pop('nosym', None)

        control_str = self._format_namelist("CONTROL", ctrl)
        system_str = self._format_namelist("SYSTEM", sys)
        electrons_str = self._format_namelist("ELECTRONS", elec)

        file_content = [
            control_str,
            system_str,
            electrons_str,
            *self._build_structure_blocks(),
            "K_POINTS (automatic)",
            k_points.strip()
        ]
        final_string = "\n".join(file_content)
        file_path = os.path.join(self.output_folder, f"{self.prefix}_{suffix}.in")

        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(final_string + "\n")

        return file_path

    def build_dos_input(self, deltae=0.01, emin=None, emax=None, outdir=None):
        outdir = outdir if outdir else self.control_params.get('outdir', './tmp/')
        dos_filename = os.path.join(self.data_folder, f"{self.prefix}.dos")

        params = {
            'prefix': self.prefix,
            'outdir': outdir,
            'fildos': dos_filename,
            'DeltaE': deltae
        }
        if emin is not None:
            params['Emin'] = emin
        if emax is not None:
            params['Emax'] = emax

        file_content = self._format_namelist("DOS", params)
        file_path = os.path.join(self.output_folder, f"{self.prefix}_dos.in")

        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(file_content + "\n")

        return file_path

    def build_fs_input(self, outdir=None):
        outdir = outdir if outdir else self.control_params.get('outdir', './tmp/')
        fs_filename = os.path.join(self.data_folder, self.prefix)

        params = {'prefix': self.prefix, 'outdir': outdir, 'file_fs': fs_filename}
        file_content = self._format_namelist("FERMI", params)
        file_path = os.path.join(self.output_folder, f"{self.prefix}_fs.in")

        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(file_content + "\n")

        return file_path

    def build_bands_input(self, k_points, outdir=None):
        ctrl, sys, elec = self._apply_overrides('bands')

        ctrl['prefix'] = self.prefix
        ctrl['outdir'] = outdir if outdir else self.control_params.get('outdir', './tmp/')

        # --- K-POINTS CONSISTENCY CHECK ---
        # Ensures the number of k-points declared in the header matches the provided coordinate lines.
        if isinstance(k_points, str):
            lines = [line.strip() for line in k_points.strip().split('\n') if line.strip()]
            if lines:
                try:
                    num_k_declared = int(lines[0])
                    num_k_lines = len(lines) - 1
                    if num_k_declared != num_k_lines:
                        raise ValueError(
                            f"K-points mismatch in BANDS calculation: "
                            f"{num_k_declared} points declared in header, but {num_k_lines} coordinates provided."
                        )
                except ValueError as e:
                    if "invalid literal" not in str(e):
                        raise e

        control_str = self._format_namelist("CONTROL", ctrl)
        system_str = self._format_namelist("SYSTEM", sys)
        electrons_str = self._format_namelist("ELECTRONS", elec)

        file_content = [
            control_str,
            system_str,
            electrons_str,
            *self._build_structure_blocks(),
            "K_POINTS (crystal_b)",
            k_points.strip()
        ]
        final_string = "\n".join(file_content)

        file_path = os.path.join(self.output_folder, f"{self.prefix}_bands.in")
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(final_string + "\n")

        return file_path

    def build_bands_pp_input(self, outdir=None):
        outdir = outdir if outdir else self.control_params.get('outdir', './tmp/')
        bands_filename = os.path.join(self.data_folder, f"{self.prefix}_bands.dat")

        params = {
            'prefix': self.prefix,
            'outdir': outdir,
            'filband': bands_filename,
            'lsym': '.false.'
        }
        file_content = self._format_namelist("BANDS", params)
        file_path = os.path.join(self.output_folder, f"{self.prefix}_bands_pp.in")

        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(file_content + "\n")

        return file_path