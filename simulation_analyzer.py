import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess


class SimulationAnalyzer:
    """
    Parses output files from Quantum ESPRESSO and generates analytical plots.
    Handles data extraction for convergence tests, structural relaxation,
    Density of States (DOS), Band Structure, and Fermi Surfaces.
    """

    def __init__(self, prefix, dos_path, bands_gnu_path, bands_out_path, fs_path=None, k_labels=None):
        self.prefix = prefix

        self.dos_path = dos_path
        self.bands_path = bands_gnu_path
        self.bands_pp_path = bands_out_path
        self.fs_path = fs_path

        self.k_labels = k_labels

        self.fermi_energy = None
        self.dos_energy = None
        self.dos_values = None
        self.k_values_matrix = []
        self.bands_energy_matrix = []
        self.high_simm_points = []

    def plot_ecut_convergence(self, out_files, ecut_list, nat, threshold=0.001):
        """
        Parses output files from a cutoff energy convergence test,
        identifies the convergence point, and generates a plot.
        """
        energies_per_atom = []
        for out_file in out_files:
            tot_e = self.extract_total_energy(out_file)
            if tot_e is None:
                raise ValueError(f"Total energy not found in file {out_file}")
            energies_per_atom.append(tot_e / nat)

        # Determine the convergence cutoff
        convergence_cutoff = None
        for i in range(1, len(energies_per_atom)):
            delta = abs(energies_per_atom[i] - energies_per_atom[i - 1])
            if delta <= threshold:
                convergence_cutoff = ecut_list[i]
                break

        if convergence_cutoff is None and len(ecut_list) > 0:
            convergence_cutoff = ecut_list[-1]

        # Generate the plot
        plt.figure(figsize=(8, 5))
        plt.plot(ecut_list, energies_per_atom, marker='o', linestyle='-', color='b', label='Energy per atom')

        if convergence_cutoff:
            plt.axvline(x=convergence_cutoff, color='red', linestyle='--', linewidth=1.5,
                        label=f'Convergence: {convergence_cutoff} Ry')

        plt.title(f'Plane Waves Saturation: {self.prefix}')
        plt.xlabel('Cutoff Energy (Ry)')
        plt.ylabel('Energy per atom (Ry)')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()

        # Save the plot
        plot_path = f"{self.prefix}_convergence_ecut.png"
        plt.savefig(plot_path, dpi=300)
        print(f"[Analyzer] Ecut convergence plot saved to {plot_path}")

        return convergence_cutoff

    def plot_kpoints_convergence(self, out_files, kpoints_list, nat, threshold=0.001):
        """
        Parses output files from a k-points convergence test,
        identifies the optimal grid, and generates a plot.
        """
        energies_per_atom = []
        for out_file in out_files:
            tot_e = self.extract_total_energy(out_file)
            if tot_e is None:
                raise ValueError(f"Total energy not found in file {out_file}")
            energies_per_atom.append(tot_e / nat)

        # Generate plot labels (e.g., '6x6x4')
        k_labels = ["x".join(k.split()[:3]) for k in kpoints_list]

        # Determine the convergence index
        convergence_idx = None
        for i in range(1, len(energies_per_atom)):
            delta = abs(energies_per_atom[i] - energies_per_atom[i - 1])
            if delta <= threshold:
                convergence_idx = i
                break

        # Fallback to the densest grid if the threshold is not met
        if convergence_idx is None and len(energies_per_atom) > 0:
            convergence_idx = len(energies_per_atom) - 1

        # Generate the plot
        plt.figure(figsize=(8, 5))
        plt.plot(k_labels, energies_per_atom, marker='o', linestyle='-', color='blue', label='Energy per atom')

        if convergence_idx is not None:
            plt.axvline(x=convergence_idx, color='red', linestyle='--', linewidth=1.5,
                        label=f'Convergence: {k_labels[convergence_idx]}')

        plt.title(f'K-points Convergence: {self.prefix}')
        plt.xlabel('K-points Grid ($N\\times N\\times M$)')
        plt.ylabel('Energy per Atom (Ry)')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()

        # Save the plot
        plot_path = f"{self.prefix}_convergence_kpoints.png"
        plt.savefig(plot_path, dpi=300)
        print(f"[Analyzer] K-points convergence plot saved to {plot_path}")

        return k_labels[convergence_idx] if convergence_idx is not None else None

    @staticmethod
    def extract_relaxed_structure(vcrelax_out_path):
        """
        Extracts final cell parameters and atomic positions from a vc-relax output.
        If convergence is reached at step 0 (no movement needed), it reconstructs
        the geometry from the corresponding input file.
        """
        if not os.path.exists(vcrelax_out_path):
            raise FileNotFoundError(f"File {vcrelax_out_path} not found.")

        with open(vcrelax_out_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

        # 1. Locate the final geometry blocks in the output
        last_cell_idx = -1
        last_atoms_idx = -1

        for i, line in enumerate(lines):
            if line.startswith("CELL_PARAMETERS"):
                last_cell_idx = i
            elif line.startswith("ATOMIC_POSITIONS"):
                last_atoms_idx = i

        # 2. Extract standard blocks if present
        if last_cell_idx != -1 and last_atoms_idx != -1:
            cell_params = []
            atomic_positions = []

            for i in range(last_cell_idx, last_cell_idx + 4):
                cell_params.append(lines[i].strip())

            i = last_atoms_idx
            while i < len(lines):
                line = lines[i].strip()
                if i > last_atoms_idx and (not line or "End final coordinates" in line or "writing output" in line):
                    break
                atomic_positions.append(line)
                i += 1

            return {
                'cell_parameters': "\n".join(cell_params),
                'atomic_positions': "\n".join(atomic_positions)
            }

        # 3. Fallback: Step 0 convergence (construct initial geometry)
        print(f"[Analyzer] Step 0 convergence detected. Reconstructing initial geometry in absolute coordinates...")

        # 3a. Extract 'alat' and crystal axes from the output
        alat = None
        cell_params = ["CELL_PARAMETERS (bohr)"]

        for i, line in enumerate(lines):
            if "lattice parameter (alat)" in line:
                alat = float(line.split('=')[1].split('a.u.')[0].strip())
            elif "crystal axes: (cart. coord. in units of alat)" in line:
                if alat is None:
                    raise ValueError("Cannot find 'alat' in the output before the crystal axes.")
                for j in range(1, 4):
                    # Parse vector components
                    right_side = lines[i + j].split('=')[1]
                    vec_raw = right_side.replace('(', '').replace(')', '').strip()
                    v_x, v_y, v_z = map(float, vec_raw.split())

                    # Scale by alat to convert to Bohr
                    cell_params.append(f"{v_x * alat:15.10f} {v_y * alat:15.10f} {v_z * alat:15.10f}")
                break

        # 3b. Retrieve ATOMIC_POSITIONS from the corresponding input file
        vcrelax_in_path = vcrelax_out_path.replace("/outputs/", "/inputs/").replace(".out", ".in")
        if not os.path.exists(vcrelax_in_path):
            vcrelax_in_path = vcrelax_out_path.rsplit('.', 1)[0] + ".in"

        atomic_positions = []
        if os.path.exists(vcrelax_in_path):
            with open(vcrelax_in_path, 'r', encoding='utf-8') as f:
                in_lines = f.readlines()

            i = 0
            while i < len(in_lines):
                line = in_lines[i].strip()
                if line.startswith("ATOMIC_POSITIONS"):
                    atomic_positions.append(line)
                    i += 1
                    while i < len(in_lines) and in_lines[i].strip() and not in_lines[i].strip().startswith("K_POINTS"):
                        atomic_positions.append(in_lines[i].strip())
                        i += 1
                    break
                i += 1

        return {
            'cell_parameters': "\n".join(cell_params),
            'atomic_positions': "\n".join(atomic_positions)
        }

    def extract_fermi_energy(self):
        with open(self.dos_path, 'r', encoding='utf-8') as f:
            for line in f:
                if "EFermi =" in line:
                    self.fermi_energy = float(line.rsplit(maxsplit=2)[1])
                    break

        if self.fermi_energy is None:
            raise ValueError(f"Fermi energy not found in file {self.dos_path}")

    def extract_dos_data(self):
        if self.fermi_energy is None:
            self.extract_fermi_energy()

        raw_energy, self.dos_values = np.loadtxt(self.dos_path, unpack=True, usecols=(0, 1))
        self.dos_energy = raw_energy - self.fermi_energy

    def extract_bands_data(self):
        i = 0
        flag = False

        self.k_values_matrix = [[]]
        self.bands_energy_matrix = [[]]

        if self.fermi_energy is None:
            self.extract_fermi_energy()

        with open(self.bands_path, 'r', encoding='utf8') as f:
            for line in f:
                if line == '\n' and not flag:
                    flag = True
                    i += 1
                    self.k_values_matrix.append([])
                    self.bands_energy_matrix.append([])
                elif line != '\n':
                    flag = False
                    row = line.split()
                    self.k_values_matrix[i].append(float(row[0]))
                    band_energy = float(row[1]) - self.fermi_energy
                    self.bands_energy_matrix[i].append(band_energy)

        if not self.k_values_matrix[-1]:
            self.k_values_matrix.pop(-1)
            self.bands_energy_matrix.pop(-1)

    def extract_high_simm_points(self):
        self.high_simm_points = []

        with open(self.bands_pp_path, 'r', encoding='utf-8') as f:
            for line in f:
                if 'Reading collected, re-writing distributed wavefunctions' in line:
                    line = next(f)
                    while line != '\n':
                        x_coord = float(line.rsplit(maxsplit=1)[1])
                        self.high_simm_points.append(x_coord)
                        line = next(f)
                    break

        if self.k_labels is not None:
            if len(self.k_labels) != len(self.high_simm_points):
                raise ValueError(
                    f"K-path mismatch: {len(self.k_labels)} labels provided, "
                    f"but {len(self.high_simm_points)} coordinates found in {self.bands_pp_path}"
                )

    def plot_bands_dos(self, energy_window=(-2, 2), dos_max=None):
        if not self.bands_energy_matrix:
            self.extract_bands_data()
        if not self.high_simm_points:
            self.extract_high_simm_points()
        if self.dos_energy is None:
            self.extract_dos_data()

        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, gridspec_kw={'width_ratios': [2, 1]},
                                       figsize=(12, 8), dpi=300)

        # --- LEFT PANEL (BANDS) ---
        for i in range(len(self.k_values_matrix)):
            ax1.plot(self.k_values_matrix[i], self.bands_energy_matrix[i], color='black', linewidth=1.2)

        ax1.axhline(0, color='red', linestyle='--', alpha=0.6)

        ax1.set_xticks(self.high_simm_points)
        if self.k_labels is not None:
            ax1.set_xticklabels(self.k_labels)

        for pt in self.high_simm_points:
            ax1.axvline(pt, color='black', linestyle='--', linewidth=0.5, alpha=0.4)

        ax1.set_xlim(self.k_values_matrix[0][0], self.k_values_matrix[0][-1])
        ax1.set_ylim(energy_window[0], energy_window[1])
        ax1.set_ylabel('Energy - E_F (eV)')
        ax1.tick_params(direction='in', top=True, right=True)

        # --- RIGHT PANEL (DOS) ---
        ax2.plot(self.dos_values, self.dos_energy, color='blue')
        ax2.fill_betweenx(self.dos_energy, 0, self.dos_values, where=(self.dos_energy < 0),
                          color='lightblue', alpha=0.5)

        ax2.axhline(0, color='red', linestyle='--', alpha=0.6)

        # Apply custom limit or calculate automatically based on the energy window
        if dos_max is not None:
            ax2.set_xlim(0, dos_max)
        else:
            mask = (self.dos_energy >= energy_window[0]) & (self.dos_energy <= energy_window[1])
            dos_in_window = self.dos_values[mask]

            if len(dos_in_window) > 0:
                calc_max = float(np.max(dos_in_window)) * 1.05
            else:
                calc_max = float(np.max(self.dos_values)) * 1.05

            ax2.set_xlim(0, calc_max)

        ax2.set_xlabel('DOS (states / eV)')
        ax2.tick_params(direction='in', top=True, right=True)

        # --- GLOBAL LAYOUT ---
        fig.suptitle(f"Band Structure and Density of States ({self.prefix})", fontsize=16, fontweight='bold')
        fig.tight_layout()

        return fig

    def plot_bands_dos_compare(self, other_analyzer, energy_window=(-4.0, 30.0), out_plot="bands_dos_compare.png",
                               dos_max=None):
        """
        Overlays the data of the current instance (Foreground, e.g., FR)
        on top of another analyzer instance (Background, e.g., SR).
        """
        from matplotlib.lines import Line2D

        # 1. Extract data for the current instance (Foreground)
        if not self.bands_energy_matrix: self.extract_bands_data()
        if not self.high_simm_points: self.extract_high_simm_points()
        if self.dos_energy is None: self.extract_dos_data()

        # 2. Extract data for the other instance (Background)
        if not other_analyzer.bands_energy_matrix: other_analyzer.extract_bands_data()
        if not other_analyzer.high_simm_points: other_analyzer.extract_high_simm_points()
        if other_analyzer.dos_energy is None: other_analyzer.extract_dos_data()

        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, gridspec_kw={'width_ratios': [2, 1]},
                                       figsize=(12, 8), dpi=300)

        # --- LEFT PANEL (BANDS) ---
        # Plot Background (other_analyzer / SR)
        for i in range(len(other_analyzer.k_values_matrix)):
            ax1.plot(other_analyzer.k_values_matrix[i], other_analyzer.bands_energy_matrix[i],
                     color='royalblue', linewidth=1.0, linestyle='--', alpha=0.7)

        # Plot Foreground (self / FR)
        for i in range(len(self.k_values_matrix)):
            ax1.plot(self.k_values_matrix[i], self.bands_energy_matrix[i], color='black', linewidth=1.2)

        custom_lines = [Line2D([0], [0], color='royalblue', lw=1.0, linestyle='--'),
                        Line2D([0], [0], color='black', lw=1.2)]
        ax1.legend(custom_lines, ['SR (Scalar)', 'FR (Spin-Orbit)'], loc='upper right', framealpha=0.9)

        ax1.axhline(0, color='red', linestyle='--', alpha=0.5, linewidth=1.0)
        ax1.set_xticks(self.high_simm_points)
        if self.k_labels is not None:
            ax1.set_xticklabels(self.k_labels)

        for pt in self.high_simm_points:
            ax1.axvline(pt, color='black', linestyle='--', linewidth=0.5, alpha=0.4)

        ax1.set_xlim(self.k_values_matrix[0][0], self.k_values_matrix[0][-1])
        ax1.set_ylim(energy_window[0], energy_window[1])
        ax1.set_ylabel('Energy - E_F (eV)')
        ax1.tick_params(direction='in', top=True, right=True)

        # --- RIGHT PANEL (DOS) ---
        # DOS Background (other_analyzer / SR)
        ax2.plot(other_analyzer.dos_values, other_analyzer.dos_energy, color='royalblue', linestyle='--', linewidth=1.0,
                 alpha=0.8)
        ax2.fill_betweenx(other_analyzer.dos_energy, 0, other_analyzer.dos_values,
                          where=(other_analyzer.dos_energy < 0),
                          facecolor='none', edgecolor='royalblue', hatch='///', alpha=0.3)

        # DOS Foreground (self / FR)
        ax2.plot(self.dos_values, self.dos_energy, color='black', linewidth=1.2)
        ax2.fill_betweenx(self.dos_energy, 0, self.dos_values, where=(self.dos_energy < 0),
                          color='gray', alpha=0.2)

        ax2.axhline(0, color='red', linestyle='--', alpha=0.5, linewidth=1.0)

        # Apply custom limit or calculate automatically based on the energy window for both datasets
        if dos_max is not None:
            ax2.set_xlim(0, dos_max)
        else:
            # Foreground (FR) max
            mask_self = (self.dos_energy >= energy_window[0]) & (self.dos_energy <= energy_window[1])
            dos_in_window_self = self.dos_values[mask_self]
            max_self = float(np.max(dos_in_window_self)) if len(dos_in_window_self) > 0 else float(
                np.max(self.dos_values))

            # Background (SR) max
            mask_other = (other_analyzer.dos_energy >= energy_window[0]) & (
                        other_analyzer.dos_energy <= energy_window[1])
            dos_in_window_other = other_analyzer.dos_values[mask_other]
            max_other = float(np.max(dos_in_window_other)) if len(dos_in_window_other) > 0 else float(
                np.max(other_analyzer.dos_values))

            # Set limit to the absolute maximum + 5% margin
            calc_max = max(max_self, max_other) * 1.05
            ax2.set_xlim(0, calc_max)

        ax2.set_xlabel('DOS (states / eV)')
        ax2.tick_params(direction='in', top=True, right=True)

        # --- GLOBAL LAYOUT ---
        # Extract base prefix for the title (e.g., 'NiTe2' from 'NiTe2_FR')
        material_name = self.prefix.split('_')[0]
        fig.suptitle(f"Band Structure and DOS of {material_name}: SR vs FR Comparison", fontsize=16)

        fig.savefig(out_plot, dpi=300, bbox_inches='tight')
        plt.close(fig)
        return out_plot

    def fix_bxsf_format(self):
        if not os.path.exists(self.fs_path):
            raise FileNotFoundError(f"File {self.fs_path} not found for Fermi surface.")

        with open(self.fs_path, 'r', encoding='utf-8') as f:
            content = f.read()

        if 'BANDGRID_3D_BANDS' in content:
            content = content.replace('BANDGRID_3D_BANDS', 'BEGIN_BANDGRID_3D')

            with open(self.fs_path, 'w', encoding='utf-8') as f:
                f.write(content)

    def show_fermi_surface(self):
        self.fix_bxsf_format()

        try:
            subprocess.Popen(['fermisurfer', self.fs_path])
        except FileNotFoundError:
            raise RuntimeError(
                "The 'fermisurfer' executable was not found. "
                "Ensure it is globally accessible in the system PATH."
            )

    @staticmethod
    def extract_total_energy(filepath):
        """Extracts the total energy (in Rydberg) from a pw.x output file."""
        if not os.path.exists(filepath):
            return None

        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith("!") and "total energy" in line:
                    return float(line.split()[4])

        return None