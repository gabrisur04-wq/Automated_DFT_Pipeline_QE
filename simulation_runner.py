import os
import shutil
import subprocess


class SimulationRunner:
    """
    Handles the execution of Quantum ESPRESSO binaries (pw.x, dos.x, bands.x, fs.x).
    Manages MPI parallelization, working directories, and verifies the physical integrity
    of the calculations upon completion.
    """

    def __init__(self, output_folder="./outputs", num_cores=1, qe_path=""):
        """
        Initializes the SimulationRunner.

        :param output_folder: Directory to save output files (.out).
        :param num_cores: Number of MPI processes to use.
        :param qe_path: Absolute path to QE binaries (e.g., '/usr/local/bin/'). Leave empty to use system PATH.
        """
        self.output_folder = output_folder
        self.num_cores = num_cores
        self.qe_path = qe_path

        os.makedirs(self.output_folder, exist_ok=True)

    def prepare_branch_dir(self, base_tmp="./tmp", target_tmp="./tmp_dos", prefix=None):
        """
        Copies the .save directory from the base outdir to a new target directory.
        This isolates the wavefunction data, allowing for independent execution of parallel
        calculation branches (e.g., DOS and BANDS) without data collision.
        """
        save_folder_name = f"{prefix}.save"
        src_save = os.path.join(base_tmp, save_folder_name)
        dst_save = os.path.join(target_tmp, save_folder_name)

        if not os.path.exists(src_save):
            raise FileNotFoundError(f"Cannot find the source directory: {src_save}")

        if os.path.exists(target_tmp):
            shutil.rmtree(target_tmp)
        os.makedirs(target_tmp, exist_ok=True)

        shutil.copytree(src_save, dst_save)
        print(f"Copied .save directory from '{base_tmp}' to '{target_tmp}'")

    def _execute_command(self, executable, input_file_path, npool=None, serial=False):
        """
        Internal method to execute a Quantum ESPRESSO binary on a specific input file.
        If serial=True, it forces execution without mpirun (required for post-processing tools).
        """
        base_name = os.path.basename(input_file_path).replace('.in', '')
        output_file_path = os.path.join(self.output_folder, f"{base_name}.out")

        exec_cmd = os.path.join(self.qe_path, executable)

        # Run with MPI if cores > 1 and serial mode is not explicitly enforced
        if self.num_cores > 1 and not serial:
            cmd = ["mpirun", "-np", str(self.num_cores), exec_cmd]
            if npool is not None:
                cmd.extend(["-npool", str(npool)])
            cmd.extend(["-in", input_file_path])
        else:
            cmd = [exec_cmd, "-in", input_file_path]

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"

        print(f"Executing {executable} for: {input_file_path}...")

        with open(output_file_path, 'w', encoding='utf-8') as out_file:
            process = subprocess.run(
                cmd,
                stdout=out_file,
                stderr=subprocess.PIPE,
                text=True,
                env=env
            )

        if process.returncode != 0:
            print(f"ERROR during the execution of {executable}!")
            print(f"Error details: {process.stderr}")
            return None

        return output_file_path

    def _check_job_done(self, output_file_path):
        """Verifies if the .out file contains the 'JOB DONE.' confirmation string."""
        if not os.path.exists(output_file_path):
            return False

        with open(output_file_path, 'r', encoding='utf-8', errors='ignore') as f:
            # Optimization: Instead of loading the entire file into memory,
            # we only scan the last 50 lines in reverse order.
            lines = f.readlines()
            for line in reversed(lines[-50:]):
                if "JOB DONE." in line:
                    return True
            return False

    def run_pw(self, input_file_path, npool=None):
        """Executes a pw.x calculation (SCF, NSCF, or BANDS) and validates the output."""
        # 1. Deduce the expected output file path
        expected_out_file = input_file_path.replace('/inputs', '/outputs').replace('.in', '.out')

        # 2. Execute the command
        out_file = self._execute_command("pw.x", input_file_path, npool=npool)

        # 3. Handle potential ungraceful MPI shutdowns (frequent in cluster environments)
        # If _execute_command returns None but the calculation finished, fallback to the expected file.
        file_to_check = out_file if out_file else expected_out_file

        # 4. Physical verification
        if self._check_job_done(file_to_check):
            if out_file is None:
                print(f"  -> [Warning] Ignored ungraceful MPI shutdown: physical calculation is intact in {file_to_check}")
            else:
                print(f"PW calculation successfully completed: {file_to_check}")

            return file_to_check
        else:
            print(f"PW calculation failed or incomplete: {file_to_check}")
            return None

    def run_dos(self, input_file_path):
        """Executes dos.x in serial mode and returns the paths to the output and .dos data files."""
        out_file = self._execute_command("dos.x", input_file_path, serial=True)
        if out_file and self._check_job_done(out_file):
            dos_data_file = os.path.join(self.output_folder, f"{self._get_prefix(input_file_path)}.dos")
            return {"out_file": out_file, "data_file": dos_data_file}

        return None

    def run_bands_pp(self, input_file_path):
        """Executes bands.x in serial mode and returns the paths to the generated data files."""
        out_file = self._execute_command("bands.x", input_file_path, serial=True)

        if out_file and self._check_job_done(out_file):
            prefix = self._get_prefix(input_file_path)

            gnu_data_file = os.path.join(self.output_folder, f"{prefix}_bands.dat.gnu")
            rap_data_file = os.path.join(self.output_folder, f"{prefix}_bands.dat.rap")

            return {
                "out_file": out_file,
                "gnu_file": gnu_data_file,
                "rap_file": rap_data_file
            }

        return None

    def run_fs(self, input_file_path):
        """Executes fs.x in serial mode and returns the path to the .bxsf file."""
        out_file = self._execute_command("fs.x", input_file_path, serial=True)
        if out_file and self._check_job_done(out_file):
            bxsf_data_file = os.path.join(self.output_folder, f"{self._get_prefix(input_file_path)}.bxsf")
            return {"out_file": out_file, "data_file": bxsf_data_file}

        return None

    def _get_prefix(self, input_path):
        """Helper method to extract the material prefix from the input file name."""
        base = os.path.basename(input_path)

        # Remove known suffixes to isolate the actual prefix (e.g., NiTe2_FR)
        known_suffixes = ['_dos.in', '_bands_pp.in', '_fs.in', '_bands.in', '_scf.in', '_nscf.in']
        for suffix in known_suffixes:
            if base.endswith(suffix):
                return base.replace(suffix, '')

        return base.split('_')[0]

    def run_batch(self, input_file_list):
        """
        Executes a sequence of pw.x calculations (e.g., for convergence tests).
        Returns a list of all successfully generated .out files.
        """
        successful_outputs = []
        print(f"Starting batch execution of {len(input_file_list)} calculations...")

        for input_file in input_file_list:
            out_file = self.run_pw(input_file)
            if out_file:
                successful_outputs.append(out_file)
            else:
                print(f"Batch interrupted: calculation on {input_file} failed.")
                break

        return successful_outputs