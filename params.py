"""Set up paths and parameters for the project."""

from pathlib import Path

random_seed = 1644321874
ncpus = 20
project = Path.cwd()

data = project / "data"
input_mols_json = data / "input_molecules.json"
mols_json = data / "molecules.json"
structures = data / "conformers"
dft_files = data / "dft_files"
opt_struct = data / "opt_struct"
