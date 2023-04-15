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

sp_basis_sets = ["def2-QZVP"]
sp_functionals = [
    "PBE0 D3BJ",
    "PW6B95 D3BJ",
    "M062X D3ZERO",
    "wB97X-D3BJ",
    "wB97X-V",
    "wB97M-V",
]
