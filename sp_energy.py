"""
Perform single point energy calculations at different theory levels.

Requires
--------
orca

"""
import json
from itertools import product
from typing import Optional

from params import (
    analysis_only,
    dft_files,
    mols_json,
    ncpus,
    opt_struct,
    sp_basis_sets,
    sp_functionals,
)

from trianglimines.dft import (
    orca_check_termination,
    orca_get_energy,
    orca_input,
)
from trianglimines.molecule import Molecule

__author__ = "Filip T. Szczypiński"


def _sp_dft(
    m: Molecule,
    funct: str,
    bset: str,
    solvent: Optional[str],
) -> None:
    workdir = (
        (
            dft_files
            / bset
            / f"{funct.replace(' ', '-')}-{solvent}"
            / m.inchikey
        )
        if solvent is not None
        else (dft_files / bset / funct.replace(" ", "-") / m.inchikey)
    )
    method = (
        f"{funct.replace(' ', '-')}-{solvent}/{bset}"
        if solvent is not None
        else f"{funct.replace(' ', '-')}/{bset}"
    )
    structure = opt_struct / f"{m.inchikey}.xyz"

    orca_input(
        workdir=workdir,
        coords_xyz=structure,
        functional=funct,
        opt=False,
        freq=False,
        solvent=solvent,
        ncpus=ncpus,
        memory=8000,
        slurm=True,
        submit=(not analysis_only),
    )

    m.DFT_energy = {method: "submitted"}


def _sp_analysis(
    m: Molecule,
    funct: str,
    bset: str,
    solvent: Optional[str],
) -> None:
    workdir = (
        (
            dft_files
            / bset
            / f"{funct.replace(' ', '-')}-{solvent}"
            / m.inchikey
        )
        if solvent is not None
        else (dft_files / bset / funct.replace(" ", "-") / m.inchikey)
    )
    method = (
        f"{funct.replace(' ', '-')}-{solvent}/{bset}"
        if solvent is not None
        else f"{funct.replace(' ', '-')}/{bset}"
    )
    orca_logfile = [*workdir.glob("orca_calc_*.out")][-1]
    m.sp_completed = orca_check_termination(orca_logfile)

    if m.sp_completed:
        m.DFT_energy[method] = orca_get_energy(orca_logfile)

    else:
        m.DFT_energy[method] = "timed_out"


if __name__ == "__main__":
    with open(mols_json, "r") as f:
        molecules = [Molecule.from_dict(m) for m in json.load(f)]

    for m in molecules:
        for funct, bset in product(sp_functionals, sp_basis_sets):
            # SMD = CHCl3
            solvent = "CHLOROFORM"
            method = f"{funct.replace(' ', '-')}-{solvent}/{bset}"
            if (
                method not in m.DFT_energy
                or m.DFT_energy[method] == "timed_out"
            ):
                _sp_dft(m, funct, bset, solvent=solvent)

            # Vacuum
            solvent = None
            method = f"{funct.replace(' ', '-')}/{bset}"
            if (
                method not in m.DFT_energy
                or m.DFT_energy[method] == "timed_out"
            ):
                _sp_dft(m, funct, bset, solvent=solvent)

    for m in molecules:
        for funct, bset in product(sp_functionals, sp_basis_sets):
            # SMD = CHCl3
            solvent = "CHLOROFORM"
            method = f"{funct.replace(' ', '-')}-{solvent}/{bset}"
            if m.DFT_energy[method] == "submitted" or analysis_only:
                _sp_analysis(m, funct, bset, solvent=solvent)

            # Vacuum
            solvent = None
            method = f"{funct.replace(' ', '-')}/{bset}"
            if m.DFT_energy[method] == "submitted" or analysis_only:
                _sp_analysis(m, funct, bset, solvent=solvent)

    with open(mols_json, "w", newline="\n") as f:
        json.dump([m.to_dict() for m in molecules], f, indent=4)
