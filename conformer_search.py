"""
Perform CREST search and remove redudant conformers.

Requires
--------
orca
xtb
spyrmsd

"""
import json
from shutil import copyfile
from typing import Type

from params import dft_files, mols_json, ncpus, opt_struct, structures
from trianglimines.dft import (
    orca_check_imag_freq,
    orca_extract_thermochemistry,
    orca_get_energy,
    orca_input,
    xtb_crest,
    xtb_opt,
)
from trianglimines.molecule import Molecule, remove_redundant_conformers

__author__ = "Filip T. Szczypiński"


def _perform_crest(m):
    odir = structures / m.inchikey
    embedded = odir / f"{m.inchikey}_ETKDG_MMFF.xyz"

    xtb_opt(
        embedded,
        basename=m.inchikey,
        gfn_version="2",
        solvent="chcl3",
        cycles=10000,
    )

    crest_odir = odir / "CREST_GFN2"

    xtb_crest(
        odir / f"{m.inchikey}_GFN2-xTB.xyz",
        crest_odir,
        gfn_version=2,
        solvent="chcl3",
        rthr=0.5,
        ncpus=ncpus,
    )

    remove_redundant_conformers(
        xyzfile=(crest_odir / "crest_conformers.xyz"),
        odir=(crest_odir / "confs"),
        basename=m.inchikey,
    )

    m.CREST_completed = True


def _optimise_conformers_dft(m: Type[Molecule]):
    confs_path = structures / m.inchikey / "CREST_GFN2" / "confs"
    opt_confs_path = dft_files / "b97-3c-chcl3" / "opt" / m.inchikey

    for conf_log in confs_path.glob("*.xyz"):
        orca_input(
            workdir=(opt_confs_path / conf_log.stem),
            coords_xyz=conf_log,
            functional="B97-3c",
            opt=True,
            freq=False,
            solvent="CHLOROFORM",
            ncpus=ncpus,
            memory=3000,
        )

    m.confs_optimised = True


def _rank_conformers_dft(m: Type[Molecule]):
    opt_confs_path = dft_files / "b97-3c-chcl3" / "opt" / m.inchikey

    lowest_e = 0
    lowest_conf = 0

    for conf_log in opt_confs_path.glob("*/orca_calc.out"):
        if (e := orca_get_energy(conf_log)) is None:
            energy = 0
        else:
            energy = e

        if energy < lowest_e:
            lowest_e = energy
            lowest_conf = conf_log.parent.name.split("_")[1]

    m.lowest_conf = lowest_conf
    m.DFT_energy = {"B97-3c (Eh)": lowest_e}

    opt_struct.mkdir(parents=True, exist_ok=True)

    copyfile(
        (opt_confs_path / f"{m.inchikey}_{lowest_conf}" / "orca_calc.xyz"),
        opt_struct / f"{m.inchikey}.xyz",
    )


def _frequency_dft(m: Type[Molecule]):
    lowest_conf = f"{m.inchikey}_{m.lowest_conf}"
    opt_confs_path = dft_files / "b97-3c-chcl3" / "opt" / m.inchikey
    freq_path = dft_files / "b97-3c-chcl3" / "freq" / m.inchikey

    orca_input(
        workdir=(freq_path),
        coords_xyz=(opt_confs_path / lowest_conf / "orca_calc.xyz"),
        functional="B97-3c",
        opt=False,
        freq=True,
        solvent="CHLOROFORM",
        ncpus=ncpus,
        memory=8000,
    )

    m.freq_completed = True


def _frequency_analysis(m: Type[Molecule]):
    freq_output_path = dft_files / "b97-3c-chcl3" / "freq" / m.inchikey
    freq_output = [*freq_output_path.glob("orca_calc_*.out")][0]

    m.imaginary_freq = orca_check_imag_freq(freq_output)

    if not m.imaginary_freq:
        m.thermochemistry = orca_extract_thermochemistry(freq_output)


if __name__ == "__main__":
    with open(mols_json, "r") as f:
        molecules = [Molecule.from_dict(m) for m in json.load(f)]

    for m in molecules:
        if not hasattr(m, "CREST_completed") or not m.CREST_completed:
            _perform_crest(m)

    for m in molecules:
        if not hasattr(m, "confs_optimised") or not m.confs_optimised:
            _optimise_conformers_dft(m)

    for m in molecules:
        if not hasattr(m, "lowest_conf"):
            _rank_conformers_dft(m)

    for m in molecules:
        if not hasattr(m, "freq_completed") or not m.freq_completed:
            _frequency_dft(m)

    for m in molecules:
        if not hasattr(m, "imaginary_freq") or m.imaginary_freq:
            _frequency_analysis(m)

    with open(mols_json, "w", newline="\n") as f:
        json.dump([m.to_dict() for m in molecules], f, indent=4)
