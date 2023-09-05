"""
Analyse single point energies to get formation energies.

Requires
--------
rdkit

"""
import json
import re

from params import mols_json, e_formation
from trianglimines.dft import to_kj
from trianglimines.molecule import Molecule

__author__ = "Filip T. Szczypiński"


def _generate_corrections(m):
    dH = (
        m.thermochemistry["Total Enthalpy (Eh)"]
        - m.thermochemistry["Electronic Energy (Eh)"]
    )
    dS_el_vib_trans = (
        m.thermochemistry["Electronic Entropy (Eh)"]
        + m.thermochemistry["Vibrational Entropy (Eh)"]
        + m.thermochemistry["Translational Entropy (Eh)"]
    )
    dS1 = dS_el_vib_trans + m.thermochemistry["Rotational Entropy, sn=1 (Eh)"]
    dS2 = dS_el_vib_trans + m.thermochemistry["Rotational Entropy, sn=2 (Eh)"]
    dS3 = dS_el_vib_trans + m.thermochemistry["Rotational Entropy, sn=3 (Eh)"]
    dS4 = dS_el_vib_trans + m.thermochemistry["Rotational Entropy, sn=4 (Eh)"]
    dS6 = dS_el_vib_trans + m.thermochemistry["Rotational Entropy, sn=6 (Eh)"]

    # G - E(el) = H - TS - E(el) = E(el) + dH - TS - E(el) = dH - TS
    m.corrections = {
        "dH": dH,  # H(corrected) - E(el)
        "TdS1": dS1,  # TdS for each symm number
        "TdS2": dS2,  # TdS for each symm number
        "TdS3": dS3,  # TdS for each symm number
        "TdS4": dS4,  # TdS for each symm number
        "TdS6": dS6,  # TdS for each symm number
        "dG1": dH - dS1,  # G - E(el) for each symm number
        "dG2": dH - dS2,  # G - E(el) for each symm number
        "dG3": dH - dS3,  # G - E(el) for each symm number
        "dG4": dH - dS4,  # G - E(el) for each symm number
        "dG6": dH - dS6,  # G - E(el) for each symm number
    }


if __name__ == "__main__":
    molecules = {}
    with open(mols_json, "r") as f:
        for m in json.load(f):
            mol = Molecule.from_dict(m)
            molecules[mol.inchikey] = mol

    for inchi, m in molecules.items():
        _generate_corrections(m)

    macrocycles = {}
    bbs = {
        "R-CHDA": "SSJXIUAHEKJCMH-PHDIDXHHSA-N",
        "S-CHDA": "SSJXIUAHEKJCMH-WDSKDSINSA-N",
        "ald": "QAQOLRCTTDVBAK-UHFFFAOYSA-N",
        "water": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
    }

    for inchi, m in molecules.items():
        if (n_imines := len(m.get_imines())) > 0:
            ald, am, chir = re.match(
                r"\[(\d+)\+(\d+)\]\-([R-S]*)", m.name
            ).groups()
            macrocycles[m.name] = {"inchikey": m.inchikey}

            for method, energy in m.DFT_energy.items():
                # Sum of electronic energies, E(el)
                am_energy = 0
                for c in chir:
                    if c == "R":
                        am_energy += molecules[bbs["R-CHDA"]].DFT_energy[
                            method
                        ]

                    elif c == "S":
                        am_energy += molecules[bbs["S-CHDA"]].DFT_energy[
                            method
                        ]
                e_form = (
                    energy
                    - am_energy
                    - int(ald) * molecules[bbs["ald"]].DFT_energy[method]
                    + n_imines * molecules[bbs["water"]].DFT_energy[method]
                )
                macrocycles[m.name][method] = {
                    "Formation energy (kJ/mol)": round(to_kj(e_form), 2),
                    "Formation energy per imine (kJ/mol)": round(
                        to_kj(e_form / n_imines), 2
                    ),
                }

                # Sum of thermally corrected enthalpies, E(el)+dH
                am_enthalpy = 0
                for c in chir:
                    if c == "R":
                        am_enthalpy += (
                            molecules[bbs["R-CHDA"]].DFT_energy[method]
                            + molecules[bbs["R-CHDA"]].corrections["dH"]
                        )

                    elif c == "S":
                        am_enthalpy += (
                            molecules[bbs["S-CHDA"]].DFT_energy[method]
                            + molecules[bbs["S-CHDA"]].corrections["dH"]
                        )
                h_form = (
                    energy
                    + m.corrections["dH"]
                    - am_enthalpy
                    - int(ald)
                    * (
                        molecules[bbs["ald"]].DFT_energy[method]
                        + molecules[bbs["ald"]].corrections["dH"]
                    )
                    + n_imines
                    * (
                        molecules[bbs["water"]].DFT_energy[method]
                        + molecules[bbs["water"]].corrections["dH"]
                    )
                )
                macrocycles[m.name][method].update(
                    {"Formation enthalpy (kJ/mol)": round(to_kj(h_form), 2)}
                )

                # Sum of thermally corrected entropies, TdS
                # sn = 2 for all aldehyde (C2v), amine (C2), and water (C2v)

                am_entropy = 0
                for c in chir:
                    if c == "R":
                        am_entropy += molecules[bbs["R-CHDA"]].corrections[
                            "TdS2"
                        ]

                    elif c == "S":
                        am_entropy += molecules[bbs["S-CHDA"]].corrections[
                            "TdS2"
                        ]

                s_form_partial = (
                    -am_entropy
                    - int(ald) * (molecules[bbs["ald"]].corrections["TdS2"])
                    + n_imines * (molecules[bbs["water"]].corrections["TdS2"])
                )

                for sn in [1, 2, 3, 4, 6]:
                    s_form = m.corrections[f"TdS{sn}"] + s_form_partial
                    g_form = h_form - s_form
                    macrocycles[m.name][method].update(
                        {
                            f"Formation TdS, sn={sn} (kJ/mol)": round(
                                to_kj(s_form), 2
                            ),
                            f"Formation Gibbs energy, sn={sn} (kJ/mol)": round(
                                to_kj(g_form), 2
                            ),
                            f"Formation Gibbs energy per imine, sn={sn}"
                            " (kJ/mol)": round(to_kj(g_form / n_imines), 2),
                        }
                    )

    with open(e_formation, "w", newline="\n") as f:
        json.dump(macrocycles, f, indent=4)
