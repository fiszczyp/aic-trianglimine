"""
Embed molecules for conformer search.

Requires
--------
rdkit >= 2021.09.5

"""
import json
import random

import rdkit.Chem.AllChem as rdchem

from params import input_mols_json, mols_json, ncpus, random_seed, structures
from trianglimines.molecule import Molecule, write_xyz_file

__author__ = "Filip T. Szczypiński"


def get_lowest_all_E_conf(opt_mol, confs_e):
    """
    Ensure all double bonds are E configuration.

    Currently RDKit does not constrain E/Z isomerism when embedding molecules
    so sometimes the imine bond ends up Z. That is never experimentally
    observed so for further calculations, only all E isomers are selected.
    This function filters through the output of `moleculeio.Molecule.embed()`
    and returns the lowest energy all E conformer.


    Parameters
    ----------
    opt_mol : rdkit.Molecule
        A molecule with embedded conformations.

    confs_e : [(int, float)]
        List of MMFF energies for the conformations.

    Returns
    -------
    int or None
        Conformer ID of the lowest energy all E conformation. If no such
        conformation is found, it returns None.

    """
    all_E = []

    for i, conf in enumerate(opt_mol.GetConformers()):
        rdchem.AssignStereochemistryFrom3D(
            opt_mol, i, replaceExistingTags=True
        )
        bonds = [
            x
            for x in conf.GetOwningMol().GetBonds()
            if x.GetBondType() == rdchem.BondType.DOUBLE
        ]
        if all([x.GetStereo() == rdchem.BondStereo.STEREOE for x in bonds]):
            all_E.append(i)

    confs_e = [(i, e[1]) for i, e in enumerate(confs_e) if i in all_E]

    return min(confs_e, key=lambda t: t[1])[0] if len(confs_e) > 0 else None


if __name__ == "__main__":
    # Random seed used for conformer generation
    random.seed(random_seed)

    with open(input_mols_json, "r") as f:
        molecules = [Molecule.from_dict(m) for m in json.load(f)]

    for m in molecules:
        odir = structures / m.inchikey
        if not odir.exists():
            odir.mkdir(parents=True)

        opt_mol, confs_e = m.embed(
            ncpus=ncpus,
            method=rdchem.ETKDGv3,
            nconf=500,
        )

        lowest_conf = get_lowest_all_E_conf(opt_mol, confs_e)
        if lowest_conf is None:
            # Try to get more and different conformers.
            # Otherwise, just embed lowest Z-containing conformer.
            print("No all E conformers found. Trying again.")

            random.seed(random_seed / 2)
            opt_mol, confs_e = m.embed(
                ncpus=ncpus,
                method=rdchem.ETKDGv3,
                nconf=1000,
            )
            random.seed(random_seed)
            lowest_conf = get_lowest_all_E_conf(opt_mol, confs_e)

        m.ETKDG_embedded = True

        if lowest_conf is not None:
            mol_xyz = rdchem.MolToXYZBlock(opt_mol, confId=lowest_conf)[2:]
            ofile_embedded = odir / f"{m.inchikey}_ETKDG_MMFF.xyz"
            write_xyz_file(
                coords=mol_xyz,
                fname=ofile_embedded,
            )
            print(f"{m.inchikey} successfully embedded.")

        else:
            lowest_conf = confs_e.index(min(confs_e, key=lambda t: t[1]))
            mol_xyz = rdchem.MolToXYZBlock(opt_mol, confId=lowest_conf)[2:]
            ofile_embedded = odir / f"{m.inchikey}_ETKDG_MMFF.xyz"
            write_xyz_file(
                coords=mol_xyz,
                fname=ofile_embedded,
            )
            print(f"{m.inchikey} embedded Z-containing conformation!")

    with open(mols_json, "w", newline="\n") as f:
        json.dump([m.to_dict() for m in molecules], f, indent=4)
