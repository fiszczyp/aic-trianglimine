"""Functions to help with reading and writing molecules."""

import random
from itertools import combinations

import numpy as np
import rdkit.Chem.AllChem as rdchem
from spyrmsd import io, rmsd

from .helpers import write_xyz_file


class Molecule:
    """
    A class used to represent a molecule.

    Attributes
    ----------
    inchi : str
        An InChI string representing the molecule.

    inchikey : str
        An InChiKey corresponding to the molecule.

    name : str
        A common name referring to the molecule.

    """

    def __init__(self, inchi, name):
        """
        Initialise Molecule.

        Parameters
        ----------
        inchi : str
            An InChI string representing the molecule.

        name : str
            A common name referring to the molecule.

        """
        self.inchi = inchi
        self.name = name
        self.inchikey = rdchem.MolToInchiKey(rdchem.MolFromInchi(inchi))

    def __str__(self):
        """Return a description of the Molecule."""
        string = (
            f"The molecule {self.name} has InChI {self.inchi}, and the"
            f" corresponding InChIKey is {self.inchikey}"
        )
        return string

    def __repr__(self):
        """Return a representation of the Molecule."""
        string = (
            f"{self.__class__.__name__}(inchi='{self.inchi}',"
            f" name='{self.name}')"
        )
        return string

    @classmethod
    def from_dict(cls, mol_dict):
        """Initialise from a dictionary."""
        mol = cls(mol_dict["inchi"], mol_dict["name"])
        mol.__dict__.update(mol_dict)
        return mol

    def to_dict(self):
        """Convert the molecule to a dictionary."""
        return self.__dict__

    def embed(
        self,
        ncpus=4,
        method=rdchem.ETKDGv3,
        nconf=500,
    ):
        """
        Embed the molecule and identify the lowest energy conformation (crude).

        This function creates `nfconf` conformations using the ETKDGv3,
        minimises them using MMFF94 force field and identifies the lowest
        energy conformation. The conformers are minimised *in place*.

        Parameters
        ----------
        ncpus : int , optional
            Number of CPUs to use for minimisation (in parallel).

        method : rdkit.ETKDG
            RDkit ETKDG implementation, by default rdkit.ETKDGv3 but might
            consider an older version or rdkit.srETKDGv3.

        nconf : int
            Number of conformers to be embedded.

        Returns
        -------
        rdmol_h : rdkit.Molecule
            A molecule with embedded conformations.

        confs_e : [(int, float)]
            List of MMFF energies for the conformations.

        """
        params = method()
        params.randomSeed = random.randint(0, 100)

        rdmol_h = rdchem.AddHs(rdchem.MolFromInchi(self.inchi))
        rdchem.EmbedMultipleConfs(rdmol_h, nconf, params=params)
        confs_e = rdchem.MMFFOptimizeMoleculeConfs(
            rdmol_h,
            maxIters=10000,
            numThreads=ncpus,
        )

        return rdmol_h, confs_e


def remove_redundant_conformers(xyzfile, odir, basename, rthr=0.5):
    """Remove redundant conformers.

    This function removes redundant conformation from a multi-conformation
    xyz file. It uses symmetric-corrected RMSD function to correctly identify
    redundant conformations in highly symmetric macrocycles.

    Parameters
    ----------
    xyzfile : Path
        A path to the xyz file containing multiple conformations (e.g., the
        output of the CREST program: `crest_conformers.xyz`).

    odir : Path
        A path where the extracted non-redudant conformers will be saved as
        separate xyz files.

    basename : str
        A basename to be used for the created files (e.g., InChIKey).

    rthr : float
        The RMSD cut-off threshold for redundancy.

    Returns
    -------
    None

    """
    confs = io.loadallmols(str(xyzfile))
    all_confs = len(confs)

    redundant_confs = set()

    for ci, cj in combinations(range(all_confs), 2):
        conf_i = confs[ci]
        conf_j = confs[cj]
        rmsd_ij = rmsd.symmrmsd(
            conf_i.coordinates,
            conf_j.coordinates,
            conf_i.atomicnums,
            conf_j.atomicnums,
            conf_i.adjacency_matrix,
            conf_j.adjacency_matrix,
            center=True,
            minimize=True,
        )
        if rmsd_ij < rthr:
            print(
                f"Conformers {ci} and {cj} are REDUNDANT"
                f" ({np.around(rmsd_ij, decimals=3)}).\n"
                f"    => Removing {cj}."
            )
            redundant_confs.add(cj)

    unique_confs = sorted(set(range(all_confs)) - redundant_confs)

    print(f"\nUnique conformers are {unique_confs}.")

    for i, conf_id in enumerate(unique_confs):
        conf = confs[conf_id]
        coords = [
            f"{natom:<3} {coords[0]:>20.10f} {coords[1]:>20.10f}"
            f"{coords[2]:>20.10f}"
            for natom, coords in zip(conf.atomicnums, conf.coordinates)
        ]

        write_xyz_file(
            coords=coords,
            fname=(odir / f"{basename}_{i}.xyz"),
            natoms=len(coords),
            title=f"{basename}",
        )
