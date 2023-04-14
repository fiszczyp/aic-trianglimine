# Modelling files for Scholes *et al.*

The main objective is to calculate relative formation energies of
a family of bromo-substituted isotrianglimines to understand their
dynamic behaviour observed in solution.

## Packages used

* python >= 3.9.
* [rdkit][1] >= 2021.09.5, used 2022.03.2.
* [xtb][2] 6.4.1 (23d549d) and [CREST][3].
* [orca][4] 5.0.1 and Orca 5.0.3.
* [spyrmsd][5], works best with `graph-tool`.

[1]: https://doi.org/10.5281/zenodo.6483170
[2]: https://doi.org/10.1021/acs.jctc.8b01176
[3]: https://doi.org/10.1039/C9CP06869D
[4]: https://doi.org/10.1002/wcms.1606
[5]: https://doi.org/10.1186/s13321-020-00455-2

## Embedding molecules: `embed_molecules.py`

The workflow starts with embedding InChI's of the molecules using the
[ETKDGv3][6] algorithm as implemented in `rdkit`. The structures are minimised
with the MMFF force field and only all-*E* conformers are kept. If none
all-*E* structures are found, then further 1000 conformations are generated
with a different random seed. Lowest energy conformation thus identified is
used for further calculations. Requires RDkit minimum 2021.09.5 to reproduce
the isomer detection workflow. Constraints migt be fixed in a future release.

[6]: https://doi.org/10.1021/acs.jcim.0c00025

## Conformer search: `conformer_search.py`

With the embedded structures at hand, CREST is used with GFN2-xTB to perform
conformational search. Redundant conformers are then removed with an RMSD
thershold of 0.5Å using sPyRMSD. Symmetry-corrected RMSD calculations are
necessary for the symmetric [2+2] and [3+3] macrocycles.

Each conformer is optimised (using Orca) with the [B97-3c][7] functional and
confirmed to be a minimum with a frequency calculation.

[7]: https://doi.org/10.1063/1.5012601