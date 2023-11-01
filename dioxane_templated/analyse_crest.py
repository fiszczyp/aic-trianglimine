"""Analyse [4+4]⊃dioxane CREST search."""

import numpy as np
from pathlib import Path
from shutil import copyfile

# Atom IDs of some of the dioxane and macrocycle atoms:
CYCLE_IDS = [16, 17, 18, 19]
DIOXANE_IDS = [130, 133]
CONFS = Path("CREST_GFN2/confs")

if __name__ == "__main__":
    templated = []

    for conf in CONFS.glob("*.xyz"):
        cycle_xyz = conf.read_text().splitlines()[2:]
        cycle_coords = np.array(
            [x.split()[1:] for x in cycle_xyz][:-1], dtype=float
        )
        com = np.average(cycle_coords, axis=0)

        cycle_dists = [
            np.linalg.norm(cycle_coords[cycle_id] - com)
            for cycle_id in CYCLE_IDS
        ]

        dioxane_dists = [
            np.linalg.norm(cycle_coords[dioxane_id] - com)
            for dioxane_id in DIOXANE_IDS
        ]

        inside = True

        for dioxane in dioxane_dists:
            for cage in cycle_dists:
                if dioxane - cage > 0:
                    inside = False

        if inside:
            templated.append(conf)

    copyfile(templated[0], "crest_dioxane_templated.xyz")
