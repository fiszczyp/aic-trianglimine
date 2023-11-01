from pathlib import Path
import argparse
import rdkit.Chem.AllChem as rdchem


def extract_confs(infile, odir):
    """
    Extract conformer xyz's from a multiple conformer file.

    Parameters
    ----------
    infile : Path
        A path to the xyz files containing multiple geometries.

    odir : Path
        A path to the directory where xyz files are to be extracted.

    """
    structs = []
    odir = Path(odir)
    infile = Path(infile)

    with open(infile, "r") as f:
        while (x := f.readline()) != "":
            natoms = int(x.split()[0])
            label = f.readline().split("\n")[0].split()[0]
            struct = [natoms, label]
            for i in range(natoms):
                struct.append(f.readline().split("\n")[0])
            structs.append(struct)

    if not odir.exists():
        odir.mkdir(parents=True)

    for i, struct in enumerate(structs):
        ofile = odir / f"conf_{i}.xyz"
        _write_xyz_file(
            natoms=struct[0],
            title=struct[1],
            coords=struct[2:],
            ofile=ofile,
        )


def _write_xyz_file(natoms, title, coords, ofile):
    xyz = "\n".join((str(natoms), f"{title}", *coords))

    with open(ofile, "w", newline="\n") as f:
        f.write(xyz)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split multi-structure xyz file."
    )
    parser.add_argument("infile", help="Multi-structure xyz file.")
    parser.add_argument("outdir", help="Output directory.")

    args = parser.parse_args()

    extract_confs(args.infile, args.outdir)
