"""General helper functions."""

import os


def readlines_reverse(fname):
    """
    Read lines of a file in a reverse order.

    Parameters
    ----------
    fname : Path
        Path to the file to read.

    """
    with open(fname) as qfile:
        qfile.seek(0, os.SEEK_END)
        position = qfile.tell()
        line = ""
        while position >= 0:
            qfile.seek(position)
            next_char = qfile.read(1)
            if next_char == "\n":
                yield line[::-1]
                line = ""
            else:
                line += next_char
            position -= 1
        yield line[::-1]


def tail(filename, nlines):
    """
    Return the last n lines of a file.

    Parameters
    ----------
    fname : Path
        Path to the file to read.

    nlines : int
        Number of lines to return.

    """
    lines = []
    a = readlines_reverse(filename)
    while len(lines) < nlines:
        lines.append(next(a))

    return lines


def write_xyz_file(coords, fname, natoms=0, title=""):
    """
    Write xyz file.

    Parameters
    ----------
    coords : [atom x y z]
        A list of atom symbols/numbers and x, y, z coords in Anstrom.
    ofile : Path
        A path where the xyz file will be written.
    natoms : int
        Number of atoms (first line of the xyz file). If not specified, it will
        be inferred from the the length of the coords list.
    title : str
        Title of the xyz file (second line of the file).

    """
    if natoms == 0:
        natoms = len(coords)

    xyz = "\n".join((str(natoms), f"{title}", *coords))

    if not fname.parent.exists():
        fname.parent.mkdir(parents=True)

    with open(fname, "w", newline="\n") as f:
        f.write(xyz)
