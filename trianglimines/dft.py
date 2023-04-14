"""Functions to help with launching xTB and DFT codes."""

import re
import subprocess
from shutil import copyfile


def orca_input(
    workdir,
    coords_xyz,
    functional="b97-3c",
    basis="",
    opt=False,
    freq=False,
    convergence="",
    solvent=None,
    ncpus=20,
    memory=8000,
    slurm=False,
    jobname="orca_calc",
    time_lim=24,
    submit=False,
):
    """
    Write ORCA input and the corresponding SLURM script.

    Parameters
    ----------
    workdir : Path
        A path to where the calculation should be executed from.

    coords_xyz : Path
        A path to the input geometry.

    functional : str
        DFT functional to use (see Orca Manual for viable options).

    basis : str
        Basis set to use (see Orca Manual for viable options).

    opt : Bool
        If True then an optimisation is performed.

    freq : Bool
        If True then a frequency calculation is performed.

    convergence : str
        Convergence criteria used for geometry optimisations. If none are
        supplied then the defaults are applied. Other common options could be
        TIGHT or VERYTIGHT.

    solvent : str
        Solvent name for CPCM SMD solvation.

    ncpus : int
        Number of threads to use for a parallel calculation.

    memory : int
        Memory per core in parallel calculations (in MB).

    slurm : Bool
        If True then SLURM submission script is also written.

    jobname : str
        SLURM jobname.

    time_lim : int
        Time limit for the SLURM job (full hours only).

    submit : Bool
        If True then the job is submitted to SLURM scheduler straight away.

    """
    opt_flag = f"{convergence}OPT " if opt else ""

    if freq:
        freq_flag = "NUMFREQ " if solvent is not None else "FREQ "
    else:
        freq_flag = ""

    if solvent is not None:
        smd_block = "\n".join(
            ["%CPCM SMD TRUE", f'    SMDSOLVENT "{solvent}"', "END"]
        )
    else:
        smd_block = ""

    orca_input = "\n".join(
        [
            f"!{functional} {basis} {opt_flag}{freq_flag} LARGEPRINT",
            f"%MAXCORE {memory}",
            f"%PAL NPROCS {ncpus} END",
            smd_block,
            "* xyzfile 0 1 coords.xyz",
            "",
        ]
    )

    orca_submit = "\n".join(
        [
            "#!/bin/bash -l",
            "# Define job name",
            f"#SBATCH -J {jobname}",
            "# Use the current working directory",
            "#SBATCH -D ./",
            "# Use the current environment for this job.",
            "#SBATCH --export=ALL",
            "# Both standard output and error are directed to a file",
            "#SBATCH -o slurm_orca.out",
            "# Request the partition",
            "#SBATCH -p nodes,cooper",
            "# Request the number of nodes",
            "#SBATCH -N 1",
            "# Request the number of cores",
            f"#SBATCH -n {ncpus}",
            "# Specify the time limit",
            f"#SBATCH -t {time_lim}:00:00",
            "",
            "module purge",
            "module load apps/anaconda3",
            "module load load apps/orca",
            "",
            "$ORCADIR/orca orca_calc.inp > orca_calc_$SLURM_JOB_ID.out",
        ]
    )

    if not workdir.exists():
        workdir.mkdir(parents=True)

    with open(workdir / "orca_calc.inp", "w", newline="\n") as f:
        f.write(orca_input)

    copyfile(coords_xyz, workdir / "coords.xyz")

    if slurm:
        with open(workdir / "orca_submit.sh", "w", newline="\n") as f:
            f.write(orca_submit)

        if submit:
            subprocess.run(
                ["sbatch", "orca_submit.sh"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=workdir,
            )


def orca_get_energy(orcalog):
    """
    Get final single point energy from an Orca logfile.

    Parameters
    ----------
    orcalog : Path
        A path to the Orca log file.

    Returns
    -------
    float | None
        Final single point energy or None if the calculation did not converge.

    """
    logfile = orcalog.read_text()
    m = re.findall(r"FINAL SINGLE POINT ENERGY[ ]+(-\d*.\d*)", logfile)
    return float(m[-1]) if m != [] else None


def xtb_opt(
    infile,
    basename,
    logfile=None,
    gfn_version="2",
    solvent="chcl3",
    cycles=10000,
    ncpus=None,
):
    """
    Perform xtb optimisation.

    Parameters
    ----------
    infile : Path
        A path to the input molecular geometry (xyz format).

    basename : str
        A basename applied to the names all of created files (e.g., InChIKey).

    logfile : str
        A path to the logfile. If None, then default will be created and used.

    gfn_version : ["ff", 1, 2]
        GFN version to use in the optimisation.

    solvent : str
        Implicit solvent to be used in the calculation. ALPB model is used.
        For a list of valid solvents, check GFN-xTB documentation.

    cycles : int
        Maximum number of cycles.

    Returns
    -------
    None

    Notes
    -----
    Check GFN-xTB documentation for valid entries.

    """
    xtbopt = infile.parent / "xtbopt.xyz"

    if gfn_version == "ff":
        out_file = xtbopt.with_name(f"{basename}_GFN-FF.out")
        gfnflag = ["--gfnff"]
        xyz_file = xtbopt.with_stem(f"{basename}_GFN-FF")
        done = "GFN-FF optimisation done."

    else:
        out_file = xtbopt.with_name(f"{basename}_GFN{gfn_version}-xTB.out")
        gfnflag = ["--gfn", str(gfn_version)]
        xyz_file = xtbopt.with_stem(f"{basename}_GFN{gfn_version}-xTB")
        done = f"GFN{gfn_version}-xTB optimisation done."

    if logfile is not None:
        out_file = logfile

    solflag = ["--alpb", solvent] if solvent is not None else []
    cycflag = ["--cycles", str(cycles)] if cycles is not None else []
    pflag = ["--parallel", str(ncpus)] if ncpus is not None else []

    process = (
        ["xtb", infile.name, "--opt"] + gfnflag + solflag + cycflag + pflag
    )

    with open(out_file, "w", newline="\n") as f:
        x = subprocess.run(
            process,
            stdout=f,
            stderr=subprocess.PIPE,
            text=True,
            cwd=infile.parent,
        )

    print(x.stderr)
    xtbopt.rename(xyz_file)

    auxfiles = [
        xyz_file.parent / ".xtboptok",
        xyz_file.parent / "charges",
        xyz_file.parent / "gfnff_charges",
        xyz_file.parent / "gfnff_topo",
        xyz_file.parent / "wbo",
        xyz_file.parent / "xtbopt.log",
        xyz_file.parent / "xtbrestart",
        xyz_file.parent / "xtbtopo.mol",
    ]

    for file in auxfiles:
        file.unlink(missing_ok=True)

    print(done)


def xtb_crest(
    infile,
    odir,
    mtd_version="3",
    gfn_version="2",
    solvent="chcl3",
    rthr=0.5,
    ethr=0.05,
    quick=None,
    prop=None,
    cluster=None,
    ncpus=4,
):
    """
    Perform CREST conformer search with xtb.

    Parameters
    ----------
    infile : Path
        A path to the input molecular geometry (xyz format).

    odir : Path
        A directory in which CREST files will be saved.

    mtd_version : [1, 2, 3, 4, "entropy"]
        CREST searching algorithm version.

    gfn_version : ["ff", 1, 2]
        GFN version to use in the optimisation.

    solvent : str
        Implicit solvent to be used in the calculation. ALPB model is used.
        For a list of valid solvents, check GFN-xTB documentation.

    rthr : float
        RMSD threshold for conformer differentiation.

    ethr : float
        Energy threshold for conformer differentiation.

    quick : [None, "quick", "squick", "mquick"]
        Perform a search with reduced settings for a crude conformer ensemble.
        No reduced settings if None.

    prop : [None, "hess", "reopt", "autoIR"]
        Defines what should be done after the search. CREST manual recommends
        --prop reopt after --quick runs.

    cluster : int
        Number of clusters for PCA / k-means clustering. No clustering if None.

    ncpus : int
        Number of threads to use for CREST.

    Returns
    -------
    None

    Notes
    -----
    Check CREST / GFN-xTB documentation for valid entries.

    """
    if not odir.exists():
        odir.mkdir()

    coords = odir / "coords.xyz"
    copyfile(infile, coords)
    out_file = odir / "CREST_log.out"

    Tflag = ["-T", str(ncpus)] if ncpus is not None else []

    if mtd_version == "entropy":
        mtdflag = ["--entropy"]

    else:
        mtdflag = [f"--v{mtd_version}"]

    gfnflag = [f"--gfn{gfn_version}"]
    solflag = ["--alpb", solvent] if solvent is not None else []
    rthrflag = ["--rthr", str(rthr)] if rthr is not None else []
    ethrflag = ["--ethr", str(ethr)] if ethr is not None else []
    quickflag = [f"--{quick}"] if quick is not None else []
    propflag = ["--prop", prop] if prop is not None else []
    clusterflag = ["--cluster", cluster] if cluster is not None else []

    done = f"CREST algorithm {mtd_version} with GFN{gfn_version} done."

    process = (
        ["crest", coords.name]
        + Tflag
        + mtdflag
        + gfnflag
        + solflag
        + rthrflag
        + ethrflag
        + quickflag
        + propflag
        + clusterflag
    )

    with open(out_file, "w", newline="\n") as f:
        x = subprocess.run(
            process,
            stdout=f,
            stderr=subprocess.PIPE,
            text=True,
            cwd=odir,
        )

    print(x.stderr)
    print(done)
