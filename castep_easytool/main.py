# -*- coding: utf-8 -*-
"""
A command-line tool to simplify common pre- and post-processing tasks
for the CASTEP plane-wave DFT code.
"""
from pathlib import Path
from typing import List, Optional

import typer
import re
import seekpath
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

app = typer.Typer(
    name="castep-easytool",
    help="⚛️ **castep-easytool**: A CLI to simplify CASTEP workflows.",
    add_completion=False,
    rich_markup_mode="markdown"
)

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _read_structure(filepath: Path) -> Structure:
    if not filepath.exists():
        typer.secho(f"❌ Error: File not found at '{filepath}'", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    try:
        struct = Structure.from_file(filepath)
    except Exception as e:
        typer.secho(f"❌ Error reading structure from '{filepath}': {e}", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    if not struct.is_ordered:
        typer.secho(f"❌ '{filepath}' contains disordered (partially occupied) sites, "
                    "which CASTEP .cell files cannot represent. Please supply an ordered structure.",
                    fg=typer.colors.RED)
        raise typer.Exit(code=1)
    return struct


def _lattice_abc_block(struct: Structure) -> str:
    a, b, c = struct.lattice.abc
    al, be, ga = struct.lattice.angles
    return (
        "%BLOCK LATTICE_ABC\n"
        f"  {a:.10f}  {b:.10f}  {c:.10f}\n"
        f"  {al:.10f}  {be:.10f}  {ga:.10f}\n"
        "%ENDBLOCK LATTICE_ABC\n"
    )


def _lattice_cart_block(struct: Structure) -> str:
    m = struct.lattice.matrix
    return (
        "%BLOCK LATTICE_CART\n"
        f"  {m[0][0]:.10f}  {m[0][1]:.10f}  {m[0][2]:.10f}\n"
        f"  {m[1][0]:.10f}  {m[1][1]:.10f}  {m[1][2]:.10f}\n"
        f"  {m[2][0]:.10f}  {m[2][1]:.10f}  {m[2][2]:.10f}\n"
        "%ENDBLOCK LATTICE_CART\n"
    )


def _positions_frac_block(struct: Structure) -> str:
    lines = ["%BLOCK POSITIONS_FRAC"]
    for site in struct:
        lines.append(f"  {site.specie.symbol:<4} {site.a:.10f}  {site.b:.10f}  {site.c:.10f}")
    lines.append("%ENDBLOCK POSITIONS_FRAC")
    return "\n".join(lines) + "\n"


def _write_param(cell_path: Path, header: str, keywords: dict) -> None:
    """
    Write the .param file that accompanies ``cell_path`` (same stem).

    Calculation-type keywords such as TASK and SPECTRAL_TASK are .param keywords;
    CASTEP rejects them in a .cell file. An existing .param file is never
    overwritten.
    """
    param_path = cell_path.with_suffix(".param")
    if param_path.exists():
        typer.secho(f"⚠️  '{param_path}' already exists and was left unchanged.", fg=typer.colors.YELLOW)
        return
    width = max(len(k) for k in keywords)
    body = "".join(f"{k:<{width}} : {v}\n" for k, v in keywords.items())
    param_path.write_text(f"# {header}\n\n{body}", encoding="utf-8")
    typer.secho(f"✅ Wrote CASTEP parameters to '{param_path}'", fg=typer.colors.GREEN)


# ---------------------------------------------------------------------------
# Command 1: convert
# ---------------------------------------------------------------------------

@app.command()
def convert(
    filepath: Path = typer.Argument(..., help="Path to the input .cif file."),
    to: str = typer.Option("primitive", "--to", help="Target cell type: 'primitive' or 'conventional'.")
):
    """
    Convert a .cif between **conventional** and **primitive** standard cells.
    """
    if to.lower() not in ("primitive", "conventional"):
        typer.secho(f"❌ Invalid --to value '{to}': use 'primitive' or 'conventional'.", fg=typer.colors.RED)
        raise typer.Exit(code=1)

    typer.echo(f"📖 Reading structure from '{filepath}'...")
    struct = _read_structure(filepath)
    sga = SpacegroupAnalyzer(struct)

    if to.lower() == "primitive":
        output_struct = sga.get_primitive_standard_structure()
        suffix = "_primitive"
    else:
        output_struct = sga.get_conventional_standard_structure()
        suffix = "_conventional"

    output_filename = filepath.parent / f"{filepath.stem}{suffix}.cif"
    CifWriter(output_struct).write_file(output_filename)
    typer.secho(f"✅ Wrote {to} cell to '{output_filename}'", fg=typer.colors.GREEN)


# ---------------------------------------------------------------------------
# Command 2: generate
# ---------------------------------------------------------------------------

@app.command()
def generate(
    cif_file: Path = typer.Argument(..., help="Path to the input .cif file."),
    task: str = typer.Option("GeometryOptimization", "--task", "-t", help="CASTEP task keyword."),
    cutoff: float = typer.Option(None, "--cutoff", "-e", help="Plane-wave cut-off energy in eV (CUT_OFF_ENERGY)."),
    xc: str = typer.Option(None, "--xc", help="XC functional, e.g. PBE, PBESOL, LDA, HSE06 (XC_FUNCTIONAL)."),
    kspacing: float = typer.Option(None, "--kspacing", "-k", help="k-point spacing in 1/A (KPOINTS_MP_SPACING)."),
    spin: bool = typer.Option(False, "--spin", help="Spin-polarised calculation (SPIN_POLARIZED : true)."),
    fix_cell: bool = typer.Option(False, "--fix-cell", help="Keep the cell fixed, optimise only atomic positions (FIX_ALL_CELL : true)."),
):
    """
    Generate a CASTEP **.cell** input file (and a matching **.param** file
    holding the task) from a .cif file.

    Optional flags add common settings:

        castep-easytool generate input.cif --cutoff 600 --xc PBESOL --kspacing 0.05

        castep-easytool generate input.cif --task SinglePoint --spin --fix-cell
    """
    typer.echo(f"📖 Reading structure from '{cif_file}'...")
    struct = _read_structure(cif_file)

    typer.echo(f"📝 Building .cell file for task '{task}'...")
    extras = ""
    if kspacing is not None:
        extras += f"\nKPOINTS_MP_SPACING : {kspacing:g}\n"
    if fix_cell:
        extras += "FIX_ALL_CELL : true\n"
    cell_content = (
        f"# CASTEP .cell file generated by castep-easytool\n\n"
        f"{_lattice_abc_block(struct)}\n"
        f"{_positions_frac_block(struct)}"
        f"{extras}"
    )

    output_filename = cif_file.parent / f"{cif_file.stem}.cell"
    output_filename.write_text(cell_content, encoding="utf-8")
    typer.secho(f"✅ Wrote CASTEP input to '{output_filename}'", fg=typer.colors.GREEN)
    params = {"TASK": task}
    if cutoff is not None:
        params["CUT_OFF_ENERGY"] = f"{cutoff:g} eV"
    if xc:
        params["XC_FUNCTIONAL"] = xc
    if spin:
        params["SPIN_POLARIZED"] = "true"
    _write_param(output_filename, "CASTEP .param file generated by castep-easytool", params)


# ---------------------------------------------------------------------------
# Command 3: bands
# ---------------------------------------------------------------------------

@app.command()
def bands(
    cif_file: Path = typer.Argument(..., help="Path to the input .cif file."),
    path: str = typer.Option(None, "--path", help="Custom k-path as high-symmetry labels, e.g. \"W-L-GAMMA-X-W-K\" "
                             "(G = GAMMA; use '|' for a jump between disconnected pieces). "
                             "Default: the SeeK-path standard path."),
):
    """
    Generate a CASTEP **bands.cell** file with a Seekpath BZ k-path, plus a
    matching **bands.param** file (TASK : SPECTRAL, SPECTRAL_TASK : BANDSTRUCTURE).

    The structure is standardised to the seekpath primitive cell automatically.
    Use `--path` to choose the points yourself, e.g. `--path W-L-GAMMA-X-W-K`.
    """
    typer.echo(f"📖 Reading structure from '{cif_file}'...")
    struct = _read_structure(cif_file)

    typer.echo("⚡️ Analysing structure and generating k-path with Seekpath...")
    try:
        path_data = seekpath.get_path(
            (struct.lattice.matrix, struct.frac_coords, struct.atomic_numbers),
            with_time_reversal=True
        )
    except Exception as e:
        typer.secho(f"❌ Seekpath could not determine a k-path for '{cif_file}': {e}", fg=typer.colors.RED)
        raise typer.Exit(code=1)

    # Primitive structure from seekpath
    prim = Structure(
        lattice=path_data['primitive_lattice'],
        species=path_data['primitive_types'],
        coords=path_data['primitive_positions'],
        to_unit_cell=True,
    )

    # Build SPECTRAL_KPOINT_PATH block.
    # Seekpath returns a list of (start, end) segments. Consecutive segments are
    # not always connected (e.g. fcc: ('X','U') is followed by ('K','GAMMA')), so
    # the start of a segment must be written whenever it differs from the end of
    # the previous one. CASTEP always joins consecutive vertices, so such a jump
    # is flagged with a "! BREAK" comment: castep-easytool bandplot uses it to
    # hide the spurious connecting segment (CASTEP itself ignores the comment).
    coords = path_data['point_coords']
    segments = path_data['path']
    if path:
        names = {k.upper(): k for k in coords}
        names.update({"G": "GAMMA", "Γ": "GAMMA"})
        segments = []
        for piece in path.split("|"):
            labels = []
            for tok in re.split(r"[-,\s]+", piece.strip()):
                if not tok:
                    continue
                key = names.get(tok.upper(), names.get(tok))
                if key is None or key not in coords:
                    typer.secho(f"❌ Unknown k-point '{tok}'. Available for this lattice: "
                                + ", ".join(coords), fg=typer.colors.RED)
                    raise typer.Exit(code=1)
                labels.append(key)
            segments += list(zip(labels[:-1], labels[1:]))
        if not segments:
            typer.secho("❌ --path needs at least two k-points.", fg=typer.colors.RED)
            raise typer.Exit(code=1)
        typer.echo("  Using custom path: " + " - ".join([segments[0][0]] + [e for _, e in segments]))
    kpath_lines = []

    def _vertex(label):
        c = coords[label]
        kpath_lines.append(f"  {c[0]:.8f}  {c[1]:.8f}  {c[2]:.8f}  ! {label}")

    for i, (start, end) in enumerate(segments):
        if i == 0:
            _vertex(start)
        elif start != segments[i - 1][1]:
            kpath_lines.append("  ! BREAK (discontinuity in the path)")
            _vertex(start)
        _vertex(end)
    kpath_block = "%BLOCK SPECTRAL_KPOINT_PATH\n" + "\n".join(kpath_lines) + "\n%ENDBLOCK SPECTRAL_KPOINT_PATH\n"

    cell_content = (
        f"# CASTEP bands.cell file generated by castep-easytool (Seekpath)\n\n"
        f"{_lattice_abc_block(prim)}\n"
        f"{_positions_frac_block(prim)}\n"
        f"{kpath_block}"
    )

    output_filename = cif_file.parent / f"{cif_file.stem}_bands.cell"
    output_filename.write_text(cell_content, encoding="utf-8")
    typer.secho(f"✅ Wrote CASTEP bands input to '{output_filename}'", fg=typer.colors.GREEN)
    _write_param(output_filename, "CASTEP bands .param file generated by castep-easytool",
                 {"TASK": "SPECTRAL", "SPECTRAL_TASK": "BANDSTRUCTURE"})


# ---------------------------------------------------------------------------
# Command 4: geomopt
# ---------------------------------------------------------------------------

@app.command()
def geomopt(
    castep_file: Path = typer.Argument(..., help="Path to the CASTEP geometry-optimisation output (.castep) file."),
    prefix: str = typer.Option(None, "--prefix", "-p", help="Output filename stem (default: input stem).")
):
    """
    Extract the final optimised geometry from a CASTEP geometry-optimisation
    output file and write a **.cell** file ready for the next calculation.

    Supports LBFGS, BFGS, DELOC, TPSD, and FIRE optimisers.
    """
    import re
    import numpy as np
    from pymatgen.core import Lattice, Structure

    if not castep_file.exists():
        typer.secho(f"❌ Error: File not found at '{castep_file}'", fg=typer.colors.RED)
        raise typer.Exit(code=1)

    typer.echo(f"📖 Parsing '{castep_file}'...")
    text = castep_file.read_text(encoding='utf-8', errors='replace')

    # Locate the last "Final Configuration" block for atomic positions
    config_re = re.compile(
        r'(?:LBFGS|BFGS|FIRE|DELOC|TPSD|GeometryOptimization): Final Configuration:',
        re.IGNORECASE
    )
    matches = list(config_re.finditer(text))
    if matches:
        tail = text[matches[-1].start():]
    else:
        unit_cell_blocks = list(re.finditer(r'^\s*Unit Cell\s*$', text, re.MULTILINE))
        if not unit_cell_blocks:
            typer.secho("❌ Cannot find any geometry block in the .castep file.", fg=typer.colors.RED)
            raise typer.Exit(code=1)
        tail = text[unit_cell_blocks[-1].start():]

    # Parse Cartesian lattice vectors from the last Unit Cell block in the file.
    # CASTEP only prints the Unit Cell block once per geometry step (not inside
    # the "Final Configuration" summary), so we search the whole file and take
    # the last occurrence.
    latt_re = re.compile(
        r'Real Lattice\(A\)[^\n]*\n'
        r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)[^\n]*\n'
        r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)[^\n]*\n'
        r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)'
    )
    latt_matches = list(latt_re.finditer(text))
    if not latt_matches:
        typer.secho("❌ Cannot parse lattice vectors from the final configuration block.", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    m = latt_matches[-1]
    latt_matrix = np.array([float(x) for x in m.groups()]).reshape(3, 3)

    # Parse fractional coordinates
    atom_re = re.compile(
        r'x\s+([A-Za-z]+)\s+\d+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+x',
        re.MULTILINE
    )
    species, frac_coords = [], []
    for am in atom_re.finditer(tail):
        species.append(am.group(1))
        frac_coords.append([float(am.group(2)), float(am.group(3)), float(am.group(4))])

    if not species:
        # Newer CASTEP versions print the positions only for the starting structure and write
        # every step to <seed>.geom, so read the final geometry from the last frame there.
        geom = castep_file.with_suffix(".geom")
        if geom.exists():
            bohr = 0.529177210903
            frame = geom.read_text(encoding="utf-8", errors="replace").split("<-- c")[-1]
            h_rows, r_rows = [], []
            for line in frame.splitlines():
                parts = line.split()
                if line.rstrip().endswith("<-- h") and len(parts) >= 4:
                    h_rows.append([float(x.replace("D", "E")) * bohr for x in parts[:3]])
                elif line.rstrip().endswith("<-- R") and len(parts) >= 6:
                    r_rows.append((parts[0], [float(x.replace("D", "E")) * bohr for x in parts[2:5]]))
            if len(h_rows) == 3 and r_rows:
                latt_matrix = np.array(h_rows)
                species = [sp for sp, _ in r_rows]
                cart = np.array([xyz for _, xyz in r_rows])
                frac_coords = (cart @ np.linalg.inv(latt_matrix)).tolist()
                typer.echo(f"  Using the final frame of '{geom.name}'.")

    if not species:
        typer.secho("❌ Cannot parse atomic positions from the final configuration block "
                    "(and no usable .geom file was found next to it).", fg=typer.colors.RED)
        raise typer.Exit(code=1)

    # Check convergence
    converged = bool(re.search(r'Geometry optimization completed successfully', text, re.IGNORECASE))
    iter_matches = re.findall(r'(?:LBFGS|BFGS|FIRE|DELOC|TPSD).*?:\s*finished iteration\s+(\d+)', text, re.IGNORECASE)
    n_iter = int(iter_matches[-1]) if iter_matches else 0

    typer.echo(f"  Species : {species}")
    typer.echo(f"  N atoms : {len(species)}")
    typer.echo(f"  Iterations: {n_iter}  |  Converged: {converged}")
    if not converged:
        typer.secho("  ⚠️  Geometry optimisation did NOT converge — using last available geometry.", fg=typer.colors.YELLOW)

    struct = Structure(
        lattice=Lattice(latt_matrix),
        species=species,
        coords=frac_coords,
        coords_are_cartesian=False,
    )

    out_stem = prefix or castep_file.stem
    output_filename = castep_file.parent / f"{out_stem}_opted.cell"

    cell_content = (
        f"# CASTEP .cell file extracted by castep-easytool from '{castep_file.name}'\n\n"
        f"{_lattice_abc_block(struct)}\n"
        f"{_positions_frac_block(struct)}\n"
    )
    output_filename.write_text(cell_content, encoding="utf-8")
    typer.secho(f"✅ Written: '{output_filename}'", fg=typer.colors.GREEN)

    cif_filename = castep_file.parent / f"{out_stem}_opted.cif"
    CifWriter(struct).write_file(str(cif_filename))
    typer.secho(f"✅ Written: '{cif_filename}'", fg=typer.colors.GREEN)


# ---------------------------------------------------------------------------
# Command 5: supercell
# ---------------------------------------------------------------------------

@app.command()
def supercell(
    cell_file: Path = typer.Argument(..., help="Path to the input structure file (.cif or any other format pymatgen can read; CASTEP .cell files are not supported)."),
    na: int = typer.Argument(..., min=1, help="Scale factor along a-axis."),
    nb: int = typer.Argument(None, min=1, help="Scale factor along b-axis (defaults to na)."),
    nc: int = typer.Argument(None, min=1, help="Scale factor along c-axis (defaults to na).")
):
    """
    Generate a CASTEP supercell **.cell** file from a structure file (.cif).

    Examples:

      castep-easytool supercell input.cif 2

      castep-easytool supercell input.cif 2 2 1
    """
    if nb is None:
        nb = na
    if nc is None:
        nc = na
    scaling_matrix = [na, nb, nc]
    suffix = f"_{na}x{nb}x{nc}"

    typer.echo(f"📖 Reading structure from '{cell_file}'...")
    struct = _read_structure(cell_file)

    typer.echo(f"🔲 Building {suffix.lstrip('_')} supercell...")
    superstruct = struct.make_supercell(scaling_matrix)

    cell_content = (
        f"# CASTEP .cell file generated by castep-easytool\n\n"
        f"{_lattice_cart_block(superstruct)}\n"
        f"{_positions_frac_block(superstruct)}"
    )

    output_filename = cell_file.parent / f"{cell_file.stem}{suffix}.cell"
    output_filename.write_text(cell_content, encoding="utf-8")
    typer.secho(f"✅ Wrote supercell to '{output_filename}'", fg=typer.colors.GREEN)
    _write_param(output_filename, "CASTEP .param file generated by castep-easytool",
                 {"TASK": "GeometryOptimization"})


# ---------------------------------------------------------------------------
# Command 6: bandplot
# ---------------------------------------------------------------------------

@app.command()
def bandplot(
    bands_file: Path = typer.Argument(..., help="Path to the CASTEP .bands file."),
    cell_file: Path = typer.Option(None, "--cell", "-c",
                                   help="Path to the .cell file (default: same stem as .bands file)."),
    emin: float = typer.Option(-6.0, "--emin", help="Plot energy minimum relative to VBM (eV)."),
    emax: float = typer.Option(6.0,  "--emax", help="Plot energy maximum relative to VBM (eV)."),
    prefix: str = typer.Option(None, "--prefix", "-p",
                               help="Output filename stem (default: .bands file stem)."),
):
    """
    Plot a CASTEP band structure and write a band gap summary.

    Reads a **.bands** file produced by a CASTEP spectral/bands calculation
    and the matching **.cell** file (auto-detected by default).

    **Outputs** (written to the same folder as the .bands file):

    - `<stem>_bands.pdf` — vector plot for publications
    - `<stem>_bands.png` — raster plot for quick viewing
    - `<stem>_bandgap.txt` — band gap, direct/indirect label, transition k-points;
      spin-resolved (spin-up / spin-down / overall) for spin-polarised calculations

    **Examples:**

        castep-easytool bandplot band.bands

        castep-easytool bandplot band.bands --emin -3 --emax 3

        castep-easytool bandplot band.bands --cell Cs2SnI6_prim.cell

        castep-easytool bandplot band.bands --prefix my_calc
    """
    from . import bandstructure as bsmod

    if not bands_file.exists():
        typer.secho(f"❌ File not found: '{bands_file}'", fg=typer.colors.RED)
        raise typer.Exit(code=1)

    # Resolve .cell file
    cell_path = cell_file or bands_file.with_suffix(".cell")
    if not cell_path.exists():
        typer.secho(f"❌ Cell file not found: '{cell_path}'. Use --cell to specify it.",
                    fg=typer.colors.RED)
        raise typer.Exit(code=1)

    out_stem = prefix or bands_file.stem
    out_dir  = bands_file.parent

    typer.echo(f"📖 Reading band structure from '{bands_file}'...")
    try:
        bs = bsmod.read_bands(bands_file)
    except (ValueError, IndexError) as e:
        typer.secho(f"❌ Could not parse '{bands_file}': {e}", fg=typer.colors.RED)
        raise typer.Exit(code=1)

    try:
        vertices, non_fractional = bsmod.read_cell_kpath(cell_path)
    except (OSError, ValueError) as e:
        typer.secho(f"❌ Could not read k-path from '{cell_path}': {e}", fg=typer.colors.RED)
        raise typer.Exit(code=1)

    if not vertices:
        typer.secho(f"⚠️  No k-point path block found in '{cell_path}'; k-points will be unlabelled.",
                    fg=typer.colors.YELLOW)
    else:
        if non_fractional:
            typer.secho("⚠️  The .cell file sets a non-default k-point unit; labels assume "
                        "fractional coordinates.", fg=typer.colors.YELLOW)
        matched = bsmod.assign_labels(bs, vertices)
        if matched < len(vertices):
            typer.secho(f"⚠️  Matched only {matched} of {len(vertices)} path points in '{cell_path}' "
                        "to the .bands k-points; check that the .cell file belongs to this .bands file.",
                        fg=typer.colors.YELLOW)

    # Band gap summary
    typer.echo("\n📊 Computing band gap...")
    lines = bsmod.bandgap_report(bs, out_stem)
    gap_path = out_dir / f"{out_stem}_bandgap.txt"
    gap_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.secho(f"✅ Written: '{gap_path}'", fg=typer.colors.GREEN)
    typer.echo("\n" + "\n".join(lines))

    # Plot
    typer.echo("\n🎨 Plotting band structure...")
    fig = bsmod.plot_band_structure(bs, emin=emin, emax=emax)
    for fmt in ("pdf", "png"):
        out_path = out_dir / f"{out_stem}_bands.{fmt}"
        fig.savefig(str(out_path), dpi=300, bbox_inches="tight")
        typer.secho(f"✅ Written: '{out_path}'", fg=typer.colors.GREEN)


# ---------------------------------------------------------------------------
# Command 7: converge
# ---------------------------------------------------------------------------

@app.command()
def converge(
    castep_file: Path = typer.Argument(..., help="Path to the CASTEP geometry-optimisation output (.castep) file."),
    prefix: str = typer.Option(None, "--prefix", "-p", help="Output filename stem (default: input stem)."),
    show: bool = typer.Option(False, "--show", help="Also open the plot in an interactive window."),
):
    """
    Analyse the convergence of a CASTEP **geometry optimisation**.

    Plots enthalpy, cell volume and every convergence criterion CASTEP reports
    (dE/ion, |F|max, |dR|max, Smax) against optimisation step, with tolerances,
    and prints which criteria are still unconverged.

    **Outputs** (written next to the .castep file):

    - `<stem>_convergence.png` / `.pdf` — convergence plots
    - `<stem>_convergence.txt` — text summary

    **Examples:**

        castep-easytool converge calculation.castep

        castep-easytool converge calculation.castep --show
    """
    from . import convergence as cvmod

    if not castep_file.exists():
        typer.secho(f"❌ File not found: '{castep_file}'", fg=typer.colors.RED)
        raise typer.Exit(code=1)

    typer.echo(f"📖 Parsing '{castep_file}'...")
    text = castep_file.read_text(encoding="utf-8", errors="replace")
    try:
        data = cvmod.parse_geomopt(text)
    except ValueError as e:
        typer.secho(f"❌ Could not analyse '{castep_file}': {e}", fg=typer.colors.RED)
        raise typer.Exit(code=1)

    lines = cvmod.summary_lines(data, castep_file.name)
    out_stem = prefix or castep_file.stem
    out_dir = castep_file.parent

    txt_path = out_dir / f"{out_stem}_convergence.txt"
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo("\n" + "\n".join(lines) + "\n")
    typer.secho(f"✅ Written: '{txt_path}'", fg=typer.colors.GREEN)

    if not data["converged"]:
        typer.secho("  ⚠️  Geometry optimisation has NOT converged (no completion message in the file).",
                    fg=typer.colors.YELLOW)

    if show:
        import matplotlib.pyplot as plt
    else:
        import matplotlib
        matplotlib.use("Agg")
    fig = cvmod.plot_convergence(data, title=out_stem)
    for fmt in ("png", "pdf"):
        out_path = out_dir / f"{out_stem}_convergence.{fmt}"
        fig.savefig(str(out_path), dpi=200, bbox_inches="tight")
        typer.secho(f"✅ Written: '{out_path}'", fg=typer.colors.GREEN)
    if show:
        plt.show()


# ---------------------------------------------------------------------------
# Command 8: summary
# ---------------------------------------------------------------------------

@app.command()
def summary(
    castep_file: Path = typer.Argument(..., help="Path to a CASTEP .castep output file."),
    save: bool = typer.Option(False, "--save", "-s", help="Also write <stem>_summary.txt."),
):
    """
    Print a one-page **summary** of a CASTEP run: task, functional, cut-off,
    k-points, smearing, system size, final energies, pressure, run time, warnings
    and whether the run finished.
    """
    from . import castepinfo as ci

    if not castep_file.exists():
        typer.secho(f"❌ File not found: '{castep_file}'", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    text = castep_file.read_text(encoding="utf-8", errors="replace")
    lines = ci.summary_lines(ci.parse_summary(text), castep_file.name)
    typer.echo("\n".join(lines))
    if save:
        out = castep_file.parent / f"{castep_file.stem}_summary.txt"
        out.write_text("\n".join(lines) + "\n", encoding="utf-8")
        typer.secho(f"\n✅ Written: '{out}'", fg=typer.colors.GREEN)


# ---------------------------------------------------------------------------
# Command 9: scf
# ---------------------------------------------------------------------------

@app.command()
def scf(
    castep_file: Path = typer.Argument(..., help="Path to a CASTEP .castep output file."),
    loop: int = typer.Option(-1, "--loop", "-l", help="SCF loop to plot, counted from 1 (default: the last one)."),
    prefix: str = typer.Option(None, "--prefix", "-p", help="Output filename stem (default: input stem)."),
    show: bool = typer.Option(False, "--show", help="Also open the plot in an interactive window."),
):
    """
    Analyse **SCF convergence**: the energy gain per cycle of one SCF loop against
    the tolerance, and the number of SCF cycles every loop needed. Flags loops that
    hit the SCF-cycle limit.

    Writes `<stem>_scf.png`, `<stem>_scf.pdf` and `<stem>_scf.txt`.
    """
    from . import castepinfo as ci

    if not castep_file.exists():
        typer.secho(f"❌ File not found: '{castep_file}'", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    text = castep_file.read_text(encoding="utf-8", errors="replace")
    try:
        loops = ci.parse_scf(text)
    except ValueError as e:
        typer.secho(f"❌ Could not analyse '{castep_file}': {e}", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    if loop == 0 or loop > len(loops) or loop < -len(loops):
        typer.secho(f"❌ --loop must be between 1 and {len(loops)} (or negative to count from the end).",
                    fg=typer.colors.RED)
        raise typer.Exit(code=1)
    index = loop - 1 if loop > 0 else loop

    info = ci.parse_summary(text)
    lines = ci.scf_lines(loops, info["scf_tol_eV_atom"], info["max_scf_cycles"], castep_file.name)
    out_stem = prefix or castep_file.stem
    out_dir = castep_file.parent
    txt = out_dir / f"{out_stem}_scf.txt"
    txt.write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo("\n" + "\n".join(lines) + "\n")
    typer.secho(f"✅ Written: '{txt}'", fg=typer.colors.GREEN)

    if show:
        import matplotlib.pyplot as plt
    else:
        import matplotlib
        matplotlib.use("Agg")
    fig = ci.plot_scf(loops, index, info["scf_tol_eV_atom"], title=out_stem)
    for fmt in ("png", "pdf"):
        out_path = out_dir / f"{out_stem}_scf.{fmt}"
        fig.savefig(str(out_path), dpi=200, bbox_inches="tight")
        typer.secho(f"✅ Written: '{out_path}'", fg=typer.colors.GREEN)
    if show:
        plt.show()


# ---------------------------------------------------------------------------
# Command 10: testgen
# ---------------------------------------------------------------------------

@app.command()
def testgen(
    kind: str = typer.Argument(..., help="'cutoff' or 'kpoints'."),
    cif_file: Path = typer.Argument(..., help="Path to the input .cif file."),
    values: str = typer.Option(None, "--values", "-v",
                               help="Comma-separated values: cut-offs in eV (default 300..800) or "
                                    "k-point spacings in 1/A (default 0.10..0.03)."),
    task: str = typer.Option("SinglePoint", "--task", "-t", help="CASTEP task keyword."),
    xc: str = typer.Option(None, "--xc", help="XC functional, e.g. PBE."),
    cutoff: float = typer.Option(600.0, "--cutoff", "-e", help="Fixed cut-off (eV) for a k-point test."),
    kspacing: float = typer.Option(0.05, "--kspacing", "-k", help="Fixed k-point spacing (1/A) for a cut-off test."),
    spin: bool = typer.Option(False, "--spin", help="Spin-polarised calculations."),
    out_dir: Path = typer.Option(None, "--out", "-o", help="Output folder (default: <stem>_<kind>_test)."),
):
    """
    Generate a series of CASTEP inputs for a **cut-off energy** or **k-point**
    convergence test, plus a `run_all.sh` script. Analyse the results afterwards
    with `castep-easytool testplot`.

    Examples:

        castep-easytool testgen cutoff input.cif

        castep-easytool testgen cutoff input.cif --values 400,500,600,700,800 --xc PBESOL

        castep-easytool testgen kpoints input.cif --cutoff 700
    """
    from . import convtest as ct

    kind = kind.lower()
    if kind not in ("cutoff", "kpoints"):
        typer.secho("❌ KIND must be 'cutoff' or 'kpoints'.", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    try:
        vals = ct.parse_values(values, kind)
    except ValueError as e:
        typer.secho(f"❌ {e}", fg=typer.colors.RED)
        raise typer.Exit(code=1)

    struct = _read_structure(cif_file)
    body = _lattice_abc_block(struct) + "\n" + _positions_frac_block(struct)
    folder = out_dir or cif_file.parent / f"{cif_file.stem}_{kind}_test"
    seeds = ct.write_series(kind, cif_file.stem, body, vals, folder, task=task, xc=xc,
                            fixed_cutoff=cutoff, fixed_kspacing=kspacing, spin=spin)
    typer.secho(f"✅ Wrote {len(seeds)} inputs and run_all.sh to '{folder}'", fg=typer.colors.GREEN)
    typer.echo(f"   Run them (e.g. `bash run_all.sh`), then: castep-easytool testplot {kind} {folder}")


# ---------------------------------------------------------------------------
# Command 11: testplot
# ---------------------------------------------------------------------------

@app.command()
def testplot(
    kind: str = typer.Argument(..., help="'cutoff' or 'kpoints'."),
    paths: List[Path] = typer.Argument(..., help="Folder(s) and/or .castep files from the test series."),
    tol: float = typer.Option(1.0, "--tol", help="Convergence threshold in meV/atom."),
    prefix: str = typer.Option(None, "--prefix", "-p", help="Output filename stem (default: <kind>_test)."),
):
    """
    Analyse a **cut-off or k-point convergence test**: reads the final energy of each
    finished .castep file, plots energy per atom against the tested setting and
    reports the smallest setting converged to `--tol` meV/atom.

    Writes `<stem>.png`, `<stem>.pdf` and `<stem>.txt` next to the first input.
    """
    from . import convtest as ct

    kind = kind.lower()
    if kind not in ("cutoff", "kpoints"):
        typer.secho("❌ KIND must be 'cutoff' or 'kpoints'.", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    files = []
    for pth in paths:
        if pth.is_dir():
            files += sorted(pth.glob("*.castep"))
        elif pth.exists():
            files.append(pth)
        else:
            typer.secho(f"❌ Not found: '{pth}'", fg=typer.colors.RED)
            raise typer.Exit(code=1)
    results, skipped = ct.collect(files)
    if len(results) < 2:
        typer.secho(f"❌ Need at least two finished .castep files; found {len(results)}.", fg=typer.colors.RED)
        for name, why in skipped:
            typer.echo(f"   skipped {name}: {why}")
        raise typer.Exit(code=1)

    rows, rec = ct.analyse(kind, results, tol)
    lines = ct.report_lines(kind, rows, rec, tol, skipped)
    typer.echo("\n".join(lines) + "\n")

    out_dir = (paths[0] if paths[0].is_dir() else paths[0].parent)
    stem = prefix or f"{kind}_test"
    (out_dir / f"{stem}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.secho(f"✅ Written: '{out_dir / (stem + '.txt')}'", fg=typer.colors.GREEN)
    import matplotlib
    matplotlib.use("Agg")
    fig = ct.plot(kind, rows, tol)
    for fmt in ("png", "pdf"):
        out_path = out_dir / f"{stem}.{fmt}"
        fig.savefig(str(out_path), dpi=200, bbox_inches="tight")
        typer.secho(f"✅ Written: '{out_path}'", fg=typer.colors.GREEN)


# ---------------------------------------------------------------------------
# Command 12: autoconv
# ---------------------------------------------------------------------------

def _echo_messages(messages):
    for level, text in messages:
        if level == "warn":
            typer.secho(f"⚠️  {text}", fg=typer.colors.YELLOW)
        else:
            typer.echo(f"ℹ️  {text}")


def _run_convergence(acmod, cfg, conv_dir, workdir, phase, ecut, plot_only):
    try:
        result = acmod.run_adaptive(cfg, conv_dir, workdir=workdir, phase=phase,
                                    ecut_fixed=ecut, plot_only=plot_only)
    except (FileNotFoundError, ValueError) as e:
        typer.secho(f"❌ {e}", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    except KeyboardInterrupt:
        typer.secho("\n⏹  Interrupted. Re-run the same command to resume (finished jobs are kept).",
                    fg=typer.colors.YELLOW)
        raise typer.Exit(code=130)
    typer.echo("\n" + result["report"])
    typer.secho(f"\n✅ Report: '{result['report_path']}'\n✅ Plot  : '{result['plot']}'", fg=typer.colors.GREEN)
    if result["problems"]:
        raise typer.Exit(code=2)


@app.command()
def autoconv(
    source: Path = typer.Argument(..., help="A structure (.cif or .cell) to set up and converge, or an existing .conv config."),
    xc: str = typer.Option("PBE", "--xc", help="XC functional for a new .param (structure input)."),
    spin: bool = typer.Option(False, "--spin", help="Spin-polarised calculation (structure input)."),
    cmd: str = typer.Option(None, "--cmd", help="CASTEP command, e.g. \"mpirun -np 24 castep.mpi <seedname>\". Remembered for next time."),
    primitive: Optional[bool] = typer.Option(None, "--primitive/--no-primitive",
                                             help="Convert to the primitive cell (default: automatic for .cif, warn only for .cell)."),
    threshold: float = typer.Option(1.0, "--threshold", help="Energy tolerance in meV/atom (structure input)."),
    out: Path = typer.Option(None, "--out", "-o", help="Folder for the generated inputs (default: <stem>_convergence)."),
    yes: bool = typer.Option(False, "--yes", "-y", help="Start running without asking for confirmation."),
    setup_only: bool = typer.Option(False, "--setup-only", help="Write the .cell/.param/.conv files and stop."),
    init: bool = typer.Option(False, "--init", help="Write a template .conv file at SOURCE and exit."),
    phase: str = typer.Option("both", "--phase", help="'both', 'ecut' (cut-off only) or 'kpoints' (k-points only)."),
    ecut: float = typer.Option(None, "--ecut", help="Fixed cut-off (eV) for a k-point-only run."),
    workdir: Path = typer.Option(None, "--workdir", "-w", help="Folder for the CASTEP jobs (default: <conv folder>/serial)."),
    plot_only: bool = typer.Option(False, "--plot-only", help="Do not run CASTEP: analyse and plot finished jobs only."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Show the analysis and planned scans, write and run nothing."),
):
    """
    **Automatic cut-off and k-point convergence** for CASTEP.

    Give it a structure and it does the rest: checks for a primitive cell (and converts
    a conventional `.cif`, telling you), writes `<seed>.cell`, `<seed>.param` and a
    tailored `<seed>.conv`, then runs
    CASTEP job by job until the energy (and stress, for optimisations) converge.

    The k-point mesh follows the cell shape (works for hexagonal, tetragonal,
    orthorhombic and low-symmetry cells) and the result reports both the converged
    `KPOINTS_MP_GRID` and the equivalent `KPOINTS_MP_SPACING`. Resumable.

    **Examples:**

        castep-easytool autoconv structure.cif --cmd "mpirun -np 24 castep.mpi <seedname>"

        castep-easytool autoconv structure.cif --setup-only

        castep-easytool autoconv structure_convergence/structure.conv

        castep-easytool autoconv structure.conv --phase kpoints --ecut 700
    """
    import tempfile
    from . import autoconv as acmod
    from . import convsetup as csmod

    if phase not in ("both", "ecut", "kpoints"):
        typer.secho("❌ --phase must be 'both', 'ecut' or 'kpoints'.", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    suffix = source.suffix.lower()

    # ---- template config -------------------------------------------------
    if init:
        if suffix != ".conv":
            typer.secho("❌ --init writes a .conv template: give a name ending in .conv.", fg=typer.colors.RED)
            raise typer.Exit(code=1)
        if source.exists():
            typer.secho(f"❌ '{source}' already exists; not overwriting.", fg=typer.colors.RED)
            raise typer.Exit(code=1)
        source.write_text(acmod.render_conv(acmod.default_config(source.stem)), encoding="utf-8")
        typer.secho(f"✅ Wrote template config to '{source}'. Edit it, then run: castep-easytool autoconv {source}",
                    fg=typer.colors.GREEN)
        return

    # ---- structure input: set up, then continue like a .conv ----------------
    if suffix in (".cif", ".cell"):
        typer.echo(f"📖 Reading '{source}'...")
        tmp = tempfile.TemporaryDirectory() if dry_run else None
        try:
            try:
                prep = csmod.prepare(source, out_dir=(Path(tmp.name) / source.stem) if dry_run else out,
                                     xc=xc, spin=spin, cmd=cmd, primitive=primitive, threshold=threshold,
                                     log=typer.echo)
            except (FileNotFoundError, ValueError, OSError) as e:
                typer.secho(f"❌ {e}", fg=typer.colors.RED)
                raise typer.Exit(code=1)
            _echo_messages(prep["messages"])
            cfg, conv_dir, conv_path = prep["cfg"], prep["out_dir"], prep["conv_path"]
            try:
                pl = acmod.plan(cfg, conv_dir, phase)
            except (FileNotFoundError, ValueError) as e:
                typer.secho(f"❌ {e}", fg=typer.colors.RED)
                raise typer.Exit(code=1)
            typer.echo("\n" + "\n".join(acmod.plan_lines(cfg, pl, conv_dir if not dry_run else "(dry run: nothing written)")))
        finally:
            if tmp:
                tmp.cleanup()
        if dry_run:
            return
        if cmd:
            csmod.save_command(cmd)
        typer.secho(f"\n✅ Set up '{conv_path}'", fg=typer.colors.GREEN)
        if setup_only:
            typer.echo(f"   Review the files, then run: castep-easytool autoconv {conv_path}")
            return
        if not yes and not typer.confirm("\nStart the CASTEP jobs now?", default=False):
            typer.echo(f"Not started. Run later with: castep-easytool autoconv {conv_path}")
            return
        _run_convergence(acmod, cfg, conv_dir, workdir, phase, ecut, plot_only)
        return

    # ---- existing .conv --------------------------------------------------
    if suffix != ".conv":
        typer.secho(f"❌ Unsupported input '{source.name}': give a .cif, .cell or .conv file.", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    if not source.exists():
        typer.secho(f"❌ Config file not found: '{source}' (create one with --init, or pass a structure).",
                    fg=typer.colors.RED)
        raise typer.Exit(code=1)
    try:
        cfg, warns = acmod.parse_config(source)
    except ValueError as e:
        typer.secho(f"❌ Invalid config '{source}': {e}", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    for w in warns:
        typer.secho(f"⚠️  {source.name}: {w}", fg=typer.colors.YELLOW)
    if cmd:
        cfg["castep_command"] = cmd
        csmod.save_command(cmd)
    conv_dir = source.resolve().parent
    if dry_run:
        try:
            pl = acmod.plan(cfg, conv_dir, phase)
        except (FileNotFoundError, ValueError) as e:
            typer.secho(f"❌ {e}", fg=typer.colors.RED)
            raise typer.Exit(code=1)
        typer.echo("\n".join(acmod.plan_lines(cfg, pl, conv_dir)))
        return
    _run_convergence(acmod, cfg, conv_dir, workdir, phase, ecut, plot_only)


if __name__ == "__main__":
    app()
