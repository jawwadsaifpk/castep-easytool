# castep-easytool

A command-line tool to simplify common pre- and post-processing tasks for the [CASTEP](http://www.castep.org/) plane-wave DFT code.

## Installation

```bash
pip install -e .
```

## Commands

### `convert` — Convert between primitive and conventional cells

```bash
castep-easytool convert input.cif --to primitive
castep-easytool convert input.cif --to conventional
```

Reads any structure file supported by pymatgen (`.cif`, etc.) and writes the standardised primitive or conventional cell as a `.cif` file.

---

### `generate` — Generate a CASTEP `.cell` input file

```bash
castep-easytool generate input.cif
castep-easytool generate input.cif --task SinglePoint
```

Reads a `.cif` file and writes a CASTEP `.cell` file with `LATTICE_ABC` and `POSITIONS_FRAC` blocks. The default task is `GeometryOptimization`.

---

### `geomopt` — Extract the final geometry from a geometry optimisation

```bash
castep-easytool geomopt calculation.castep
castep-easytool geomopt calculation.castep --prefix my_output
```

Parses the final optimised geometry from a CASTEP `.castep` output file (supports LBFGS, BFGS, and FIRE optimisers) and writes:

- `<stem>_opted.cell` — ready-to-use `.cell` file for the next calculation
- `<stem>_opted.cif` — CIF file of the optimised structure

Warns if the geometry optimisation did not converge, but still writes the last available geometry.

---

### `supercell` — Build a supercell

```bash
castep-easytool supercell input.cif 2
castep-easytool supercell input.cif 2 2 1
```

Scales the unit cell by `na x nb x nc` along the three lattice directions and writes a CASTEP `.cell` file. If only one scale factor is given, it is applied to all three axes.

---

### `bands` — Generate a k-path `.cell` file for a band structure calculation

```bash
castep-easytool bands input.cif
```

Uses [SeeK-path](https://github.com/giovannipizzi/seekpath) to generate a high-symmetry Brillouin zone path and writes a `<stem>_bands.cell` file with the `SPECTRAL_KPOINT_PATH` block. The structure is automatically standardised to the seekpath primitive cell.

---

### `bandplot` — Plot a CASTEP band structure

```bash
castep-easytool bandplot calculation.bands
castep-easytool bandplot calculation.bands --emin -3 --emax 3
castep-easytool bandplot calculation.bands --cell my_structure.cell
castep-easytool bandplot calculation.bands --prefix my_output
```

Reads a `.bands` file from a completed CASTEP spectral/bands calculation and writes:

- `<stem>_bands.pdf` — publication-quality vector plot
- `<stem>_bands.png` — raster plot for quick viewing
- `<stem>_bandgap.txt` — band gap, direct/indirect character, and transition k-points (spin-resolved for spin-polarised calculations)

Requires [sumo](https://github.com/SMTG-Bham/sumo): `pip install sumo`

---

## Typical workflow

```bash
# 1. Start from a CIF file and generate a CASTEP input
castep-easytool generate structure.cif

# 2. Run the geometry optimisation in CASTEP, then extract the result
castep-easytool geomopt structure.castep

# 3. Build a supercell from the optimised geometry
castep-easytool supercell structure_opted.cif 2 2 2

# 4. Generate a k-path for a band structure calculation
castep-easytool bands structure_opted.cif

# 5. After the bands calculation, plot the result
castep-easytool bandplot structure_opted.bands
```

## Dependencies

| Package | Purpose |
|---|---|
| [typer](https://typer.tiangolo.com/) | CLI framework |
| [pymatgen](https://pymatgen.org/) | Structure reading, writing, and manipulation |
| [seekpath](https://github.com/giovannipizzi/seekpath) | BZ k-path generation (`bands` command) |
| [sumo](https://github.com/SMTG-Bham/sumo) | Band structure parsing and plotting (`bandplot` command) |
