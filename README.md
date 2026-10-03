# RasMol 2.6-beta-2 — Modernized

This is [Roger Sayle](mailto:rasmol@ggr.co.uk)'s **RasMol 2.6-beta-2**
(`RasMol26b2.tar.gz`), the classic molecular graphics visualisation
program for proteins, nucleic acids and small molecules, patched up to
build and run on current Linux/macOS with a minimal dependency
footprint, plus a native mmCIF/PDBx reader it never had.

It is **not** affiliated with, and should not be confused with,
[OpenRasMol](http://www.openrasmol.org/) — the long-running, separately
maintained GPL continuation of RasMol (currently in the 2.7.x series).
This repository starts from a different, earlier snapshot (2.6-beta-2)
and has a much narrower goal: keep that specific codebase buildable and
correct on modern systems, not extend it into a new release line.

## What RasMol is

RasMol reads molecular coordinate files and interactively displays the
molecule as wireframe, sticks, alpha-carbon trace, CPK spacefilling
spheres, ribbons/cartoons, hydrogen bonds and dot surfaces, with
rotation, zoom, slabbing, a scripting language, and image export to
several raster/vector formats. The full original description, feature
list and supported-platform history are preserved verbatim in
[`README.RASMOL-ORIGINAL`](README.RASMOL-ORIGINAL).

## What's changed from stock 2.6-beta-2

- **SDL2 display backend** (`sdlwin.c`) replacing the original X11/Xlib
  backend (`x11win.c`) — no X server required. Mouse-driven
  rotate/translate/zoom/slab, click-to-pick and window resizing all
  still work; the hand-drawn X11 menu/scrollbar/dials-box chrome was
  dropped in favour of RasMol's own command language.
- **64-bit correctness**: `Card`/`Long`/`Pixel` are pinned to 32 bits
  via `_LONGLONG` (`rasmol.h`), fixing pixel-packing corruption that
  plain `long` being 64 bits on modern LP64 systems silently caused.
- **mmCIF / PDBx reader** (`LoadCIFMolecule`, in `infile.c`) — RasMol
  2.6-beta-2 could only *write* CIF (`SaveCIFMolecule`), never read it.
  This adds a small dependency-free tag/loop tokeniser for the subset
  of the mmCIF grammar the PDB and AlphaFold DB actually emit, and
  structurally parses the `_atom_site` loop by column *name* rather
  than fixed column *position* (mmCIF is tag/loop based, not
  column-based like PDB). Atom/element typing is built to land in the
  same internal representation the PDB parser produces, so CIF-loaded
  structures get the same backbone tracing, cartoons and element
  colouring as PDB ones.
- **Format auto-detection**: a bare `load foo.cif` or `rasmol foo.cif`
  now correctly picks the mmCIF reader from the `.cif`/`.mmcif`
  extension. Previously (and this bit without the fix) it silently fed
  the file to the fixed-column PDB parser, which misreads free-form
  mmCIF text as 80-column fields and produces garbage, visibly
  flattened coordinates. Explicit `-cif` / `load cif ...` still work
  as before.
- **Bug fixes** found while testing the above:
  - A double-free that aborted (`SIGABRT`) almost any `-script` run —
    `main()` was closing the script file handle a second time after
    `LoadScriptFile()` had already closed it on every return path.
  - An infinite, 100%-CPU busy loop when stdin reaches EOF without an
    `exit`/`quit` ever being issued (piped/automated/headless use, or
    Ctrl-D at a real terminal) — the main command loop never checked
    `ReadCharacter()`'s EOF, and on this build's ABI the loop variable
    being `char` (unsigned here) silently turned EOF (`-1`) into `255`
    before any check could even see it.
  - Several implicit-function-declaration, maybe-uninitialized and
    other build warnings; `make` is now warning-free.
- Full dated details of every change — both this modernization work
  and the original development history back to 1992 — are in
  [`ChangeLog`](ChangeLog).

## Provenance

- Base source: `RasMol26b2.tar.gz` ("RasMol 2.6-beta-2"), from Roger
  Sayle's distribution via the University of Edinburgh, as indexed on
  the [UMass RasMol source page](https://www.umass.edu/microbio/rasmol/srccode.htm).
  Fetch that tarball directly from the link above if you want to diff
  this repository against the pristine upstream source.
- Original author: **Roger Sayle**, Biomolecular Structures Group,
  Glaxo Wellcome Research & Development (RasMol was originally
  developed at the University of Edinburgh's Biocomputing Research
  Unit). See [`README.RASMOL-ORIGINAL`](README.RASMOL-ORIGINAL) and
  [`doc/`](doc/) for the original README, manual and reference card.

## Building

Requires a C compiler and SDL2 development headers.

```sh
# Debian/Ubuntu
sudo apt install libsdl2-dev pkg-config
make
```

For other platforms, see [`INSTALL`](INSTALL) (original build notes;
the SDL2/Makefile path is specific to this fork and targets
Linux/macOS).

## Running

```sh
./rasmol data/1crn.pdb                     # PDB
./rasmol somestructure.cif                 # mmCIF, auto-detected from the extension
./rasmol -cif somefile.txt                 # force mmCIF when the extension doesn't give it away
./rasmol -nodisplay -script myscript.rsc   # headless / scripted use
```

See `doc/rasmol.txt` / `rasmol.hlp` for the full command language.

### Known mmCIF reader limitations

- Chain identifiers are truncated to one character — `Chain.ident` is
  a plain `char` throughout the whole codebase, a PDB-era limitation
  this fork hasn't lifted.
- No `_struct_conf` (HELIX/SHEET) parsing yet, so secondary structure
  reads as "No Assignment" until you ask RasMol to calculate it.
- Multi-model mmCIF files load only the first `pdbx_PDB_model_num`.

## Test fixtures

`data/` contains a couple of small PDB-format structures used to
exercise the original parser. The mmCIF reader was validated during
development against real PDB and AlphaFold DB entries (a liganded
structure, a two-chain structure, and a large AlphaFold prediction) —
see `ChangeLog` for specifics; grab any `.cif` from
[RCSB PDB](https://www.rcsb.org/) or the
[AlphaFold DB](https://alphafold.ebi.ac.uk/) to exercise it yourself.

## License

RasMol's original terms, which this entire repository — including all
modifications — continues to follow:

> The source code is public domain and freely distributable provided
> that the original author is suitably acknowledged.

See [`LICENSE`](LICENSE) for the full statement and provenance.
