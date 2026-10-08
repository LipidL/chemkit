# chemkit

A small, type-safe Rust toolkit for reading, writing, and manipulating
atomistic structures used in computational chemistry.

`chemkit` parses common structure file formats into a single canonical
in-memory representation, so downstream code does not have to care which
format a structure came from. Units are checked at compile time via
[`uom`](https://docs.rs/uom), and geometry is backed by
[`nalgebra`](https://docs.rs/nalgebra).

## Features

- **Compile-time unit safety.** Lengths and angles are `uom` quantities;
  you choose the unit you want at the call site, and dimension mismatches
  fail to compile instead of producing silent bugs.
- **Generic precision.** The public API is generic over `T`, so parsing and
  writing work with `f32`, `f64`, or any other `nalgebra::RealField` that
  supports the required unit conversions.
- **Multi-frame support.** Formats that hold several structures per file
  (such as ARC) parse to a `Vec<Structure<T>>` and round-trip losslessly.

## Supported formats

| Format            | Read | Write | Notes |
| ----------------- | :--: | :---: | ----- |
| ARC (Biosym/MS)   |  ✓  |  ✓   | Multi-frame; reads and writes `PBC` cell parameters, atom records, and per-frame `index` / `Q` / `energy` metadata. |
| POSCAR (VASP)     |  ✓  |  ✗    | Single frame. Supports direct and volume-based scaling factors, species names or bare ion counts (with a matching POTCAR), `Selective Dynamics`, and Cartesian or fractional coordinates. |

## The data model

All parsers converge on `Structure<T>`:

```text
Structure<T>
├── atoms: Vec<Atom<T>>
│   ├── element: Element (name, atomic number, mass, radii, …)
│   └── position: CoordinateSystem<T>   // Cartesian or Spherical
├── cell: Option<CellParameters<T>>     // None for non-periodic structures
│   └── a, b, c, alpha, beta, gamma
└── properties: HashMap<String, String> // format-specific metadata
```

Key points:

- **Positions are always stored as Cartesian or spherical `uom` values.**
  Fractional coordinates from POSCAR are converted to Cartesian during
  parsing, using the cell that is read from the same file.
- **`properties` carries format-specific metadata** as string key/value
  pairs, so the common model stays format-neutral. The ARC parser writes
  `arc_index`, `arc_q_value`, and `arc_energy`; the POSCAR parser stores
  the comment line under `comment`. The ARC writer reads those keys back
  and falls back to `0` / `0.0` when they are absent.

## Notes and limitations

- Writing POSCAR is not implemented at the current version; only ARC can be serialised.
- The ARC parser silently skips unrecognised lines and returns an empty
  vector for truncated input. The only hard error is an atom line whose
  element symbol is not in the periodic table.
- `PeriodicTable` covers the elements needed by the parsers; reverse lookup
  by atomic number is a linear scan, which is fine for typical structures.

## Development

```sh
cargo test      # run the test suite (unit tests + ARC round-trip)
cargo build     # build the library
```
