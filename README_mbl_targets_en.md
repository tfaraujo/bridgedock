# Target table generation — `mbl_targets.csv`

This document describes how `mbl_targets.csv`, the target table consumed by
`mbl_dock.py`, was produced for the metallo-β-lactamase (MBL) docking study; which
of its fields are geometric measurements and which are human curation; and the
known limitations of the scripts included here.

The goal is honest reproducibility: anyone with the PDB files should be able to
recover the published numbers and identify unambiguously which fields depended on
author judgement.

---

## 1. What the table contains

One row per target, holding the parameters `mbl_dock.py` needs to prepare the
receptor and centre the search box:

| Column | Meaning |
|---|---|
| `target` | Enzyme name (NDM-1, VIM-2, SPM-1, …) |
| `pdb` | PDB code of the selected structure |
| `chain` | Chain used for docking |
| `keep_zn` | `resSeq` of the two catalytic Zn²⁺ ions retained |
| `drop_zn` | `resSeq` of additional Zn²⁺ ions (crystallographic/adventitious) removed |
| `remove_het` | HETATM records stripped from the receptor (solvents, buffers, cryoprotectants, co-crystallised inhibitors) |
| `center_x/y/z` | Search-box centre — midpoint between the two retained Zn²⁺ ions |
| `znzn_A` | Zn1–Zn2 distance in Å |
| `status` | Curation note (rationale for chain, structure, and removals) |

Targets included: NDM-1 (8B1W), NDM-5 (6MGY), NDM-7 (7AEZ), VIM-1 (5N5G),
VIM-2 (7A5Z), VIM-7 (2Y87), IMP-1 (7YH9), SPM-1 (5NDB). All subclass B1,
dinuclear MBLs.

Note that the values in the `status` column of the published CSV are written in
Portuguese, as recorded during curation.

---

## 2. How the table was actually produced

The table was assembled target by target using ad hoc structural inspection
scripts: for each PDB, these listed the Zn²⁺ ions present, the distances between
them, the donor atoms (N, O, S) within 2.9 Å of each ion, and the HETATM records in
the file. From that survey the author decided the chain, the catalytic pair, and
the zinc ions and heteroatoms to remove, recording the rationale in `status`.

**Provenance notice.** The set of scripts that generated the final table was not
preserved. What survived is included here:

- **`gerador_5NDB_SPM-1.py`** — surviving fragment of that workflow, in the form in
  which it was applied to SPM-1 (5NDB), with residual checks for VIM-1 and VIM-2.
  It prints the coordination sphere of each Zn²⁺ and the midpoint of the pair. It
  writes no CSV and has hard-coded filenames. It is documentary, not production
  code.
- **`make_mbl_targets.py`** — a **later reconstruction**, written by reverse
  engineering the published table. It automates the purely geometric part of the
  original procedure and applies it uniformly across targets.

`make_mbl_targets.py` **is not** the script that generated `mbl_targets.csv`, and
it does not reproduce it in full. The differences are listed in section 4; they are
known, not accidental.

---

## 3. Provenance by column

| Column | Origin | Reproducible by `make_mbl_targets.py` |
|---|---|---|
| `target`, `pdb` | Literature-based selection | — (fixed mapping in the script) |
| `keep_zn` | Geometry: Zn pair in the same chain with 2.8 ≤ d ≤ 5.0 Å | **Yes** |
| `center_x/y/z` | Geometry: midpoint of the pair | **Yes** |
| `znzn_A` | Geometry: pair distance | **Yes** |
| `chain` | Heuristic + human review | Partly (see 4.1) |
| `drop_zn` | Curation | **No** — column not emitted |
| `remove_het` | Curation over the file inventory | **No** — see 4.2 |
| `status` | Literature and author judgement | **No** — by nature |

---

## 4. Known divergences of the reconstruction

Running

```bash
python make_mbl_targets.py --pdb_dir ./pdb --out mbl_targets_auto.csv
```

over the eight PDB files reproduces the published geometric columns for **seven of
the eight targets**. The divergences are as follows.

### 4.1 SPM-1 (5NDB) — chain selection

`make_mbl_targets.py` breaks ties between valid Zn pairs alphabetically and returns
chain A. The published table uses chain **B**:

| | Chain | Zn | Centre (x, y, z) |
|---|---|---|---|
| Published (curated) | B | 401, 402 | −43.106 / 16.023 / 8.747 |
| `make_mbl_targets.py` | A | 401, 402 | −45.690 / 21.040 / −12.390 |

5NDB contains two copies in the asymmetric unit, each with three Zn²⁺ ions
(A: 401, 402, 407; B: 401, 402, 406). Chain B was chosen by inspection of the site,
not by alphabetical ordering. Alphabetical ordering is an arbitrary criterion and is
flagged as such — the correct alternative (actual site quality: number of N/O/S
donors within 2.9 Å of each Zn and mean B-factor of the first coordination shell)
is planned for v2.

**To reproduce the published row, the SPM-1 chain must be set manually.**

### 4.2 `remove_het` — fixed list instead of inventory

`make_mbl_targets.py` writes the same list of common solvents and buffers for every
target (`HOH, GOL, EDO, SO4, PO4, MES, TRS, ACT, DMS, PEG`). The published table
lists the HETATM records actually present in each structure, which includes the
co-crystallised inhibitors: `OQU` (8B1W), `R8W` (7AEZ), `BCN` (5N5G), `QZH` (7A5Z),
`IT0` (7YH9), `8TW` (5NDB).

**Practical consequence:** used without review, the automated output leaves the
crystallographic ligand inside the receptor, and the search box lands on an
occupied site. This is not cosmetic — it invalidates docking for the affected
target.

### 4.3 `drop_zn` — not emitted

Zinc ions beyond the catalytic pair are not identified by the reconstruction:
VIM-1 (303), VIM-2 (302) and SPM-1 (406/407) were annotated manually. Without this
column, a third ion remains in the receptor.

### 4.4 Column name

The reconstruction writes `zn_zn_dist`; the published table uses `znzn_A`. Rename
before feeding the file to `mbl_dock.py`.

---

## 5. Current usage

1. Download the PDB files listed in section 1 into a directory.
2. Run `make_mbl_targets.py` to obtain the geometric columns.
3. Check the output against the published `mbl_targets.csv`.
4. Fill in `drop_zn`, `remove_het` and `status`, and correct the SPM-1 chain,
   following the published table.

For the target set of this study, the recommended path is simply to **use the
`mbl_targets.csv` versioned in this repository**, which is the table actually used
in the calculations. `make_mbl_targets.py` is there to audit the geometric columns
and to extend the procedure to new targets, with manual review of the curated
fields.

---

## 6. Future work (`scaffold_targets.py`, v2)

This version documents the procedure as it was carried out. A generalised version,
to be run before `mbl_dock.py` on arbitrary targets, is planned for after
publication:

- **`--fetch 8B1W,6MGY`** — direct download from the RCSB, so that PDB codes are
  the only required input.
- **`--mononuclear`** — support for subclass B2 MBLs (CphA, Sfh-I), which carry a
  single catalytic Zn²⁺ and are currently reported as failures by the dinuclear
  pair criterion.
- **Chain selection by site quality** — counting N/O/S donors within 2.9 Å of each
  Zn and the mean B-factor of the first shell, replacing alphabetical ordering.
- **HETATM inventory and classification** — solvent/buffer vs. ligand, from a
  reference list combined with heavy-atom count and distance to the site.
- **Automatic `drop_zn`** — Zn²⁺ ions in the selected chain outside the catalytic
  pair.
- **`status` as a factual diagnostic** — number of chains, surplus Zn ions, bulky
  HET groups near the site, incomplete coordination. The field remains
  human-authored: notes such as *"5NDB replaces oxidised 2FHX"* come from the
  literature, not from the file, and no heuristic should fabricate them.

---

## 7. Stated limitations

- Chain selection and HETATM classification are heuristics and require review. No
  future version will remove that requirement; v2 will only supply the diagnostic
  that makes it faster and traceable.
- Apo structures, structures with partial Zn occupancy, and metal-substituted forms
  (Cd, Co) are not handled. The 2.8–5.0 Å criterion assumes an intact dinuclear
  site.
- Alternate conformers: only blank or `A` `altLoc` records are read.
- The table reflects the state of the structures on the date of access to the PDB;
  deposition revisions may change residue numbering.
