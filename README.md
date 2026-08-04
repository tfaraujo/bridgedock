# bridgedock

A docking pipeline and a **coordination-aware bimetallic bridging filter** for di-nuclear metalloenzymes, applied here to subclass B1 metallo-β-lactamases (MBLs).

The point of this repository is a single methodological claim: **geometric proximity alone is not an acceptable operational definition of metal coordination in docking.** A criterion that asks only "is some ligand atom within 3 Å of both metals?" classifies saturated hydrocarbons — which cannot coordinate Zn²⁺ at all — as zinc-bridging ligands. Restricting the criterion to heteroatoms that can actually form inner-sphere bonds to Zn²⁺ (S, N, O) removes every one of those false positives while retaining genuine coordination.

This code and data accompany the manuscript *"Heteroatom-restricted bimetallic bridging criteria reveal that garlic-derived organosulfur compounds, but not saturated hydrocarbons, are predicted to span the di-zinc centre of subclass B1 metallo-β-lactamases"* (submitted to *BMC Research Notes*).

---

## What the filter does

A pose is classified as **bridging** when at least one ligand atom capable of coordinating Zn²⁺ lies within 3.0 Å of Zn1, and at least one such atom lies within 3.0 Å of Zn2 (the same atom may satisfy both conditions).

Coordinating atoms are restricted to **S, N and O**. Carbon and hydrogen are excluded because they do not form inner-sphere bonds to Zn²⁺. The 3.0 Å threshold is permissive relative to observed coordination distances (Zn–S ≈ 2.3 Å, Zn–N ≈ 2.1 Å, Zn–O ≈ 2.0 Å).

Pass `--coord_elements any` to reproduce the proximity-only behaviour for comparison.

### Validation

The filter was tested in both directions, which is the part that matters:

| Control | Result |
|---|---|
| **Negative** — 3 saturated hydrocarbons × 8 enzymes | 0/24 pairs bridge (0/215 poses) |
| **Positive, experimental** — crystallographic poses of 5 co-crystallised ligands | the 2 that coordinate both zincs are accepted; the 3 that contact one zinc or none are rejected |
| **Positive, docked** — 4 thiol inhibitors × 8 enzymes | 27/32 pairs bridge; captopril bridges in all 8, in 5 through sulfur |
| Study compounds — 6 organosulfur × 8 enzymes | 37/48 pairs bridge (109/430 poses) |

Thiol inhibitors and organosulfur compounds bridge at statistically indistinguishable rates (25.7% vs 25.3%, Fisher p = 0.93); both differ from hydrocarbons at p ≈ 10⁻²⁰.

Redocking of the 5 co-crystallised ligands reproduced the deposited pose within 2.0 Å in 3 cases (median 1.74 Å). Both failures involve ligands with weak or absent zinc anchoring, and in both the near-correct pose was generated but ranked second — sampling found the crystallographic mode, scoring failed to rank it first.

---

## Contents

### Code and environment

| File | Description |
|---|---|
| `mbl_dock.py` | the whole pipeline: receptor preparation, docking, bridging filter, redocking validation, reanalysis |
| `environment.yml` | conda environment (Python 3.10, numpy, scipy, RDKit, Open Babel) |
| `LICENSE` | MIT |
| `CITATION.cff` | citation metadata |

### Inputs

| File | Description |
|---|---|
| `mbl_targets.csv` | the 8 receptors: PDB code, chain, catalytic Zn residue numbers, HETATM codes removed, grid centre, Zn–Zn distance, preparation notes |
| `ligands.csv` | the 9 study compounds (6 organosulfur + 3 hydrocarbon negative controls), as SMILES |
| `ligands_positive_controls.csv` | the 4 thiol positive controls (captopril, thiorphan, tiopronin, dimercaprol), as SMILES |

Receptor coordinates are **not** included — download them from the PDB using the codes in `mbl_targets.csv`: 8B1W (NDM-1), 6MGY (NDM-5), 7AEZ (NDM-7), 5N5G (VIM-1), 7A5Z (VIM-2), 2Y87 (VIM-7), 7YH9 (IMP-1), 5NDB (SPM-1).

### Results

| File | Description |
|---|---|
| `results_main_panel.csv` | primary output, 72 compound–target pairs. Columns ending `_FIXED` are the heteroatom-restricted filter; `_ORIGINAL` are the proximity-only filter, retained so the difference between the two criteria can be recomputed |
| `results_thiol_controls.csv` | the 4 thiol inhibitors × 8 enzymes, same format as the pipeline's standard output |
| `redock_validation.csv` | redocking of the 5 co-crystallised ligands: RMSD of the top-ranked and closest poses, RMSD method, and whether the crystallographic pose itself satisfies the bridging criterion |

### Manuscript tables and figures

| File | Description |
|---|---|
| `Table1_docking_MBL.csv` | Table 1 as published: ΔG_bridge and bridging pose counts per compound and target |
| `Table_S1_full_docking_results.csv` | Additional file, Table S1: all 72 pairs in long format, with both criteria side by side |
| `Table_S4_redocking_validation.csv` | Additional file, Table S4: redocking validation, formatted |
| `Table_S5_thiol_positive_controls.csv` | Additional file, Table S5: thiol positive controls, formatted |
| `Figure1_heatmap.png` / `.svg` | Figure 1: bridging coordination across all compounds and targets |
| `Figure2_allicin_IMP1_2D.png` / `_3D.png` | Figure 2: allicin in IMP-1, best bridging pose |
| `Figure2_thiophene_NDM7_2D.png` / `_3D.png` | Figure 2: thiophene-2-carboxylic acid in NDM-7, best bridging pose |

Two-dimensional interaction diagrams were generated in BIOVIA Discovery Studio Visualizer. Note that its default metal-interaction cutoff is shorter than the 3.0 Å used here, so a coordination contact reported in the tables may not be drawn in the diagram; heteroatom–Zn distances in the tables are authoritative.

---

## Requirements

External tools, none bundled:

| Tool | Notes |
|---|---|
| [AutoDock Vina 1.2](https://github.com/ccsb-scripps/AutoDock-Vina) | run in AutoDock4 scoring mode (`--scoring ad4`) |
| AutoGrid4 | from the AutoDock 4.2 suite |
| [MGLTools](https://ccsb.scripps.edu/mgltools/) | provides `prepare_receptor4.py` and `pythonsh` |
| [AutoDock4Zn](https://autodock.scripps.edu/resources/autodock4zn/) | provides `zinc_pseudo.py` and `AD4Zn.dat` — **not** part of standard MGLTools |
| [Open Babel 3.x](https://openbabel.org) | `obabel` for conversion, `obrms` for symmetry-corrected RMSD |

```bash
conda env create -f environment.yml
conda activate mbl
```

`obrms` matters more than it looks: without it, RMSD is computed by atom order and is **overestimated** for ligands with symmetric groups (carboxylates, phenyl rings). A correct redocking can look like a failure. The script logs which method it used and records it in the `rmsd_method` column.

---

## Usage

### Docking

```bash
python mbl_dock.py \
  --targets mbl_targets.csv --ligands ligands.csv \
  --pdb_dir . --outdir mbl_out \
  --prepare_receptor4 /path/to/prepare_receptor4.py \
  --engine vina --exhaustiveness 16
```

Writes `mbl_out/mbl_results.csv`: one row per compound–target pair with the best score, the best **bridged** score, the two shortest heteroatom–Zn distances, which elements coordinate each zinc, and the number of bridging poses.

The output filename is fixed, so point `--outdir` at a fresh directory unless you want the previous run overwritten.

### Redocking (protocol validation)

```bash
python mbl_dock.py --redock \
  --targets mbl_targets.csv --pdb_dir . --outdir mbl_out_redock \
  --prepare_receptor4 /path/to/prepare_receptor4.py
```

Extracts each co-crystallised ligand listed in `remove_het` (skipping solvents, buffers and `UNX`/`UNK` placeholders), docks it back into its own receptor and reports symmetry-corrected RMSD against the deposited pose. Targets without a co-crystallised ligand are skipped with a log message. Use `--redock_resname` to force a specific HETATM code.

### Reanalysis (no docking)

```bash
python mbl_dock.py --reanalyze --outdir mbl_out --targets mbl_targets.csv
```

Re-applies the bridging filter to existing `out_*.pdbqt` files and writes `mbl_results_fixed.csv` with `_ORIGINAL` columns for comparison. Useful for cutoff sensitivity, since no docking is repeated:

```bash
python mbl_dock.py --reanalyze --outdir mbl_out --bridge_cut 2.5 \
  --reanalyze_out results_cut2.5.csv
```

### Reproducing the published results

```bash
PREP=/path/to/prepare_receptor4.py

# 1. main panel: 6 organosulfur compounds + 3 hydrocarbon controls
python mbl_dock.py --targets mbl_targets.csv --ligands ligands.csv \
  --pdb_dir . --outdir mbl_out --prepare_receptor4 $PREP \
  --engine vina --exhaustiveness 16

# 2. protocol validation by redocking
python mbl_dock.py --redock --targets mbl_targets.csv \
  --pdb_dir . --outdir mbl_out_redock --prepare_receptor4 $PREP

# 3. thiol positive controls
python mbl_dock.py --targets mbl_targets.csv \
  --ligands ligands_positive_controls.csv \
  --pdb_dir . --outdir mbl_out_controls --prepare_receptor4 $PREP \
  --engine vina --exhaustiveness 32 --control Captopril
```

Docking is stochastic and each pair was run once, so pose counts will differ slightly between runs. The qualitative result — zero bridging for hydrocarbons, bridging for organosulfur compounds and thiols — is not sensitive to the seed.

---

## Two implementation notes worth knowing

Both were bugs found during development. They are documented because anyone reimplementing this filter will hit them.

**1. The element field sits at a different column in PDB and PDBQT.** In PDB it occupies columns 77–78; in PDBQT the AutoDock atom type occupies columns 78–79. Reading the PDBQT offset from a PDB file turns `CL` into `L` and — much worse — turns `ZN` into `N`, so every receptor zinc read from a `.pdb` is silently classified as nitrogen. `element_of()` tries both offsets and prefers the two-letter match.

**2. AutoDock atom types are not chemical elements.** `OA`, `NA` and `SA` are the hydrogen-bond-acceptor forms of oxygen, nitrogen and sulfur. A filter that tests membership in `{"S","N","O"}` against the raw type field discards almost every acceptor oxygen and sulfur in the ligand, and therefore silently under-counts bridging. Types are mapped to elements explicitly (`OA→O`, `NA→N`, `SA→S`, `A→C`, `HD→H`).

One thing that is **not** a bug: `zinc_pseudo.py` places one TZ pseudo-atom per *vacant* tetrahedral coordination site, not four per zinc. An MBL Zn1 (3 His) and Zn2 (Asp/Cys/His) each have one vacant site, so two TZ atoms for a di-nuclear site is correct.

---

## Scope and limitations

Docking here is **non-covalent**. It does not model the redox and covalent chemistry — thiol–disulfide exchange, S-allylmercapto modification of protein cysteines — through which organosulfur compounds such as allicin exert their documented antimicrobial activity. Results indicate the plausibility of an interaction with the catalytic site and must not be read as predictions of Kᵢ or MIC.

Further caveats, discussed in the manuscript: the catalytic μ-hydroxide and all waters are removed before docking, so bridging frequencies are upper bounds; each enzyme is a single rigid crystal structure; AD4Zn binding-free-energy errors (~2 kcal mol⁻¹) exceed most differences between the compounds studied, so no ranking among them should be inferred; and the number of bridging poses out of nine describes the output modes of a stochastic search run once per pair, not a Boltzmann-weighted probability.

---

## Citation

If you use this code, please cite the manuscript and the archived release:

```
[manuscript citation once published]
[Zenodo DOI]
```

Please also cite the underlying tools: AutoDock4Zn (Santos-Martins et al., *J Chem Inf Model* 2014, 54:2371), AutoDock Vina 1.2 (Eberhardt et al., *J Chem Inf Model* 2021, 61:3891), AutoDockTools (Morris et al., *J Comput Chem* 2009, 30:2785) and Open Babel (O'Boyle et al., *J Cheminform* 2011, 3:33).

## License

MIT.
