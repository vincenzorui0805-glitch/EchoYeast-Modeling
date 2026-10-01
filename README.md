# EchoYeast-Modeling

Multi-scale computational pipeline for H1N1-specific nanobody screening.

**Team:** EchoYeast | **iGEM 2026**

## Pipeline

ClusPro rigid docking → HADDOCK 2.5 semi-flexible refinement → PyMOL/SPICE interface analysis → GROMACS 100 ns MD → gmx_MMPBSA binding free energy + per-residue decomposition.

## Key Result

**R1aB6** is the top candidate for H1N1 HA:

- ClusPro specificity index: **-93.25**
- HADDOCK score: **-93.2** (BSA 1890.1 Å²)
- 100 ns MD: stable (RMSD 0.49 nm, interface RMSF < 0.33 nm)
- MM-GBSA ΔG_bind: **-14.67 kcal/mol**
- Hotspot residue: **ARG-B:31 (-5.46 kcal/mol)**

Positive control: 9VH4 redocking, i-RMSD = **2.113 Å**.

## Directory

| Folder | Contents |
| :--- | :--- |
| `mdp_files/` | GROMACS MDP parameter files |
| `gromacs/` | 100 ns MD outputs and analysis curves |
| `mmgbsa/` | MM-GBSA energy and per-residue decomposition |
| `scripts/` | Python analysis scripts |
| `figures/` | Figures used in the iGEM Wiki |
| `positive_control/` | 9VH4 redocking validation |
| `data/` | Sequences, docking scores, contact maps |
| `structures/` | Complex PDB + 3Dmol interactive HTML |

## Reproducibility

- GROMACS 2025.2, Amber99sb-ildn, TIP3P, gen_seed = -1
- gmx_MMPBSA 1.5.0.3 (AmberTools 20), 51 frames (every 2 ns)
- HADDOCK 2.5 web server, default semi-flexible protocol

## References

See iGEM Wiki Model page for full reference list.
