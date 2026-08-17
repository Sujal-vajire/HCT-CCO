# HCT-CCO: Hierarchy CT-Driven Constrained Constructive Optimization

HCT-CCO generates subject-specific 3D coronary arterial trees by growing them as
anatomically grounded extensions of a segmented epicardial **base tree**, confined to a
CT-derived myocardial domain. It is the reference implementation accompanying the HCT-CCO
article.

A single command generates the tree **and** runs the full post-processing, so one run
produces the final artifacts directly (3D meshes, morphometry, and a hemodynamic check),
with no intermediate steps to run by hand.

---

## How to cite

If you use this code, please cite the HCT-CCO article:

> S. L. Vajire, J. S. Choy, G. S. Kassab and L.-C. Lee, "Hierarchy CT-Driven Constrained
> Constructive Optimization: A Framework for Coronary Network Generation." *(citation to be
> updated upon publication.)*

This work builds on the OpenCCO implementation:

> B. Kerautret, P. Ngo, N. Passat, H. Talbot and C. Jaquet, "OpenCCO: An Implementation of
> Constrained Constructive Optimization for Generating 2D and 3D Vascular Trees,"
> *Image Processing On Line*, vol. 13, pp. 258–279, 2023.

---

## Repository layout

```
HCT-CCO/
├── CMakeLists.txt
├── LICENSE.md
├── bin/
│   ├── CMakeLists.txt
│   └── HctCco.cpp        # command-line generator (3D)
├── src/
│   ├── HctccoModel.h/.ih  # tree model, hydraulics, Kamiya bifurcation optimization
│   ├── MyocardialDomain.h        # myocardial-domain / distance-field constraints
│   └── helpers/
│       ├── CoronaryGeometry.h         # geometry and vessel-intersection tests
│       ├── CoronaryGrowth.h   # growth loop and candidate sampling
│       └── CoronaryTreeIO.h          # base-tree import and XML export
├── ext/
│   └── CLI11.hpp                 # bundled command-line parser
├── samples/
│   ├── case0.vol                 # sample myocardial domain (mask volume)
│   └── case0.xml                 # sample segmented epicardial base tree (GXL)
└── postprocessing/
    └── coronary_postprocess.py            # ordering, morphometry, hemodynamics, VTP/OBJ export
```

---

## Requirements

**Generator (C++):**
- CMake ≥ 3.5 and a C++11 compiler
- [DGtal](https://dgtal.org/) (`find_package(DGtal)`)
- [Ceres Solver](http://ceres-solver.org/) 2.x (`find_package(Ceres)`)

**Post-processing (Python 3):** `numpy`, `vtk`, `matplotlib`

The build records the path to `postprocessing/coronary_postprocess.py`, and the generator invokes
it automatically after writing the tree. The generator calls `python3` from your `PATH`,
so make sure that interpreter has `numpy`, `vtk`, and `matplotlib` installed (or pass
`--no-postprocess` and run the script yourself).

---

## Build

```bash
mkdir build && cd build
cmake ..
make
# -> build/bin/HctCco
```

---

## Run the sample case

The sample grows a coronary tree from the segmented epicardial **base tree**
(`case0.xml`) into the **myocardial domain** (`case0.vol`). Using a subject-specific base
tree is the central feature of HCT-CCO, so both inputs are provided together:

```bash
./build/bin/HctCco \
    -d samples/case0.vol \
    -i samples/case0.xml \
    -n 40000 \
    -a 10000 \
    -m 1 \
    -p -26 41 13 \
    -x final_tree.xml \
    -w
```

The post-processing runs automatically and writes all outputs next to `final_tree.xml`.

> **Runtime:** `-n 40000` is a full-scale run and takes on the order of hours (runtime grows
> super-linearly with the terminal count). For a quick end-to-end check, use a small count
> such as `-n 500`, which completes in seconds.

### Outputs (prefix = export-XML basename, here `final_tree`)

| File | Contents |
|---|---|
| `final_tree.xml` | generated tree graph (GXL) with per-segment radius and flow |
| `final_tree_tree.vtp` | tree for ParaView, with Strahler `final_order` and `radius` |
| `final_tree_tree.obj` | 3D tube mesh of the tree |
| `final_tree_diameter_vs_order.png` | vessel diameter vs Strahler order, overlaid on literature |
| `final_tree_hemodynamics.png` | wall shear stress and pressure vs diameter (resting + stress) |
| `final_tree_hemodynamics.csv` | per-order shear, pressure, and flow |
| `final_tree_connectivity_matrix.csv` | order-to-order connectivity |
| `final_tree_diameter_stats.csv` / `..._diameter_ranges.csv` | per-order diameter statistics |

The transmural wall-distribution analysis (`-w`) reports the inner-50% vs outer-50% split
of terminal vessels.

To generate the raw tree only, add `--no-postprocess`. To run post-processing on an
existing tree XML directly:

```bash
python3 postprocessing/coronary_postprocess.py final_tree.xml final_tree samples/case0.vol
```

---

## Command-line options

| Option | Meaning |
|---|---|
| `-d, --organDomain` | myocardial-domain mask volume (`case0.vol`) |
| `-i, --importXML` | segmented epicardial base tree; growth starts from it and extends beyond it (`case0.xml`) |
| `-n, --nbTerm` | number of terminal vessels |
| `-a, --aPerf` | perfusion volume |
| `-m, --minDistanceToBorder` | clearance from the wall (use `1` for `case0`) |
| `-p, --posInit` | initial root position (x y z); default is the domain center |
| `-g, --gamma` | bifurcation exponent (default 3.0) |
| `--rootRadiusMax` | cap on root radius in mm (default 1.588) |
| `-x, --exportXML` | export tree XML (also the post-processing input) |
| `-w, --wallAnalysis` | report transmural (inner/outer 50%) distribution of terminals |
| `--no-postprocess` | generate the raw tree only, skip post-processing |

Run `HctCco --help` for the complete list.

---

## What the pipeline does

1. **Generation** — grows the tree from the base tree into the myocardial domain, adding one
   terminal at a time with flow-based bifurcation optimization and hierarchy-aware
   branching constraints.
2. **Ordering** — assigns each vessel a diameter-based Strahler order.
3. **Morphometry** — produces the diameter-vs-order curve against literature, the
   order-to-order connectivity matrix, and per-order diameter statistics.
4. **Hemodynamics** — a steady lumped-parameter resistive-network solve reports wall shear
   stress and pressure at both resting and stress perfusion.

---

## Notes

- `samples/case0` requires `-m 1`; the default clearance is too restrictive for this domain.
- Runtime grows super-linearly with `-n`; large runs take hours, a few hundred terminals
  complete in seconds.
- Make sure the `python3` on your `PATH` is the one with `numpy`, `vtk`, and `matplotlib`.

---

## Acknowledgments

This work was supported by NIH R01HL160997 and NSF 2222066.

---

## License and attribution

This project is a derivative work of **OpenCCO** and is distributed under the **GNU
General Public License v3** (see `LICENSE.md`, which is unchanged from the original).

- **Base code (OpenCCO):** Copyright (C) 2023 B. Kerautret, P. Ngo, N. Passat, H. Talbot
  and C. Jaquet. Source: https://github.com/OpenCCO-team/OpenCCO
- **Coronary-specific extension (HCT-CCO):** Copyright (C) 2024–2026 S. L. Vajire,
  J. S. Choy, G. S. Kassab and L.-C. Lee (Michigan State University; California Medical
  Innovations Institute).

Modifications made to the original OpenCCO code, released under the same GPL v3:

- Growth resumed from a segmented patient-specific epicardial base tree and a CT-derived
  myocardial domain, instead of a synthetic root.
- Boundary-aware growth confined to the myocardial wall via a Euclidean distance field.
- Flow-weighted bifurcation initialization with HK-type daughter-radius seeding for the
  Kamiya-type local optimization.
- Hierarchy-aware, stage-wise branching schedule controlling where new vessels attach.
- Additional physiological constraints: root-radius cap, taper floor, HK/Murray radius
  scaling, minimum radius, maximum length, and a radius-aware vessel-overlap test.
- Incremental (path-only) hydraulic resistance and flow updates.
- New in-house post-processing (`postprocessing/coronary_postprocess.py`): diameter-based Strahler
  ordering, order–order connectivity, diameter–order morphometry, a lumped-parameter
  hemodynamic check, and VTP/OBJ export.
