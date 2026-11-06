# Graphene under true uniaxial stress: nonlinear elasticity, internal relaxation and M-point band gap (DFT)

Data and code for the article

> M. Taleb, M. Attia, C. Derbal, F. Djefaflia, *Orientation-resolved nonlinear elasticity and internal relaxation of graphene under uniaxial stress: does strain open a band gap before the lattice fails?* (submitted).

All numbers, tables and figures in the article can be regenerated from the data in this repository with the scripts in `analysis/` (no DFT run needed). The DFT drivers in `dft/` regenerate the data themselves.

## Contents

```
dft/
  driver.py            uniaxial-stress driver: for each applied strain eps along n(theta), iterates the
                       transverse strain et and shear g (decoupled secant) until sigma_tt = sigma_nt = 0,
                       with full atomic relaxation (Quantum ESPRESSO pw.x)
  trimgap.py           scf + nscf at the three M points + projwfc for each relaxed structure
  pz.py                parser of projwfc output: band energies and C 2p_z weights
  run/th{0,10,20,30}/results.json   relaxed structures, strains and stresses (one entry per eps)
  gap/trim_th{0,10,20,30}.json      M-point pi-band splittings and sigma states near E_F
  pseudo/              place C.UPF here (see below)
analysis/
  analyze.py           loads results, nominal stress, relaxed and affine bond strains, TB helper
  final.py             ideal strain/strength, global symmetry-constrained fits -> summary.json
  gapfit.py            signed M-point gap, gap-opening strain, TB (t0, beta) fit -> gapsummary.json
  figs.py              Figures 1-5 of the article -> analysis/fig/
  ga.py                graphical abstract -> analysis/out/
  summary.json, gapsummary.json     derived quantities used in the article
examples/scf_example.in  example pw.x input (2% zigzag strain)
```

### Data format

`dft/run/th<theta>/results.json` is a list; each entry has

| key | meaning | unit |
|---|---|---|
| `theta` | loading angle from the zigzag axis | deg |
| `eps` | applied engineering strain along **n** | – |
| `et`, `g` | transverse strain and tensorial shear strain in the stretch **U** | – |
| `E` | total energy | eV |
| `snn`, `stt`, `snt` | 2D stresses in the loading frame (tension positive) | N/m |
| `cell` | deformed lattice vectors (rows) | Å |
| `frac`, `pos` | fractional and Cartesian atomic positions | –, Å |
| `iters` | secant iterations needed | – |

`dft/gap/trim_th<theta>.json` maps each strain (key `"0.240"` etc.) to the Fermi energy, the pi-band splitting at M1 = b1/2, M2 = b2/2, M3 = (b1+b2)/2 (unsigned, eV), and the sigma-type states (p_z weight < 0.5) within 3 eV of E_F (energy relative to E_F, eV).

## Reproducing the analysis

Requirements: Python ≥ 3.9 with `numpy`, `scipy`, `matplotlib` (`pip install -r requirements.txt`).

```bash
cd analysis
python final.py     # writes summary.json
python gapfit.py    # writes gapsummary.json, prints fitted t0, beta and gap-opening strains
python figs.py      # writes fig/fig1..fig5 (.pdf, .png)
python ga.py        # writes out/Graphical_abstract.*
```

## Re-running the DFT calculations

* Quantum ESPRESSO 6.6 (`pw.x`, `projwfc.x`); set `PW` and `PROJWFC` if they are not on the `PATH`.
* Pseudopotential: SSSP (efficiency) carbon, `C.pbe-n-kjpaw_psl.1.0.0.UPF` (PSlibrary 1.0.0), saved as `dft/pseudo/C.UPF` (or set `PSEUDO_DIR`).
* Settings (defaults in `driver.py`): PBE, 50/400 Ry, 24×24×1 Monkhorst–Pack, Marzari–Vanderbilt smearing 0.01 Ry, vacuum Lz = 12 Å, a0 = 2.4669 Å; forces < 5e-5 Ry/bohr; |sigma_tt|, |sigma_nt| < 0.02 N/m.

```bash
cd dft
python driver.py 10 0.005,0.01,0.02,0.03,0.04,0.06,0.08,0.10,0.12,0.15,0.18,0.21,0.24,0.27,0.30
python trimgap.py 10 0.12,0.15,0.18,0.21,0.24,0.27
```
`driver.py` resumes from an existing `results.json`, so interrupted runs can be restarted with the same command.

## License

Code: MIT License (see `LICENSE`). Data (`*.json`): CC BY 4.0.
