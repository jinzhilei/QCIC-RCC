# MDD-QCIC: Multi-source Data-Driven Quantitative Cancer-Immunity Cycle Model for Advanced RCC

This repository contains the MATLAB source code accompanying the paper:

> Lei Du, Chenghang Li, Jinzhi Lei. *A multi-source data-driven quantitative cancer-immunity cycle model predicts clinical responses to combination therapy in advanced renal cell carcinoma.* Preprint (2026).

The code implements the **MDD-QCIC** framework, which integrates immunogenomic data,
population pharmacokinetics (PopPK), RECIST v1.1 response criteria, and overall-survival
records to generate a virtual patient cohort and perform *in silico* clinical trials for
atezolizumab + bevacizumab combination therapy in advanced renal cell carcinoma (RCC).

---

## Repository structure

```
RCC/
├── QCIC/                     # Baseline QCIC ordinary differential equation model
│   ├── main_basic.m          #   Entry: integrates the ODE system with ode45
│   ├── QCIC_basic.m          #   RHS of the 28-variable QCIC ODE system
│   ├── parameter_basic.m     #   Fixed model parameters (migration, volume, chemotaxis, ...)
│   └── initialzation_basic.m #   Initial cell densities across the 5 compartments
│
├── VP_Immunogenomic/         # Virtual patient (VP) cohort generation & calibration
│   ├── main_VP.m             #   Entry: LHS sampling → parallel ODE simulation → selection
│   ├── main_parallel.m       #   Per-patient ODE simulation worker (called via parfor)
│   ├── QCIC_VP.m             #   VP-version of the QCIC ODE (heterogeneous parameters)
│   ├── parameter_VP.m        #   VP model parameters
│   ├── initialzation_VP.m    #   VP initial conditions (includes seeded tumor cells)
│   ├── LHS.m                 #   Latin Hypercube Sampling of plausible patients
│   ├── Data_processing.m     #   Loads Immunogenomic.xlsx, log-transforms clinical ratios
│   ├── select.m              #   Probabilistic inclusion filtering (Allen et al., 2016)
│   ├── OptimizeVPgeneration.m#   Goodness-of-fit via simulated annealing (β optimization)
│   ├── nsphereVolume.m       #   N-sphere volume for kNN density estimation
│   ├── KL_JS_div.m           #   KL / JS divergence between VP and clinical distributions
│   ├── picture_select.m      #   Plots the VP vs. clinical immune-ratio distributions
│   ├── ColorMatrix.m         #   Plot color palette
│   └── Immunogenomic.xlsx    #   Clinical immune-subset ratio data (iAtlas KIRC/KIRP)
│
└── picture/                  # Figure generation scripts (Figs. 2–8 of the paper)
    ├── fig_2.m ... fig_8.m   #   One script per manuscript figure
    ├── ColorMatrix.m         #   Color palette
    ├── Color_RdYlBu.m        #   Red-Yellow-Blue colormap
    ├── data/                 #   Pre-computed .mat result files (input to fig_*.m)
    └── picture/              #   Rendered output figures (PDF/PNG)
```

---

## Model overview

The MDD-QCIC model spans **five physiological compartments**: bone marrow & thymus (A),
peripheral blood (B), tumor-draining lymph nodes (C), lymphatic vessels (E), and the tumor
microenvironment (D). It tracks **28 state variables** — dendritic cells, naïve/helper/regulatory/
cytotoxic T cells, macrophages, tumor-associated macrophages, and four tumor-cell subtypes
(drug-sensitive, atezolizumab-resistant, bevacizumab-resistant, double-resistant) — plus the
spatio-temporal dynamics of cytokines and chemokines.

The temporal evolution of each cell type is governed by the balance of eight biological
processes (source, differentiation, proliferation, transition, migration, chemotaxis, killing,
apoptosis), following the general equation:

```
d[i]_n / dt = S_i + D_i + P_i + T_i + V_i + X_i + I_i + A_i
```

Pharmacological mechanisms of **atezolizumab** (PD-L1 blockade enhancing T-cell killing) and
**bevacizumab** (anti-angiogenic reduction of tumor-bearing capacity) is explicitly embedded.

---

## Requirements

- **MATLAB** (R2018b or later recommended)
- **Toolboxes:** Parallel Computing Toolbox (`parpool`, `parfor`), Optimization Toolbox
  (`simulannealbnd`, `optimoptions`), Statistics and Machine Learning Toolbox
  (`kstest2`, `ksdensity`, `readtable`)
- The `QCIC/` and `VP_Immunogenomic/` directories must both be on the MATLAB path
  (the VP pipeline calls `QCIC_VP`, `parameter_VP`, `initialzation_VP`).

---

## How to run

### 1. Baseline QCIC model (sanity check)
```matlab
cd QCIC
main_basic        % integrates the 28-variable ODE for t = 0..365 days
```
Outputs the time series `x` of all cell/molecule populations.

### 2. Generate the virtual patient cohort
```matlab
cd VP_Immunogenomic
main_VP           % ~100,000 plausible patients → calibrated cohort (~1,664 VPs)
```
This script:
1. Draws `sample = 100000` plausible patients via LHS (`LHS.m`) over 7 heterogeneity
   parameters (T-cell / dendritic-cell / macrophage source & proliferation rates).
2. Simulates each patient in parallel (`main_parallel.m` → `QCIC_VP.m`, `ode45`).
3. Loads clinical immune ratios from `Immunogenomic.xlsx` (`Data_processing.m`).
4. Filters patients by probabilistic inclusion against the observed distributions
   (`select.m`), optimizing the scaling factor β with simulated annealing
   (`OptimizeVPgeneration.m`).
5. Computes KL/JS divergence (`KL_JS_div.m`) and plots the matched distributions
   (`picture_select.m`).

> **Note:** Generating 100,000 patients is computationally intensive and requires the
> Parallel Computing Toolbox. Intermediate results are written to `./VP_parallel/` and
> `./VP_immune/`. Adjust `sample` and `parpool(20)` in `main_VP.m` to your hardware.

The selected cohort indices are saved to `./VP_immune/Number_VP_patients.mat`.

### 3. Reproduce the manuscript figures
```matlab
cd picture
fig_2   % VP vs. clinical immunogenomic distributions (paper Fig. 2)
fig_3   % PopPK profiles of atezolizumab / bevacizumab (paper Fig. 3)
fig_4   % In silico RECIST response validation (paper Fig. 4)
fig_5   % Spider / waterfall tumor trajectories (paper Fig. 5)
fig_6   % Alternative dosing strategies (paper Fig. 6)
fig_7   % Immune biomarker ROC / AUC analysis (paper Fig. 7)
fig_8   % Overall-survival Kaplan-Meier stratification (paper Fig. 8)
```
Each `fig_*.m` reads pre-computed result files from `picture/data/` and writes figures to
`picture/picture/`. PopPK, in-silico trial, dosing, ROC, and survival results are expected
as `.mat` inputs under the corresponding `picture/data/fig.N/` subfolders.

---

## Data note

`VP_Immunogenomic/Immunogenomic.xlsx` contains the three clinical immune-subset ratios
(CD8/CD4, CD4/Treg, Treg/TAM) derived from the "Immune Feature Trend" module of the
[iAtlas portal](https://www.itolonc.org/) for clear-cell (KIRC) and papillary (KIRP) RCC.
The intermediate/result `.mat` files produced by the VP pipeline and consumed by the figure
scripts are **not** included in this repository due to size; regenerate them by running the
pipeline above (or use your own calibrated cohort).

---

## Citation

If you use this code, please cite:

```
Du L, Li C, Lei J. A multi-source data-driven quantitative cancer-immunity cycle
model predicts clinical responses to combination therapy in advanced renal cell
carcinoma. Preprint(2026).
```

## License

For academic/research use. Please contact the corresponding author
(Jinzhi Lei, jzlei@tiangong.edu.cn) for permission regarding other uses.
