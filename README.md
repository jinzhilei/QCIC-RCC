# MDD-QCIC: Multi-source Data-Driven Quantitative Cancer-Immunity Cycle Model for Advanced RCC

This repository contains the MATLAB source code accompanying the paper:

> Lei Du, Chenghang Li, Jinzhi Lei. *A multi-source data-driven quantitative cancer-immunity cycle model predicts clinical responses to combination therapy in advanced renal cell carcinoma.* J. R. Soc. Interface (2026).

The code implements the **MDD-QCIC** framework, which integrates immunogenomic data,
population pharmacokinetics (PopPK), RECIST v1.1 response criteria, and overall-survival
records to generate a virtual patient cohort and perform *in silico* clinical trials for
atezolizumab + bevacizumab combination therapy in advanced renal cell carcinoma (RCC).

---

## Repository structure

```
RCC/
├── QCIC/                     # Baseline QCIC ordinary differential equation model (sanity check)
│   ├── main_basic.m          #   Entry: integrates the 28-variable ODE with ode45
│   ├── QCIC_basic.m          #   RHS of the baseline QCIC ODE system
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
├── Pharmacokinetics/         # Population PK: two-compartment model for both drugs
│   ├── PK_model_atezolizumab/
│   │   ├── main_atezolizumab.m          #   Entry: Beta sampling + Q3W dosing simulation
│   │   ├── PK_model_atezolizumab.m      #   Two-compartment ODE (central/peripheral)
│   │   ├── PK_model_atezolizumab_parameter.m # PK-parameter Beta sampling
│   │   ├── Beta_distribution.m          #   Beta-distributed random draws
│   │   └── Numerical_integration.m      #   AUC by trapezoidal rule
│   └── PK_model_bevacizumab/            #   (analogous files for bevacizumab)
│
├── RECIST/                   # In silico clinical trial under 3 treatment arms
│   ├── 2. NoDrug/            #   Control arm (no therapy)
│   │   ├── main_VP.m                  #   Entry: simulate cohort & classify by RECIST
│   │   ├── QCIC_VP_NoDrug.m           #   ODE (no drug terms)
│   │   ├── parameter_VP_NoDrug.m / initialzation_VP_NoDrug.m
│   │   ├── solve_one_patient.m        #   Single-patient integrator
│   │   ├── tumor_parameter.m          #   Tumor-heterogeneity parameter sampling
│   │   ├── PSO.m / e.m                #   PSO calibration to RECIST data
│   │   ├── model_parameter_select.m   #   Immune-parameter sampling
│   │   └── evaluating.m               #   RECIST counting (CR / PR / SD / PD)
│   ├── 3. A/                 #   Atezolizumab monotherapy arm (analogous scripts)
│   └── 4. AB/                #   Atezolizumab + bevacizumab combination arm
│       ├── main_VP_AB.m              #   Entry for the combination arm
│       ├── QCIC_VP_AB.m              #   ODE with both drugs (PK-coupled killing/angiogenesis)
│       └── ...                        #   (same helper scripts as above)
│
├── Treat/                    # Alternative dosing strategies & recurrence-rate analysis
│   ├── treat_AB.m            #   Entry: simulate 3 schedules × dose ratios
│   ├── QCIC_VP_AB.m          #   Combination ODE (shared with RECIST/4. AB)
│   ├── solve_one_patient.m   #   Single-patient sim with impulse dosing
│   ├── tumor_parameter.m     #   Tumor sub-population Beta sampling
│   ├── parameter_VP_AB.m / initialzation_VP_AB.m
│   ├── model_parameter_select.m # reads VP_Immunogenomic/VP_immune/PDF_immune.mat
│   ├── evaluating.m          #   Evaluates response metrics per schedule
│   ├── counts_recurrence_rate.m   # aggregates recurrence / response rates
│   ├── Recurrence_rate/      #   Output .mat files per schedule & dose pair
│   └── VP_immune/            #   Intermediate immune / tumor initialization data
│
├── OS/                       # Overall-survival analysis & calibration
│   ├── OS.m                  #   Entry: survival-index fit to clinical OS curves
│   ├── OS_P.m                #   OS stratified by immune biomarkers
│   ├── PSO.m                 #   PSO optimizer for survival-index params (μ, σ)
│   └── data/                 #   Clinical OS refs (Rini 2019, Rini 2021, McGregor 2020)
│
├── ROC/                      # ROC / AUC immune-biomarker analysis (paper Fig. 7)
│   └── ROC.m                 #   Computes responder labels & immune profiles → VP_immune/
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
**bevacizumab** (anti-angiogenic reduction of tumor carrying capacity) are embedded explicitly.

---

## Requirements

- **MATLAB** (R2018b or later recommended)
- **Toolboxes:** Parallel Computing Toolbox (`parpool`, `parfor`), Optimization Toolbox
  (`simulannealbnd`, `optimoptions`), Statistics and Machine Learning Toolbox
  (`kstest2`, `ksdensity`, `readtable`, `betarnd`, `cumtrapz`)
- The `QCIC/`, `VP_Immunogenomic/`, and `RECIST/` directories must be on the MATLAB path
  (the downstream pipelines call `QCIC_VP*`, `parameter_VP*`, `initialzation_VP*`,
  and the shared `Beta_distribution` / `solve_one_patient` helpers).
- The modules read/write intermediate `.mat` files through **relative paths**, so run each
  entry script from *inside its own folder* (e.g. `cd RECIST/'4. AB'`) and execute the
  pipeline in the order below (Step 2 must finish before Steps 4–6).

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

### 3. Population pharmacokinetics (PopPK)
```matlab
cd Pharmacokinetics/PK_model_atezolizumab
main_atezolizumab      % 1000 VPs, Q3W 1200 mg; two-compartment ODE → ./VP/

cd ../PK_model_bevacizumab
main_bevacizumab       % analogous run for 1050 mg bevacizumab
```
Each `main_*.m` draws PK parameters from Beta distributions (`PK_model_*_parameter.m`),
simulates the central/peripheral compartment ODE with impulse (bolus) dosing every 21 days,
computes the 95% prediction interval and AUC, and saves results under `./VP/`.
(These profiles feed the drug terms in the RECIST and Treat simulations.)

### 4. In silico clinical trial (RECIST v1.1)
```matlab
cd RECIST/'2. NoDrug';  main_VP       % control arm
cd ../'3. A';           main_VP_A     % atezolizumab monotherapy
cd ../'4. AB';          main_VP_AB    % atezolizumab + bevacizumab
```
Each arm calibrates tumor-heterogeneity parameters to clinical RECIST data via particle
swarm optimization (`PSO.m` / `e.m`), simulates every virtual patient with impulse dosing,
and classifies the outcome into CR / PR / SD / PD (`evaluating.m`). The combination arm
(`4. AB`) writes per-patient trajectories to `RECIST/4. AB/VP_AB/VPsample_i.mat`, which
are later consumed by the OS and Treat modules.

### 5. Alternative dosing strategies & recurrence rate (Treat)
```matlab
cd Treat
treat_AB                % 3 schedules (maximum / average / waning) × 4 dose ratios
counts_recurrence_rate  % aggregate recurrence / response rates across all arms
evaluating              % response metrics per schedule
```
`treat_AB.m` reuses the combination ODE (`QCIC_VP_AB.m`) with impulse dosing defined by
`programme` (dosing interval) and `ratio` (fractionation scheme), sampling patient-specific
tumor parameters via `tumor_parameter.m` and immune parameters via `model_parameter_select.m`
(reads `VP_Immunogenomic/VP_immune/PDF_immune.mat`). Outputs go to `Treat/Recurrence_rate/`.

### 6. Overall-survival analysis (OS)
```matlab
cd OS
OS          % fit survival index χ to clinical OS (Rini 2019); writes OS.mat
OS_P        % biomarker-stratified OS (Tc, Tr/Tc, Tc/TAM, Th/Tc); writes OS*.mat
```
`OS.m` maps the tumor-burden survival index to a logistic mortality probability and fits
`μ, σ` with the custom PSO optimizer (`OS/PSO.m`), using the clinical OS curves in `OS/data/`
(`1.Rini.2019.dat`, `2.Rini.2021.dat`, `3.McGregor.2020.dat`). `OS_P.m` reads the
combination-arm trajectories from `RECIST/4. AB/VP_AB/` to stratify patients.

### 7. Immune-biomarker ROC / AUC analysis (ROC)
```matlab
cd ROC
ROC          % responder (PD vs. others) labels + immune profiles at week 12
```
`ROC.m` reads the combination-arm trajectories from `RECIST/4. AB/VP_AB/` (note: the script
references this path as `../RECISE/...` — a typo for `RECIST`), derives the Tumor Response
Index at day 84, labels responders/non-responders, and saves `response.mat` / `immune.mat`
to `ROC/VP_immune/`. These feed the ROC curves and AUC values in paper Fig. 7.

### 8. Reproduce the manuscript figures
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
carcinoma. J R Soc Interface. 2026.
```

## License

For academic/research use. Please contact the corresponding author
(Jinzhi Lei, jzlei@tiangong.edu.cn) for permission regarding other uses.
