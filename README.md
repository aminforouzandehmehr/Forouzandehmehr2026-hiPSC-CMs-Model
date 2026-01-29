# Forouzandehmehr2026-hiPSC-CMs-Model-hiMCES
A computational model of human iPSC cardiomyocytes electro-mechano-energetics with a new oxygen consumption component

This repository contains code associated with:

Li W, Forouzandehmehr MA, McLeod D, Pozo MR, Heinson YW, Kay MW, Li Z, Morotti S, Entcheva E.
Response of human iPSC-cardiomyocytes to adrenergic drugs assessed by high-throughput pericellular oxygen measurements and computational modeling.
bioRxiv (2025), https://doi.org/10.1101/2025.06.27.662066.

### Key scripts

- `hiMCES2025.m` — Model ODEs  
- `MasterCompute_hiMCES.m` — Integrates the model with standard plots
- `passiveForces.m` — The function calculating the passive forces for the myofilament
- `fpca.m` — Simulates Tension-Ca response calling myofilament.m
- `myofilament.m` — The myofilament component of the model
- `extractBiomarkers1.m` — Calculates AP and contractile biomarkers
- `extractBiomarkers2.m` — Calculates CaT biomarkers

### MAT data

- `spnt_ss_control_yfinal.mat` — Control-model steady-state initial conditions from spontaneous pacing 
- `1Hz_ss_control_yfinal.mat` — Control-model steady-state initial conditions from 1 Hz pacing
- `ISO_ss_spnt_yfinal.mat` — Isoproterenol steady-state initial conditions from spontaneous pacing
- `ISO_ss_1Hz_yfinal.mat` — Isoproterenol steady-state initial conditions from 1 Hz pacing
- `BLEB_ss_spnt_yfinal.mat` — Blebbistatin steady-state initial conditions from spontaneous pacing
- `BLEB_ss_1Hz_yfinal.mat` — Blebbistatin steady-state initial conditions from 1 Hz pacing
