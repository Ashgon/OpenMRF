![OpenMRF](OpenMRF_banner.png)

# OpenMRF

## Introduction

OpenMRF is a modular and vendor-neutral framework for Magnetic Resonance Fingerprinting (MRF) built on the open-source [Pulseq](https://pulseq.github.io) standard. It is built upon the MATLAB version of Pulseq by Kelvin J. Layton and Maxim Zaitsev ([doi:10.1002/mrm.26235](https://doi.org/10.1002/mrm.26235)). OpenMRF unifies all core components of the MRF workflow within a single MATLAB-based toolbox: flexible sequence generation, automated Bloch-based dictionary simulation, and low-rank image reconstruction. The provided tools support a wide range of contrast preparations and readouts (e.g., spiral, radial, rosette) and include integrated solutions for trajectory calibration, spin-lock modeling, slice profile simulation, and metadata storage. Designed for reproducibility and portability, OpenMRF enables easy deployment of MRF protocols across multiple scanner platforms, including Siemens and GE systems.

## Contents

- `include_miitt/`: Contains the low-rank reconstruction code provided by the MIITT group and Jeffrey Fessler's [MIRT toolbox](https://web.eecs.umich.edu/~fessler/code/). Includes an installation script. **Do not** add this folder manually to your MATLAB path; use the `install_OpenMRF.m` script.
- `include_misc/`: Miscellaneous utilities and helper functions.
- `include_pulseq/`: Copy of the official Pulseq repo ([GitHub link](https://github.com/pulseq/pulseq), v1.5, 01.04.2025). Includes minor modifications to the plotting functions for visualizing trigger inputs/outputs.
- `include_pulseq_toolbox/`: Contains standard imaging readouts (cartesian, radial, spiral, rosette) combined with various preparation modules (inversion recovery, saturation, MLEV-T2, spin-lock, adiabatic spin-lock, CEST). Also includes simulation tools for MRF dictionary generation.
- `include_pre_sim_library/`: Library of optimized RF pulse waveforms (including `sigpy`-generated pulses) and flip angle patterns for MRF. Also used to store pre-simulated slice profiles, adiabatic efficiencies and compressed dictionaries.
- `pulseq_sequences/`: Example Pulseq sequences and reconstruction scripts.

## System Requirements

- **MATLAB** tested with R2024b on Win11 and Ubuntu 22.04
- **Python** with the `sigpy` package — required for designing SLR and adiabatic RF pulses (e.g., BIR-4)

## Getting Started

1. **Clone the repository**:
   ```bash
   git clone https://https://github.com/HarmonizedMRI/OpenMRF
   ```

2. **Create your pulseq backup directory** (e.g., `C:/Users/YourName/Pulseq`).  
   This will be used to:
   - Store `.mat` files containing k-space trajectories, sequence parameters and header information
   - Save a copy of toolbox functions for reproducibility

3. **Start MATLAB** and run:
   ```matlab
   install_OpenMRF.m
   ```

   During setup:
   - Enter a Pulseq username (single word, no whitespace)
   - Enter your lab name (e.g., University of Wuerzburg, Department of Physics, EP5, Wuerzburg, Germany)
   - Choose your Pulseq path (backup directory; e.g., `C:/Users/YourName/Pulseq`)
   - Select your MRI system (or define a new one using a `.csv` file in `include_pulseq_toolbox/system_specifications`)

4. **Generate an example sequence**:
   - Navigate to `pulseq_sequences/fingerprinting/example_mrf_sequences/`
   - Run `pulseq_mrf.m` with `flag_backup = 1`
   - Check your backup folder for:
     - The generated `.seq` file
     - A `.mat` file with all relevant parameters (trajectories, timings, etc.)

## Note

This repository includes third-party software distributed under their respective licenses. Please consult the [NOTICE](./NOTICE) file before use. Not all code is MIT-licensed. Users intending commercial use must review the [NOTICE](./NOTICE) file and seek appropriate permissions from the original authors.

_Maximilian Gram: 07.08.2025_