# Perovskite Polariton Analysis

MATLAB code used for coupled-oscillator modelling of the
angle-resolved polariton dispersion reported in the manuscript
"Stable Continuous-Wave Perovskite Lasers by Harvesting Triplet Excitons".

## Requirements

- MATLAB R20XXx
- Operating system: Windows 10/11
- No non-standard hardware is required.

## Description

The code fits the experimentally measured angle-resolved dispersion
using a coupled-oscillator model and calculates the lower and upper
polariton branches, Rabi splitting, and Hopfield coefficients.

## Demo data

Example angle-resolved experimental data are provided in the
`data` folder.

## Usage

1. Open MATLAB.
2. Set the repository folder as the current working directory.
3. Run `fit_polariton_lp_cavity.m`.
4. The calculated LP and UP dispersions and fitting parameters will
   be generated.

## Expected output

- Lower-polariton (LP) dispersion
- Upper-polariton (UP) dispersion
- Cavity dispersion
- Rabi splitting
- Polariton dispersion figure

## Runtime

Typical runtime: <1 min on a standard desktop computer.