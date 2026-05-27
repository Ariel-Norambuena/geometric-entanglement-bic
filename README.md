# Geometric control of maximal entanglement via BICs

MATLAB scripts used to reproduce the numerical and analytical results shown
in Figs. 2 and 3 of the manuscript:

**Geometric control of maximal entanglement via bound states in the continuum**

The code studies two giant atoms coupled to a one-dimensional waveguide. It
connects the geometry of the coupling points with the concurrence, fidelity,
and dynamical robustness of bound states in the continuum (BICs).

## Repository contents

- `Codes/Figure2.m` reproduces the concurrence curve and fidelity landscape.
- `Codes/Figure3_Left_Panel.m` reproduces the exact concurrence dynamics and
  Markovian exponential benchmarks.
- `Codes/Figure3_Right_Panel.m` reproduces the on-shell decay-rate maps.

All generated figures are written to `outputs/`.

## Requirements

- MATLAB R2022a or newer recommended.
- No external datasets are required.
- The scripts were prepared and checked in MATLAB R2025b.

The exact dynamics script diagonalizes dense finite Hamiltonians and is the
most computationally demanding part of the repository. The default parameters
match the manuscript settings.

## Reproducing the figures

From the repository root, run:

```matlab
cd Codes
Figure2
Figure3_Left_Panel
Figure3_Right_Panel
```

The scripts export PDF and PNG files to `../outputs`.

For a quick smoke test of the exact dynamics script without the full
publication grid, set:

```matlab
setenv('FIG3_FAST_TEST','1')
Figure3_Left_Panel
```

Unset the variable, or set it to `0`, before regenerating the manuscript
figure.

## Numerical method

For Fig. 3, the non-Markovian memory-kernel equation is solved through its
equivalent finite-mode Hamiltonian representation in the single-excitation
sector. The waveguide is discretized using positive quasi-momenta in
`[0, pi)` and explicit right/left propagation directions, which retain the
resonant pair `+/- k`. The Hamiltonian is diagonalized once for each parameter
set, and the atomic amplitudes are reconstructed from the spectral expansion.
The concurrence is then computed as `C(t) = 2*abs(c1(t)*c2(t))`.

## Contact

For questions about the code or manuscript, please contact Ariel Norambuena.
