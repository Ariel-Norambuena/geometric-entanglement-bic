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
- `scripts/reproduce_previous_paper_baseline.m` writes or verifies a compact
  production-grid regression baseline for the previous Fig. 3 dynamics.
- `src/matlab/` contains the first development modules for the new
  Floquet-controlled Bell-BIC project.
- `docs/` records the numerical conventions, validation status, and scientific
  caveats for the Floquet extension.

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

## Baseline regression

The previous-paper Fig. 3 left-panel dynamics are frozen in
`data/baseline_previous_paper.json`. To verify that later Floquet development
does not silently change the old result, run:

```matlab
run('tests/test_previous_paper_regression.m')
```

This test recomputes selected production-grid observables and checks the
concurrence against the stored baseline.

## Floquet Bell-BIC development

The Floquet extension currently uses a full-Brillouin-zone finite-mode model
with state ordering `psi = [c1, c2, phi_k]`. The implementation keeps the
laboratory-frame driven model, rotating-frame micromotion, central-sideband
effective model, and Markov diagnostics conceptually separate.

Useful development checks:

```matlab
run('tests/test_floquet_controls.m')
run('tests/test_effective_rhs_static_limit.m')
run('tests/test_counterdiabatic_atomic_limit.m')
run('tests/test_lab_frame_cd_micromotion.m')
addpath('scripts')
run_dark_state_passage
scan_passage_duration
run_counterdiabatic_passage
run_two_stage_bell_bic_protocol
scan_k_robustness_comparison
validate_floquet_effective_model
```

The reduced-grid development scan in `data/development/` is not a manuscript
claim. It currently shows that the simple unassisted `u1/u2` passage does not
yet prepare a high-fidelity Bell-BIC; the remaining scientific step is
validating the full instantaneous BIC branch, including photonic dressing.

The counterdiabatic development script adds the minimal atomic correction
`H_CD = thetaDot(t) sigma_y` and exports concurrence, protocol fidelity, and
control plots to `outputs/floquet/`.

The two-stage protocol script first prepares `|eg>` from `|gg>` with a local
single-atom `pi` pulse and then runs the CD-assisted Floquet passage in one
continuous dynamics.

To regenerate the paper-ready two-stage figure, run:

```matlab
addpath('scripts')
run_two_stage_bell_bic_protocol("paper")
make_two_stage_bell_bic_paper_figure
```

The final figure is exported to `figures/Figure_TwoStage_Bell_BIC.pdf`,
`figures/Figure_TwoStage_Bell_BIC.png`, and
`figures/Figure_TwoStage_Bell_BIC.tif`.

To compare the \(K=\pi/2\) sensitivity of a passive PRA-style benchmark with
the two-stage Floquet-CD protocol, run:

```matlab
addpath('scripts')
scan_k_robustness_comparison("paper")
```

This exports `figures/Figure_K_Robustness_Comparison.pdf` and matching
source data under `data/paper/`.

To validate the central-sideband Floquet model against the exact finite-mode
laboratory-frame dynamics, run:

```matlab
addpath('scripts')
validate_floquet_effective_model("paper")
```

This exports `figures/Figure_Floquet_Effective_Validation.pdf` and matching
source data under `data/paper/`.

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
