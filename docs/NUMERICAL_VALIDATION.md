# Numerical Validation Plan

## Completed Baseline Checks

- The previous-paper Fig. 3 left-panel dynamics were restored to the paper
  layout and exact spectral reconstruction.
- `data/baseline_previous_paper.json` stores production-grid observables at
  selected times.
- `tests/test_previous_paper_regression.m` recomputes the baseline and checks
  concurrence against a tolerance of `5e-8`.
- MATLAB `checkcode` reports no messages for the current baseline and Floquet
  development files.

## Completed Floquet Unit Checks

- Bessel inversion on the first branch satisfies `J0(beta_i)=u_i`.
- Endpoint derivatives of the smooth schedule are finite and vanish.
- The two-quadrature control identity is verified by numerical integration.
- The effective-model RHS preserves norm instantaneously for a random state.
- The sideband diagnostic reports dimensionless separation ratios.
- The counterdiabatic atomic-limit test verifies that `H_CD =
  thetaDot(t) sigma_y` drives the closed atomic state exactly along the target
  path.

## Development Effective-Passage Scan

The first reduced-grid effective-model scan used:

```text
Nc = 204
nu/xi = 8
u0 = 0.9
xi*T in {2, 5, 10, 20, 50, 100}
rtol = 1e-9
atol = 1e-11
```

The final conditional Bell-plus fidelity remains close to `0.51` across this
development scan, while norm conservation stays below `4e-15`. This is not a
production result. It is an early warning that unassisted variation of
`u1/u2` does not automatically imply finite-time following of the intended
atomic direction. The next necessary step is to construct and validate the
full instantaneous BIC branch, including photonic dressing, before assigning a
physical interpretation to the poor transfer.

## Additional Validation and Future Refinements

- Production-grid dark-state-passage duration scans.
- Frozen instantaneous BIC branch construction, including photonic dressing.
- Real-space separation of localized photonic dressing from outgoing radiation.
- Solver tolerance convergence.
- Finite-size convergence using grids that retain `k = +/- pi/2`.

These items remain useful for a deeper microscopic characterization of the
instantaneous dressed BIC branch. The paper-level preparation claims below are
based on the CD-assisted two-stage dynamics, the laboratory-frame Floquet
validation, and the robustness scans.

## Development Counterdiabatic Result

The reduced-grid CD comparison used the same parameters as the first
effective-passage scan, with `xi*T = 20`. The unassisted run ends at
`C = 0.018966` and `F_D^(a) = 0.509620`, while the CD-assisted run reaches
`C = 0.995871` and `F_D^(a) = 0.999990`. Norm conservation remains below
`6e-15`, and the peak CD amplitude is
`max[Omega_CD(t)]/xi = 0.073631`.

The corresponding plot is written to:

```text
outputs/floquet/counterdiabatic_passage_dev.png
outputs/floquet/counterdiabatic_passage_dev.pdf
```

This is strong evidence that the failure of the unassisted reduced-grid
passage is a nonadiabatic-following problem in the atomic direction. The
laboratory-frame comparison reported below validates the central-sideband
description used for the CD-assisted passage.

## Development Two-Stage Protocol

The full two-stage diagnostic starts from `|gg,0>`, applies a smooth local
`pi` pulse on atom 1 for `xi*T_pi=2`, and then applies the CD-assisted Floquet
passage for `xi*T=20`. The run is performed as one continuous ODE propagation
in the enlarged basis `{|gg,0>, |eg,0>, |ge,0>, |gg,1_k>}`.

For the same reduced grid, the preparation stage reaches `P_eg=1.000000` at
the stage boundary. The final CD-assisted result is
`C=0.995871`, `F_{Psi+}=0.995881`, and
`F_{Psi+}^{(a)}=0.999990`, with norm error below `3e-13`. The comparison
without CD remains at `C=0.018966`.

The corresponding plot is written to:

```text
outputs/floquet/two_stage_bell_bic_protocol_dev.png
outputs/floquet/two_stage_bell_bic_protocol_dev.pdf
```

## Production K-Mismatch Robustness Comparison

The robustness scan in `scripts/scan_k_robustness_comparison.m` reproduces the
PRA-style passive Markovian benchmark from `Effect_dk_Bell_states.m` and
compares it with the two-stage Floquet-CD protocol. The production run uses:

```text
Nc = 2004
nu/xi = 8
u0 = 0.9
xi*T_pi = 2
xi*T = 20
deltaK/K0 in {-0.10, -0.075, ..., 0.075, 0.10}
```

In the passive benchmark, the fidelity at `xi*t=500` drops from `0.940` at
`deltaK/K0=0.005` to `0.782` at `0.01`, `3.8e-3` at `0.05`, and `1.1e-7` at
`0.10`. In the effective Floquet-CD benchmark, the final Bell fidelity remains
above `0.992` for every tested mismatch with `|deltaK|/K0 <= 0.10`. The same
scan gives final concurrence and atomic population above `0.992` over that
window.

The paper-ready figure and source data are written to:

```text
figures/Figure_K_Robustness_Comparison.pdf
figures/Figure_K_Robustness_Comparison.png
figures/Figure_K_Robustness_Comparison.tif
data/paper/k_robustness_comparison_paper_source_data.csv
data/paper/k_robustness_comparison_paper_pra_fidelity_curves.csv
```

This comparison is a robustness design benchmark, not yet a full microscopic
claim. The passive and driven calculations use different effective
descriptions, so the next validation step is a laboratory-frame robustness scan
including drive errors, finite sideband separation, and parameter noise.

## Production Laboratory-Frame Floquet Validation

The central-sideband effective passage was validated against the microscopic
laboratory-frame finite-mode dynamics in
`scripts/validate_floquet_effective_model.m`. The validation isolates the
Floquet passage after the local `|eg>` preparation and includes the CD term in
both descriptions. The laboratory-frame state is transformed to the slow frame
with the local micromotion phases before comparison.

The production run uses:

```text
Nc = 2004
nu/xi = 8
u0 = 0.9
xi*T = 20
g/xi = 0.1
```

The final effective and laboratory-frame slow-frame Bell fidelities are
`0.995881` and `0.995179`, respectively. The maximum pointwise Bell-fidelity
difference is `1.17e-3`, and the maximum slow-frame state infidelity is
`1.70e-3`. The sideband diagnostic gives
`Delta_F/max|g_m| = 2308` and `Delta_F/max|betaDot| = 56.4`.

The paper-ready validation figure and source data are written to:

```text
figures/Figure_Floquet_Effective_Validation.pdf
figures/Figure_Floquet_Effective_Validation.png
figures/Figure_Floquet_Effective_Validation.tif
data/paper/floquet_effective_validation_paper_source_data.csv
data/paper/floquet_effective_validation_paper_summary.json
```

## Production Control-Parameter Robustness Maps

The final Bell fidelity of the full two-stage Floquet-CD protocol was scanned
against four pairs of control errors in
`scripts/scan_protocol_parameter_robustness.m`. The calculation uses the
validated central-sideband effective model and reports the final unconditional
fidelity `F_{Psi+}(t_f)`.

The production run uses:

```text
Nc = 804
nu/xi = 8
u0 = 0.9
xi*T_pi = 2
xi*T_0 = 20
g/xi = 0.1
```

The final fidelity remains above `0.995` over the full `+/-10%` calibration
window for the Floquet-renormalized couplings `u1` and `u2`. For loading-pulse
errors, the minimum fidelity over `|delta A_pi/A_pi| <= 0.10` and
`|delta phi_pi| <= 0.1*pi` is `0.971510`. For the CD-assisted passage, the
scan remains above `0.995` for `|delta alpha_CD| <= 0.025` over the full
duration window `-0.5 <= delta T/T_0 <= 0.5`. The most restrictive parameter
is the relative atomic detuning: over all tested common detunings
`|Delta_c|/xi <= 0.05`, the final fidelity remains above `0.99` for
`|Delta_r|/xi <= 0.0125`.

The paper-ready figure and source data are written to:

```text
figures/Figure_Protocol_Parameter_Robustness.pdf
figures/Figure_Protocol_Parameter_Robustness.png
figures/Figure_Protocol_Parameter_Robustness.tif
data/paper/protocol_parameter_robustness_paper_source_data.csv
data/paper/protocol_parameter_robustness_paper_summary.json
```
