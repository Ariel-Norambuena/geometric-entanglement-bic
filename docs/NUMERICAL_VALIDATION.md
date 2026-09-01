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

## Pending Validation Before Scientific Claims

- Production-grid dark-state-passage duration scans.
- Full laboratory-frame versus effective-model comparison after removing
  micromotion.
- Frozen instantaneous BIC branch construction, including photonic dressing.
- Real-space separation of localized photonic dressing from outgoing radiation.
- Solver tolerance convergence.
- Finite-size convergence using grids that retain `k = +/- pi/2`.
- Robustness scans over amplitude, phase, geometry, detuning, and timing.

No manuscript-level claim about successful Bell-BIC preparation should be made
until the pending items above are completed.

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
passage is a nonadiabatic-following problem in the atomic direction. It is not
yet proof of a laboratory-frame Bell-BIC protocol; the lab-frame comparison
and instantaneous dressed BIC validation remain pending.
