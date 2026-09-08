# PRR revision: reproducibility and interpretation

This analysis uses connections (0,6),(2,8), g/xi=0.1, u0=0.9, and
the full Brillouin-zone convention implemented in `giant_atom_static_model`.
The original static-BIC scripts and archived figures are preserved.

The boxed-figure update additionally recovers A. R. Legon's geometric
amplitude construction and tests its exact and approximate limits. See
`GEOMETRIC_DERIVATION_LIMITS.md` and `data/prr_revision/geometric_appendix/`.

## New checks

- `compact_bic_branch.m` gives the normalized compact eigenstate, its
  parameter derivative, atomic weight, and the localized photon direction.
- `run_prr_revision_validation.m` propagates atomic CD, ideal dressed CD,
  fixed-guide exchange, Bell-then-attachment, isolated exchange, and no CD.
  It checks a 100/xi hold with CD off, reversed-path retrieval, finite grids
  Nc=804,2004,4004, matched K errors, durations, and carrier convergence.
- `scan_prr_beta_errors.m` scans a 21-by-21 grid of multiplicative beta
  errors after ideal loading, followed by the same hold. Errors in beta
  shift Bessel zeros; multiplying u does not.
- `giant_atom_markov_generator.m` supplies a separate corrected static
  reference with the tight-binding 1/sin(K) factor and a trace-preserving
  Lindblad generator. It does not generate the revised comparison figure.

The selected compact BIC has an atomic-CD residual |etaDot|, where
eta=atan[(g*u0/(sqrt(2)*xi))*sin(2*theta)]. The ideal reference adds
i*etaDot*(|p><D|-|D><p|). This requires photonic access in addition to
the two-body atomic exchange. It is not asserted to follow from local
longitudinal modulation alone.

## Results and scope

Nominal preparation (xi*T=20, after ideal loading):

| Quantity | Value |
|---|---:|
| Unconditional Bell fidelity | 0.99588074205 |
| Unconditional concurrence | 0.99587079906 |
| Fidelity to the selected dressed BIC | 0.99992176357 |
| Atomic weight of the exact BIC | 0.99596633634 |
| Bell fidelity after xi*t_hold=100 | 0.99588754379 |
| Retrieved eg population after reversed passage | 0.99986382230 |

The nominal numbers agree across Nc=804,2004,4004 within 1e-10.
The full-state dressed reference tracks the target within 1e-12.
All new wavefunction propagations conserve norm to their documented
tolerances. Grid spacing is never used as a physical adiabatic gap.

Across the matched |deltaK|/K0 <= 0.1 scan, the proposed preparation
fidelity is >=0.99210, but the final hold fidelity drops to 0.97620.
The hold result is 0.94289 for fixed-guide exchange and 0.97328 for
Bell-then-attachment at the worst sampled K values. The reference
attachment is a sudden switch, not a bandwidth-matched control.

For beta errors of +/-10% applied after ideal loading, the minimum
preparation fidelity is 0.99439 and the minimum final hold fidelity
is 0.99458. This scan does not validate loading in the presence of
residual sideband emission or beta errors during the local pulse.

The carrier sweep (nu/xi=4,6,8,12,16) compares 401 sampled times in the
passage only. The archived fine trace at nu/xi=8 is retained separately.
Maxima on the coarser sweep grid should not be confused with continuous
time maxima or the finer archived trace.

## Replotting and tests

Python dependencies are NumPy, SciPy, and Matplotlib. Figures are vector
PDFs with rasterized heatmap cells and vector axes/text/curves. Grid cells
are displayed directly, without smoothing. Titles, legends, and axes use
a common font scale; legends occupy reserved space outside the data axes.

```matlab
run('tests/test_compact_bic_branch.m')
run('tests/test_markov_generator.m')
run('tests/test_counterdiabatic_atomic_limit.m')
run('tests/test_floquet_controls.m')
run('tests/test_effective_rhs_static_limit.m')
run('tests/test_lab_frame_cd_micromotion.m')
```

The artifact verifier checks finite-grid agreement, dressed-state
normalization and tracking, loading-phase invariance, pulse-area scaling,
and bounds of all generated fidelities. Legacy static-paper regression is
separate because those scripts and Hamiltonians were not changed.

## Interpretation corrections

The frozen retarded resolvent is a stationary family at independent t0;
it is not an ordinary Laplace solution of the nonstationary memory equation.
Both resonant emission directions and the real-energy eigenvalue equation
must hold. A regulator is not a physical momentum window, and cancellation
of coupling derivatives is an additional design constraint. The full
self-energy is retained without assuming it is normal.

The atomic exchange creates a Bell state even for g=0. The guide-control
benefit is therefore assessed by loading, field dressing, retention, and
retrieval. Photon population includes localized dressing; spatially outgoing
population is reported separately. Loading-phase invariance is a symmetry
check, not relative Bell-phase robustness.

## Remaining physical limitations

The new hold and retrieval calculations are effective and lossless.
Atomic relaxation, dephasing, resonator loss, laboratory-frame validation
of the entire sequence, multi-excitation effects during imperfect loading,
and device-specific exchange bandwidth remain outside the demonstrated
scope. The duration scan is not a global optimum under fixed resources.

The revised manuscript replaces the earlier unbalanced passive comparison.
Archived `scripts/scan_k_robustness_comparison.m` and its data are retained
for provenance, but their Markovian implementation and interpretation are
superseded by the corrected reference and matched finite-mode comparison.
