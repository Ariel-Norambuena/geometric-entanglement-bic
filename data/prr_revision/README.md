# September 2026 source data

All quantities use xi=1. Unless stated otherwise, the input is |eg,0>,
g=0.1, u0=0.9, connections (0,6),(2,8), and a passage T=20 followed by
a hold of 100 with the atomic exchange off. These are lossless
central-sideband calculations; carrier_convergence.csv is the exception.

## File guide

- `revision_summary.json`: production parameters, tolerances, and key checks.
- `nominal_trajectories.mat`: the `nominal` structure contains sampled time,
  F (unconditional Bell fidelity), C, FB (selected dressed-state fidelity),
  atomic/photon populations, and Poutside (photons outside sites 0 through 8).
  Preparation is MATLAB index 401; the rest is the hold, with time measured
  from the start of the passage. `prepF` and `prepFB` are endpoint values.
- `nominal_<method>_states.mat`: the `trajectory` structure also contains
  `psi`, with rows (c1,c2,phi_k) and columns corresponding to sampled times.
  Splitting protocols keeps every file below common Git hosting size limits.
- `retrieval.mat`: the `retrieval` structure contains the reverse passage
  after the nominal atomic-CD hold, including the full state and eg population.
- `matched_k_benchmark.csv`: methods 1=atomic CD, 2=fixed-guide exchange,
  3=isolated exchange then attachment. deltaK_fraction is deltaK/K0, not percent.
  All methods share the same detuned final Hamiltonian during the hold.
- `matched_final_states.mat`: prepared and end-of-hold state columns follow
  the row order in the corresponding CSV.
- `duration_benchmark.csv`: methods 1=atomic CD, 2=ideal dressed control,
  3=fixed-guide exchange, 4=attachment, 5=no CD. Durations share exchange
  area, not fixed peak amplitude. Full endpoint states are in the matching MAT.
- `grid_convergence.csv`: the `minimum_F_Bell` field covers the entire
  preparation-plus-hold trace, including the initial Bell overlap of 1/2.
  It must not be reported as the minimum during the hold.
- `carrier_convergence.csv`: microscopic passage versus the effective model,
  with maxima taken on 401 output times. The separate archived nu=8 trace in
  `data/paper/` uses a denser grid and therefore slightly different maxima.
- `beta_error_maps.csv`: independent fractional errors in beta1,beta2 after
  ideal loading. The 441 rows are ordered by beta2, then beta1. Final states
  are columns of `states` in beta_error_maps.mat. Errors during loading,
  local-pulse leakage, and higher-excitation sectors are not covered.
- `beta_map_grid_checks.csv`: the nominal beta point and four corners repeated
  with Nc=804,2004,4004; no claim of checking every map point at all three sizes.
- `artifact_checks.json`: output of the saved-data consistency verifier.

The nominal methods `atomic`, `dressed`, `fixed`, `attach`, `isolated`, `none`
are distinguished in the generation script. In particular, `dressed` requires
photonic control, and `isolated` remains disconnected even during the hold.
Its overlap with the nominal guide BIC is a target overlap, not an eigenstate
claim for the disconnected Hamiltonian.

## Reproduce

From the repository root in MATLAB R2025b:

```matlab
addpath('scripts')
run_prr_revision_validation('paper')
scan_prr_beta_errors('paper')
```

Then run `python scripts/verify_prr_revision_artifacts.py` and
`python scripts/make_prr_revision_figures.py`. The figure script also reads
three archived CSVs in `data/paper/`; they are kept with their original
parameters and precision. No fitted or smoothed fidelity values are substituted
for simulated data. See `docs/PRR_REVISION_20260908.md` for interpretation.
