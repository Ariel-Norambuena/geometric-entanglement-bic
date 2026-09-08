# Scientific Caveats

This file tracks statements that must remain carefully qualified while the new
Floquet Bell-BIC manuscript is being developed.

The September 2026 revision is documented in `PRR_REVISION_20260908.md`.
It supersedes the initial reduced-grid interpretation with an analytic compact
BIC, finite-grid checks, a CD-free hold, retrieval, and matched benchmarks.
The geometric amplitude construction is additionally recovered and tested in
`GEOMETRIC_DERIVATION_LIMITS.md`, with explicit attribution to A. R. Legon.
Its conditional atomic identities are exact in the stated family; vacuum-field
and frozen-memory approximations have separate, numerically quantified limits.

## Bell-BIC Language

Use "Bell-BIC" only for a bound state whose normalized atomic component is a
Bell state. A high Bell fidelity of the atomic component alone does not prove
that the full state is a BIC.

## Photonic Population

Photonic population is not automatically leakage. A true BIC can carry a
localized, normalizable photonic dressing. Leakage should be associated with
radiation leaving the bound region or reaching numerical boundaries.

## Effective Model

The central-sideband Hamiltonian is an approximation to the driven microscopic
model. Its validity must be supported by sideband-separation diagnostics and
by direct comparison with the laboratory-frame model after micromotion removal.

## Instantaneous Branch

The vector

```text
d(t) = [u2(t), u1(t)] / sqrt(u1(t)^2 + u2(t)^2)
```

is the atomic direction of the selected compact dressed BIC for connections
(0,n),(2,n+2). The September revision verifies its full eigenstate, derivative,
and the n=6 degeneracy in `test_compact_bic_branch.m`. Other geometries still
require an independent construction and validation.

## Adiabatic Following

Slow passage is not guaranteed to be monotonic or conventionally adiabatic in
an embedded, gapless setting. Nonmonotonic leakage or finite-size artifacts are
scientifically meaningful outcomes, not failures to hide.

## Counterdiabatic Control

The minimal atomic counterdiabatic term can be used as a benchmark only. It
should not be described as the exact counterdiabatic gauge potential of the
full BIC unless the photonic component is included.

The September revision derives the residual left by the atomic CD and compares
it with a state-selective dressed control that also requires photonic access.
Atomic exchange alone is already an entangling resource; the guide-path benefit
is assessed through dressing, retention, and retrieval. The full laboratory
loading, hold, and retrieval sequence has not been validated.

## K-Mismatch Robustness

The earlier PRA-style Markov/passage comparison was unbalanced and is now
archival only. Its generator is superseded by a separately tested static
reference. The revised figure instead compares finite-mode calculations with
identical geometry, input state, exchange pulse, and hold time. It resolves
preparation from storage sensitivity; sudden attachment is not bandwidth
matched. Neither this comparison nor the post-loading beta-error map establishes
full laboratory-frame robustness, dissipation tolerance, or optimality.
