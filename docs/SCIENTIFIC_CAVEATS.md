# Scientific Caveats

This file tracks statements that must remain carefully qualified while the new
Floquet Bell-BIC manuscript is being developed.

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

is a target atomic direction for the chosen symmetric geometry. The code must
still verify the full instantaneous BIC or BIC subspace, including any photonic
dressing and finite-grid degeneracies.

## Adiabatic Following

Slow passage is not guaranteed to be monotonic or conventionally adiabatic in
an embedded, gapless setting. Nonmonotonic leakage or finite-size artifacts are
scientifically meaningful outcomes, not failures to hide.

## Counterdiabatic Control

The minimal atomic counterdiabatic term can be used as a benchmark only. It
should not be described as the exact counterdiabatic gauge potential of the
full BIC unless the photonic component is included.

The first reduced-grid result shows that this atomic CD term almost perfectly
tracks the intended atomic path in the effective model. In manuscript language,
this should be framed as evidence for the control principle and as a benchmark
for nonadiabatic error, not yet as a complete microscopic protocol.
