# Previous-Paper Baseline for Floquet Bell-BIC Work

This document freezes the numerical baseline inherited from the previous
giant-atom BIC manuscript before adding Floquet control. The purpose is to
protect the old Fig. 3 dynamics from accidental changes while the new
time-dependent code is developed.

## Scope

The baseline concerns the exact finite-mode single-excitation dynamics used
for the left part of Fig. 3 in the previous paper. The right-panel decay-rate
maps are analytic on-shell Markovian diagnostics and are kept separate.

## Numerical Convention

The published MATLAB figure script uses positive quasi-momenta
`k in [0, pi)` and represents the two propagation branches explicitly as
right- and left-moving modes. With `NcPositive = 4*501 = 2004`, the Hamiltonian
dimension is

```text
2 + 2*NcPositive = 4010.
```

This is equivalent to a full Brillouin-zone representation with the resonant
pair `+/- k` retained symmetrically, but the basis ordering differs from the
new Floquet brief's proposed full-zone state vector.

## Parameters

```text
xi = 1
wc = 0
g = 0.1 xi
n1 = n2 = 6
x1 = 0
kstar0 = pi/2
Omega0 = wc - 2 xi cos(kstar0) = 0
lambda detuning = 10%
theta detuning = 0.1 * 2*pi
kstar detuning = 10%
```

For the symmetric Bell sector, `Dx = 2`, so `kstar0 * Dx = pi`. For the
antisymmetric Bell sector, `Dx = 0`.

## Frozen Observables

The machine-readable baseline is stored in:

```text
data/baseline_previous_paper.json
```

It records concurrence and norm at:

```text
xi*t = 0, 100, 200, 400, 600, 800, 1000.
```

The current reference values include:

```text
plus_ideal:            C(1000) = 0.9796771010812017
plus_lambda_detuned:   C(1000) = 0.09801716326134342
minus_ideal:           C(1000) = 1.0000000000000000
minus_lambda_detuned:  C(1000) = 0.1428722342897660
```

The maximum norm error in the stored production-grid cases is below
`6e-15`.

## Regression Test

Run:

```matlab
cd '<repository root>'
run('tests/test_previous_paper_regression.m')
```

The test recomputes the selected observables and fails if any concurrence
value changes by more than `5e-8`, or if the norm error exceeds `1e-8`.

## Important Caveat

The baseline should not be interpreted as a Floquet result. It only validates
the static previous-paper dynamics. The new Floquet implementation must keep a
strict distinction between the laboratory-frame driven model, the rotating
frame with micromotion, the central-sideband effective model, and any
Markovian comparison.
