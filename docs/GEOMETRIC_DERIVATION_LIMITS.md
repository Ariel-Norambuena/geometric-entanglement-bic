# Original geometric derivation: recovered conditions and numerical limits

The original amplitude-selection construction was developed by Alexis R. Legon.
The revised appendix retains this contribution and identifies the exact
conditional-atomic results separately from two approximations.

## Recovered statements

At K*n_j=(2*l_j+1)*pi, the on-shell couplings vanish. Canceling their first
momentum derivative additionally selects

```text
a2/a1 = -lambda_F exp(-i*K*Dx), lambda_F = n1*u1/(n2*u2).
```

The two roots require exp(2*i*K*Dx)=1 for nonzero components. For K=pi/2,
n_j=4*l_j+2 and even Dx, the selected direction belongs to an exact dressed
BIC family. The conditional concurrence is exactly 2*abs(lambda_F)/(1+lambda_F^2).
Thus local Bessel factors can compensate unequal connection lengths in this
specified family. The full real-energy equation is also checked.

Equal lengths and equal frozen couplings permit exact unitary Bell-channel
diagonalization. The stationary Laplace construction is exact at fixed controls;
it is not a solution of the nonstationary memory equation at a running time.
Identical coincident connection sets provide a further exact case: the selected
atomic state is photon free, with Z=1, even at finite coupling.

## What is approximate

For braided connections, the vacuum-field ansatz neglects localized photons.
Its full-state infidelity is 1-Z. At g/xi=0.1, u0=0.9, lambda_F=1, and
(n1,n2,Dx)=(6,6,2), this is 0.00403366. The atomic conditional formula remains
exact, while unconditional concurrence and Bell fidelity acquire a factor Z.
The weak-dressing limit is explicitly scanned and approaches zero error as g^2.

Replacing u_j(t-tau) by u_j(t) in the memory is a different approximation.
`frozen_envelope_rhs.m` integrates that equation with auxiliary memory variables;
these variables are not independent physical photon modes. The comparison
does not normalize the approximate atomic vector or treat its enlarged norm
as a probability.

At xi*T=20, the maximum atomic-amplitude error is 1.3882e-4 for g/xi=0.01,
and 1.35996e-2 for g/xi=0.1. At g/xi=0.1, increasing xi*T to 80 or 160 does
not remove the accumulated error. The Bell-overlap expression differs by up
to 0.00349684 at the nominal coupling and duration. The original stationary
analysis is therefore retained, but the paper's precision claims continue
to use the full finite-mode propagation without envelope freezing.

## Reproduction

```matlab
addpath('scripts')
validate_geometric_bic_approximation('paper')
check_geometric_appendix_convergence()
run('tests/test_geometric_bic_approximation.m')
```

Data are in `data/prr_revision/geometric_appendix/`; the new result is
`figures/prr_revision/Figure_Geometric_BIC_Approximation.pdf`.
The amplitude-rule CSV covers both (6,6) and (6,10), with 41 lambda values each.
The radiation CSV resolves both roots, not just the positive wave vector.
The dressing CSV contains 41 coupling values for each geometry. The memory
comparison covers six couplings and three durations; its MAT file preserves
the compared atomic trajectories. The convergence CSV repeats four cases at
Nc=804,2004,4004 and tighter tolerances with a smaller maximum solver step.

The tests check real-space versus momentum-space state reconstruction,
Hamiltonian eigenstate residuals, conditional entanglement, the photon-free
coincident case, symmetric Bell-channel diagonalization, and equivalence of
the complete memory-variable state reconstruction at constant controls.
The saved-data verifier also checks near-root powers, the analytic dressing
weights, weak-coupling scaling, and finite-grid/solver convergence.

Original manuscript excerpts remain in the private manuscript package for
provenance. They are not overwritten or inserted unchanged as exact claims
outside the identified regime. The earlier static-paper scripts and data are
untouched. This appendix does not assert arbitrary-geometry compensation or
automatic adiabatic transport in a degenerate BIC sector.
