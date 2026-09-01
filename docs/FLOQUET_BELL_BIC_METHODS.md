# Floquet Bell-BIC Numerical Methods

This document describes the numerical implementation planned for the new
Floquet-controlled Bell-BIC manuscript.

## Model Separation

The code keeps four descriptions conceptually separate:

1. Microscopic laboratory-frame driven model.
2. Rotating-frame model with the exact local micromotion phase.
3. Central-sideband effective Hamiltonian.
4. Markovian or on-shell rates, used only as comparison diagnostics.

The principal dynamics are exact within the finite-mode single-excitation
representation. No Lindblad or Markovian dynamics are used for the controlled
passage.

## Full-Brillouin-Zone Convention

The new Floquet code uses

```text
k_m = -pi + 2*pi*m/Nc,  m = 0,...,Nc-1.
```

The state vector is ordered as

```text
psi = [c1, c2, phi_0, ..., phi_(Nc-1)].
```

This differs from the previous Fig. 3 MATLAB figure script, which uses
positive `k in [0, pi)` plus explicit right/left modes. The two conventions
are physically equivalent when normalized consistently, but they are kept
separate in code and documentation.

## Effective Central-Sideband Dynamics

The central-sideband effective model uses

```text
f_1k(t) = u_1(t) g_1k
f_2k(t) = u_2(t) g_2k
u_i(t) = J_0[beta_i(t)].
```

The RHS is implemented matrix-free:

```text
dc1  = -i * (Omega0*c1 + sum_k f_1k(t) phi_k)
dc2  = -i * (Omega0*c2 + sum_k f_2k(t) phi_k)
dphi = -i * (wk*phi_k + conj(f_1k)c1 + conj(f_2k)c2)
```

Current implementation:

```text
src/matlab/floquet_effective_rhs.m
```

## Laboratory-Frame Driven Dynamics

The laboratory model keeps the bare couplings and time-dependent atomic
frequencies:

```text
Omega_i(t) = Omega0 + nu beta_i(t) cos(nu t) + betaDot_i(t) sin(nu t).
```

Current implementation:

```text
src/matlab/floquet_lab_rhs.m
```

Laboratory and effective states must be compared only after removing the local
atomic micromotion:

```text
c_i^rot(t) = exp(+i beta_i(t) sin(nu t)) c_i^lab(t).
```

## Dark-State-Passage Controls

The default effective coupling path is

```text
u1(t) = u0 sin(theta(t))
u2(t) = u0 cos(theta(t))
theta(t) = (pi/4) * (10s^3 - 15s^4 + 6s^5),  s=t/T.
```

The inverse Bessel map is solved with a bracketed root solver on the first
monotonic branch:

```text
0 <= beta <= 2.404825557695773.
```

Current implementation:

```text
src/matlab/floquet_controls.m
src/matlab/besselj0_inverse_first_branch.m
```

## Observables

The current observables include norm, atomic population, photonic population,
unconditional concurrence, conditional concurrence, unconditional Bell-plus
fidelity, and conditional Bell-plus fidelity.

Current implementation:

```text
src/matlab/observables_single_excitation.m
```

Photonic population is not automatically interpreted as leakage. Future
real-space diagnostics must separate localized bound dressing from emitted
radiation.
