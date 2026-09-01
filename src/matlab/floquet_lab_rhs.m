function dpsi = floquet_lab_rhs(t, psi, model, controls)
%FLOQUET_LAB_RHS Matrix-free microscopic laboratory-frame dynamics.
%
% This model keeps the physical time-dependent atomic frequencies Omega_i(t)
% and the bare atom-waveguide couplings g_ik. It must be compared with the
% effective model only after applying the local micromotion transformation.

c1 = psi(1);
c2 = psi(2);
phi = psi(3:end);

omegaAtoms = controls.atomicFrequencies(t, model.Omega0);

dc1 = -1i * (omegaAtoms(1) * c1 + transpose(model.g1k) * phi);
dc2 = -1i * (omegaAtoms(2) * c2 + transpose(model.g2k) * phi);
dphi = -1i * (model.wk .* phi + conj(model.g1k) * c1 + conj(model.g2k) * c2);

dpsi = [dc1; dc2; dphi];
end
