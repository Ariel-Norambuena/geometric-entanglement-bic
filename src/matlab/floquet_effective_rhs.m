function dpsi = floquet_effective_rhs(t, psi, model, controls)
%FLOQUET_EFFECTIVE_RHS Matrix-free central-sideband effective dynamics.
%
% State ordering:
%   psi = [c1; c2; phi_0; ...; phi_(Nc-1)].
%
% The effective couplings are f_ik(t) = u_i(t) g_ik.

c1 = psi(1);
c2 = psi(2);
phi = psi(3:end);

[u1, u2] = controls.u(t);
f1k = u1 * model.g1k;
f2k = u2 * model.g2k;

dc1 = -1i * (model.Omega0 * c1 + transpose(f1k) * phi);
dc2 = -1i * (model.Omega0 * c2 + transpose(f2k) * phi);
dphi = -1i * (model.wk .* phi + conj(f1k) * c1 + conj(f2k) * c2);

dpsi = [dc1; dc2; dphi];
end
