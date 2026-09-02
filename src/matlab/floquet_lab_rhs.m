function dpsi = floquet_lab_rhs(t, psi, model, controls, options)
%FLOQUET_LAB_RHS Matrix-free microscopic laboratory-frame dynamics.
%
% This model keeps the physical time-dependent atomic frequencies Omega_i(t)
% and the bare atom-waveguide couplings g_ik. It must be compared with the
% effective model only after applying the local micromotion transformation.
%
% Optional counterdiabatic correction:
%   H_CD^rot = alpha * thetaDot(t) * sigma_y
% is implemented in the laboratory frame with the local micromotion phases.
% The default alpha=0 leaves the bare Floquet laboratory model unchanged.

if nargin < 5
    options = struct();
end

options = with_default(options, 'counterdiabaticScale', 0);

c1 = psi(1);
c2 = psi(2);
phi = psi(3:end);

omegaAtoms = controls.atomicFrequencies(t, model.Omega0);
omegaCD = options.counterdiabaticScale * controls.thetaDot(t);

[beta1, beta2] = controls.beta(t);
chi1 = beta1 .* sin(controls.nu * t);
chi2 = beta2 .* sin(controls.nu * t);
cdPhase12 = exp(1i * (chi2 - chi1));
cdPhase21 = exp(1i * (chi1 - chi2));

dc1 = -1i * (omegaAtoms(1) * c1 + transpose(model.g1k) * phi);
dc2 = -1i * (omegaAtoms(2) * c2 + transpose(model.g2k) * phi);
dphi = -1i * (model.wk .* phi + conj(model.g1k) * c1 + conj(model.g2k) * c2);

dc1 = dc1 - omegaCD * cdPhase12 * c2;
dc2 = dc2 + omegaCD * cdPhase21 * c1;

dpsi = [dc1; dc2; dphi];
end

function options = with_default(options, name, value)
if ~isfield(options, name)
    options.(name) = value;
end
end
