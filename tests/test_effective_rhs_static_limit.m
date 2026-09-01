%% test_effective_rhs_static_limit.m
% Basic structural test for the matrix-free effective-model RHS.

clear; clc;

testDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(testDir);
addpath(fullfile(projectRoot, 'src', 'matlab'));

params.Nc = 4 * 25;
params.xi = 1.0;
params.wc = 0.0;
params.g = 0.1;
params.n1 = 6;
params.n2 = 6;
params.x1 = 0;
params.Dx = 2;
params.Omega0 = 0;

model = giant_atom_static_model(params);
controls = floquet_controls(0.9, 50, 8);

rng(7);
psi = randn(model.dimension, 1) + 1i * randn(model.dimension, 1);
psi = psi / norm(psi);

dpsi = floquet_effective_rhs(0, psi, model, controls);
assert(isequal(size(dpsi), size(psi)));

% For a Hermitian finite-mode Hamiltonian, d/dt <psi|psi> = 0.
normDerivative = 2 * real(psi' * dpsi);
assert(abs(normDerivative) < 1e-12);

obs = observables_single_excitation(psi);
assert(abs(obs.norm - 1) < 1e-12);
assert(obs.atomicPopulation + obs.photonicPopulation > 1 - 1e-12);

report = sideband_validity(model, controls, linspace(0, 50, 11), 4);
assert(report.deltaF > 0);
assert(report.deltaFOverSidebandCoupling > 1);
assert(report.deltaFOverBetaDot > 1);

fprintf('Effective RHS static-limit tests passed.\n');
