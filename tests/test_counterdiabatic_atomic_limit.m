%% test_counterdiabatic_atomic_limit.m
% Verify the sign and normalization of the atomic counterdiabatic term.

clear; clc;

testDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(testDir);
addpath(fullfile(projectRoot, 'src', 'matlab'));

params.Nc = 4 * 5;
params.xi = 1.0;
params.wc = 0.0;
params.g = 0.0;
params.n1 = 6;
params.n2 = 6;
params.x1 = 0;
params.Dx = 2;
params.Omega0 = 0;

model = giant_atom_static_model(params);
T = 10;
controls = floquet_controls(0.9, T, 8);
tEval = linspace(0, T, 501);

psi0 = zeros(model.dimension, 1);
psi0(1) = 1;

solverOptions = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, 'MaxStep', T/1000);

cdOptions.counterdiabaticScale = 1;
[tOut, psiOutRows] = ode113(@(t, psi) floquet_effective_rhs(t, psi, model, controls, cdOptions), ...
    tEval, psi0, solverOptions);

psi = transpose(psiOutRows);
theta = controls.theta(tOut.');
target = [cos(theta); sin(theta)];
targetError = vecnorm(psi(1:2, :) - target, 2, 1);
obs = observables_single_excitation(psi);

assert(max(targetError) < 5e-8);
assert(max(abs(obs.norm - 1)) < 5e-10);
assert(abs(obs.conditionalBellPlusFidelity(end) - 1) < 5e-8);

adiabaticOptions.counterdiabaticScale = 0;
[~, psiAdiabaticRows] = ode113(@(t, psi) floquet_effective_rhs(t, psi, model, controls, adiabaticOptions), ...
    tEval, psi0, solverOptions);
psiAdiabatic = transpose(psiAdiabaticRows);

assert(abs(psiAdiabatic(2, end)) < 1e-10);

fprintf('Counterdiabatic atomic-limit test passed.\n');
