%% test_lab_frame_cd_micromotion.m
% Verify the laboratory-frame CD term after removing Floquet micromotion.

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
nu = 8;
controls = floquet_controls(0.9, T, nu);
tEval = linspace(0, T, 501);

psi0 = zeros(model.dimension, 1);
psi0(1) = 1;

rhsOptions.counterdiabaticScale = 1;
solverOptions = odeset( ...
    'RelTol', 1e-10, ...
    'AbsTol', 1e-12, ...
    'MaxStep', min(T / 2000, 2*pi / nu / 100));

[tOut, psiLabRows] = ode113(@(t, psi) floquet_lab_rhs(t, psi, model, controls, rhsOptions), ...
    tEval, psi0, solverOptions);

psiLab = transpose(psiLabRows);
psiSlow = micromotion_to_rotating_frame(tOut.', psiLab, controls);

theta = controls.theta(tOut.');
target = [cos(theta); sin(theta)];
targetError = vecnorm(psiSlow(1:2, :) - target, 2, 1);
obs = observables_single_excitation(psiSlow);

assert(max(targetError) < 2e-7);
assert(max(abs(obs.norm - 1)) < 5e-10);
assert(abs(obs.conditionalBellPlusFidelity(end) - 1) < 2e-7);

fprintf('Laboratory-frame CD micromotion test passed.\n');
