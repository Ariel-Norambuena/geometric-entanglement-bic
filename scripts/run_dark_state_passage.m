function summary = run_dark_state_passage(runMode)
%RUN_DARK_STATE_PASSAGE Development entry point for the effective passage.
%
% This script runs the exact finite-mode central-sideband effective model.
% It is meant as the first reproducible development calculation, not as a
% production claim. Production runs should increase Nc and perform solver,
% sideband, and finite-size convergence scans.

if nargin < 1
    runMode = "dev";
else
    runMode = string(runMode);
end

if runMode ~= "dev"
    error('Only runMode="dev" is implemented at this stage.');
end

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
addpath(fullfile(projectRoot, 'src', 'matlab'));

dataDir = fullfile(projectRoot, 'data', 'development');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

params.Nc = 4 * 51;
params.xi = 1.0;
params.wc = 0.0;
params.g = 0.1;
params.n1 = 6;
params.n2 = 6;
params.x1 = 0;
params.Dx = 2;
params.Omega0 = params.wc;

u0 = 0.9;
T = 20 / params.xi;
nu = 8 * params.xi;

model = giant_atom_static_model(params);
controls = floquet_controls(u0, T, nu);

tEval = linspace(0, T, 1001);
psi0 = zeros(model.dimension, 1);
psi0(1) = 1;

options = odeset( ...
    'RelTol', 1e-9, ...
    'AbsTol', 1e-11, ...
    'MaxStep', min(T/1000, 0.05));

timer = tic;
[tOut, psiOutRows] = ode113(@(t, psi) floquet_effective_rhs(t, psi, model, controls), ...
    tEval, psi0, options);
runtimeSeconds = toc(timer);

psi = transpose(psiOutRows);
obs = observables_single_excitation(psi);
target = instantaneous_atomic_direction(tOut.', controls);
targetOverlap = abs(conj(target(1, :)) .* psi(1, :) + conj(target(2, :)) .* psi(2, :)).^2;
targetConditionalFidelity = targetOverlap ./ obs.atomicPopulation;

sidebandReport = sideband_validity(model, controls, linspace(0, T, 101), 6);

summary = struct();
summary.runMode = runMode;
summary.model = "central-sideband effective finite-mode dynamics";
summary.Nc = model.Nc;
summary.dimension = model.dimension;
summary.T = T;
summary.nu = nu;
summary.u0 = u0;
summary.rtol = 1e-9;
summary.atol = 1e-11;
summary.maxStep = min(T/1000, 0.05);
summary.runtimeSeconds = runtimeSeconds;
summary.finalNormError = abs(obs.norm(end) - 1);
summary.maxNormError = max(abs(obs.norm - 1));
summary.finalAtomicPopulation = obs.atomicPopulation(end);
summary.finalPhotonicPopulation = obs.photonicPopulation(end);
summary.finalConcurrence = obs.concurrence(end);
summary.finalConditionalConcurrence = obs.conditionalConcurrence(end);
summary.finalBellPlusConditionalFidelity = obs.conditionalBellPlusFidelity(end);
summary.finalTargetConditionalFidelity = targetConditionalFidelity(end);
summary.sideband = sidebandReport;

matPath = fullfile(dataDir, 'dark_state_passage_effective_dev.mat');
jsonPath = fullfile(dataDir, 'dark_state_passage_effective_dev_summary.json');
save(matPath, 'tOut', 'psi', 'obs', 'target', 'targetConditionalFidelity', ...
    'summary', 'params', 'u0', 'T', 'nu');
write_json(jsonPath, summary);

fprintf('Development passage data written to:\n  %s\n  %s\n', matPath, jsonPath);
fprintf('Final conditional Bell fidelity: %.12f\n', summary.finalBellPlusConditionalFidelity);
fprintf('Final conditional target fidelity: %.12f\n', summary.finalTargetConditionalFidelity);
fprintf('Maximum norm error: %.3e\n', summary.maxNormError);
end

function target = instantaneous_atomic_direction(t, controls)
[u1, u2] = controls.u(t);
normalization = sqrt(u1.^2 + u2.^2);
target = [
    u2 ./ normalization;
    u1 ./ normalization
];
end

function write_json(path, value)
text = jsonencode(value, 'PrettyPrint', true);
fid = fopen(path, 'w');
if fid < 0
    error('Could not open %s for writing.', path);
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', text);
delete(cleaner);
end
