function results = scan_passage_duration(runMode)
%SCAN_PASSAGE_DURATION Development scan for effective dark-state passage.
%
% The scan is intentionally labeled as development data. It uses a reduced
% full-Brillouin-zone grid to probe qualitative behavior before production
% convergence.

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

durationsXi = [2 5 10 20 50 100];
u0 = 0.9;
nu = 8 * params.xi;
numOutputTimes = 1001;

model = giant_atom_static_model(params);
results = repmat(empty_result(), numel(durationsXi), 1);

for idx = 1:numel(durationsXi)
    T = durationsXi(idx) / params.xi;
    controls = floquet_controls(u0, T, nu);
    tEval = linspace(0, T, numOutputTimes);

    psi0 = zeros(model.dimension, 1);
    psi0(1) = 1;

    maxStep = min(T/1000, 0.05);
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-11, 'MaxStep', maxStep);

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

    results(idx).durationXi = durationsXi(idx);
    results(idx).Nc = model.Nc;
    results(idx).dimension = model.dimension;
    results(idx).nuOverXi = nu / params.xi;
    results(idx).u0 = u0;
    results(idx).rtol = 1e-9;
    results(idx).atol = 1e-11;
    results(idx).maxStep = maxStep;
    results(idx).runtimeSeconds = runtimeSeconds;
    results(idx).maxNormError = max(abs(obs.norm - 1));
    results(idx).finalAtomicPopulation = obs.atomicPopulation(end);
    results(idx).finalPhotonicPopulation = obs.photonicPopulation(end);
    results(idx).finalConcurrence = obs.concurrence(end);
    results(idx).finalConditionalConcurrence = obs.conditionalConcurrence(end);
    results(idx).finalBellPlusConditionalFidelity = obs.conditionalBellPlusFidelity(end);
    results(idx).finalTargetConditionalFidelity = targetConditionalFidelity(end);
    results(idx).worstTargetConditionalInfidelity = max(1 - targetConditionalFidelity);
    results(idx).sideband = sidebandReport;

    fprintf(['T*xi=%6.1f | F_Bell_cond=%.6f | C_cond=%.6f | ' ...
        'P_at=%.6f | max norm err=%.2e\n'], ...
        results(idx).durationXi, ...
        results(idx).finalBellPlusConditionalFidelity, ...
        results(idx).finalConditionalConcurrence, ...
        results(idx).finalAtomicPopulation, ...
        results(idx).maxNormError);
end

matPath = fullfile(dataDir, 'duration_scan_effective_dev.mat');
jsonPath = fullfile(dataDir, 'duration_scan_effective_dev_summary.json');
save(matPath, 'results', 'params', 'durationsXi', 'u0', 'nu', 'numOutputTimes');
gridSummary = struct( ...
    'Nc', model.Nc, ...
    'dimension', model.dimension, ...
    'kConvention', model.kConvention);

write_json(jsonPath, struct( ...
    'runMode', runMode, ...
    'model', "central-sideband effective finite-mode dynamics", ...
    'grid', gridSummary, ...
    'results', results));

fprintf('Development duration scan written to:\n  %s\n  %s\n', matPath, jsonPath);
end

function result = empty_result()
result.durationXi = [];
result.Nc = [];
result.dimension = [];
result.nuOverXi = [];
result.u0 = [];
result.rtol = [];
result.atol = [];
result.maxStep = [];
result.runtimeSeconds = [];
result.maxNormError = [];
result.finalAtomicPopulation = [];
result.finalPhotonicPopulation = [];
result.finalConcurrence = [];
result.finalConditionalConcurrence = [];
result.finalBellPlusConditionalFidelity = [];
result.finalTargetConditionalFidelity = [];
result.worstTargetConditionalInfidelity = [];
result.sideband = [];
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
