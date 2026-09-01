function summary = run_counterdiabatic_passage(runMode)
%RUN_COUNTERDIABATIC_PASSAGE Compare adiabatic and CD-assisted passage.
%
% This development script propagates the same finite-mode central-sideband
% effective model with and without the minimal atomic counterdiabatic term
% H_CD = thetaDot(t) sigma_y. It exports concurrence, protocol fidelity, and
% control plots for rapid scientific inspection.

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
outputDir = fullfile(projectRoot, 'outputs', 'floquet');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
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
numOutputTimes = 1501;

model = giant_atom_static_model(params);
controls = floquet_controls(u0, T, nu);
tEval = linspace(0, T, numOutputTimes);

caseSpecs = struct( ...
    'name', {'adiabatic', 'counterdiabatic'}, ...
    'plotLabel', {'adiabatic', 'adiabatic + CD'}, ...
    'counterdiabaticScale', {0, 1});

results = repmat(empty_result(), numel(caseSpecs), 1);
psiAll = cell(numel(caseSpecs), 1);

for idx = 1:numel(caseSpecs)
    psi0 = zeros(model.dimension, 1);
    psi0(1) = 1;

    rhsOptions.counterdiabaticScale = caseSpecs(idx).counterdiabaticScale;
    solverOptions = odeset( ...
        'RelTol', 1e-9, ...
        'AbsTol', 1e-11, ...
        'MaxStep', min(T/1000, 0.05));

    timer = tic;
    [tOut, psiOutRows] = ode113(@(t, psi) floquet_effective_rhs(t, psi, model, controls, rhsOptions), ...
        tEval, psi0, solverOptions);
    runtimeSeconds = toc(timer);

    psi = transpose(psiOutRows);
    obs = observables_single_excitation(psi);
    target = instantaneous_atomic_direction(tOut.', controls);
    targetOverlap = abs(conj(target(1, :)) .* psi(1, :) + conj(target(2, :)) .* psi(2, :)).^2;
    targetConditionalFidelity = targetOverlap ./ obs.atomicPopulation;

    results(idx).name = caseSpecs(idx).name;
    results(idx).counterdiabaticScale = caseSpecs(idx).counterdiabaticScale;
    results(idx).runtimeSeconds = runtimeSeconds;
    results(idx).maxNormError = max(abs(obs.norm - 1));
    results(idx).finalAtomicPopulation = obs.atomicPopulation(end);
    results(idx).finalPhotonicPopulation = obs.photonicPopulation(end);
    results(idx).finalConcurrence = obs.concurrence(end);
    results(idx).finalConditionalConcurrence = obs.conditionalConcurrence(end);
    results(idx).finalBellPlusConditionalFidelity = obs.conditionalBellPlusFidelity(end);
    results(idx).finalTargetConditionalFidelity = targetConditionalFidelity(end);
    results(idx).minimumTargetConditionalFidelity = min(targetConditionalFidelity);
    results(idx).obs = obs;
    results(idx).targetConditionalFidelity = targetConditionalFidelity;

    psiAll{idx} = psi;

    fprintf(['%-16s | C=%.6f | C_cond=%.6f | F_target_cond=%.6f | ' ...
        'P_at=%.6f | max norm err=%.2e\n'], ...
        caseSpecs(idx).plotLabel, ...
        results(idx).finalConcurrence, ...
        results(idx).finalConditionalConcurrence, ...
        results(idx).finalTargetConditionalFidelity, ...
        results(idx).finalAtomicPopulation, ...
        results(idx).maxNormError);
end

tXi = params.xi * tOut;
theta = controls.theta(tOut.');
targetConcurrence = sin(2 * theta);
[u1, u2] = controls.u(tOut.');
omegaCD = controls.thetaDot(tOut.');

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 22]);
tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on; box on; grid on;
plot(tXi, targetConcurrence, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 2.0);
plot(tXi, results(1).obs.concurrence, 'LineWidth', 1.9);
plot(tXi, results(2).obs.concurrence, 'LineWidth', 2.2);
xlabel('$\xi t$', 'Interpreter', 'latex');
ylabel('$\mathcal{C}(t)$', 'Interpreter', 'latex');
ylim([-0.03 1.03]);
legend({'target atomic path', 'adiabatic', 'adiabatic + CD'}, ...
    'Interpreter', 'latex', 'Location', 'southeast');
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
hide_axes_toolbar(gca);

nexttile; hold on; box on; grid on;
plot(tXi, results(1).targetConditionalFidelity, 'LineWidth', 1.9);
plot(tXi, results(2).targetConditionalFidelity, 'LineWidth', 2.2);
xlabel('$\xi t$', 'Interpreter', 'latex');
ylabel('$F_D^{(a)}(t)$', 'Interpreter', 'latex');
ylim([-0.03 1.03]);
legend({'adiabatic', 'adiabatic + CD'}, 'Interpreter', 'latex', 'Location', 'southeast');
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
hide_axes_toolbar(gca);

nexttile; hold on; box on; grid on;
yyaxis left;
plot(tXi, u1, 'LineWidth', 2.0);
plot(tXi, u2, 'LineWidth', 2.0);
ylabel('$u_i(t)$', 'Interpreter', 'latex');
ylim([-0.03 1.03]);
yyaxis right;
plot(tXi, omegaCD / params.xi, 'k--', 'LineWidth', 2.0);
ylabel('$\Omega_{\rm CD}(t)/\xi$', 'Interpreter', 'latex');
xlabel('$\xi t$', 'Interpreter', 'latex');
legend({'$u_1(t)$', '$u_2(t)$', '$\Omega_{\rm CD}(t)$'}, ...
    'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex');
hide_axes_toolbar(gca);

pngPath = fullfile(outputDir, 'counterdiabatic_passage_dev.png');
pdfPath = fullfile(outputDir, 'counterdiabatic_passage_dev.pdf');
exportgraphics(fig, pngPath, 'Resolution', 300);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');

sidebandReport = sideband_validity(model, controls, linspace(0, T, 101), 6);
jsonCases = rmfield(results, {'obs', 'targetConditionalFidelity'});

summary = struct();
summary.runMode = runMode;
summary.model = "central-sideband effective finite-mode dynamics";
summary.counterdiabaticHamiltonian = "H_CD = thetaDot(t) sigma_y in the atomic subspace";
summary.Nc = model.Nc;
summary.dimension = model.dimension;
summary.T = T;
summary.nu = nu;
summary.u0 = u0;
summary.rtol = 1e-9;
summary.atol = 1e-11;
summary.maxStep = min(T/1000, 0.05);
summary.maxOmegaCDOverXi = max(abs(omegaCD / params.xi));
summary.sideband = sidebandReport;
summary.cases = jsonCases;
summary.outputs = struct('png', string(pngPath), 'pdf', string(pdfPath));

matPath = fullfile(dataDir, 'counterdiabatic_passage_dev.mat');
jsonPath = fullfile(dataDir, 'counterdiabatic_passage_dev_summary.json');
save(matPath, 'tOut', 'psiAll', 'results', 'caseSpecs', 'target', ...
    'targetConcurrence', 'u1', 'u2', 'omegaCD', 'summary', 'params', 'u0', 'T', 'nu');
write_json(jsonPath, summary);

fprintf('Counterdiabatic passage outputs written to:\n  %s\n  %s\n  %s\n  %s\n', ...
    matPath, jsonPath, pngPath, pdfPath);
end

function result = empty_result()
result.name = "";
result.counterdiabaticScale = [];
result.runtimeSeconds = [];
result.maxNormError = [];
result.finalAtomicPopulation = [];
result.finalPhotonicPopulation = [];
result.finalConcurrence = [];
result.finalConditionalConcurrence = [];
result.finalBellPlusConditionalFidelity = [];
result.finalTargetConditionalFidelity = [];
result.minimumTargetConditionalFidelity = [];
result.obs = [];
result.targetConditionalFidelity = [];
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

function hide_axes_toolbar(ax)
if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end
end
