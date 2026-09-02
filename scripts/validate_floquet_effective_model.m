function summary = validate_floquet_effective_model(runMode)
%VALIDATE_FLOQUET_EFFECTIVE_MODEL Compare lab-frame and effective dynamics.
%
% The validation isolates the Floquet passage stage used in the two-stage
% Bell-BIC protocol. Starting from |eg,0>, it propagates:
%
%   1. the microscopic laboratory-frame Hamiltonian with the full time-
%      dependent atomic frequencies Omega_i(t), and
%   2. the central-sideband effective Hamiltonian with couplings
%      g_ik^eff(t)=J0[beta_i(t)] g_ik.
%
% The laboratory state is compared after removing the local micromotion
% phases c_i -> exp[i beta_i(t) sin(nu t)] c_i. The same CD term is included
% in the slow frame and represented in the laboratory frame with the
% corresponding micromotion phases.

if nargin < 1
    runMode = "dev";
else
    runMode = string(runMode);
end

if ~ismember(runMode, ["dev", "paper"])
    error('runMode must be "dev" or "paper".');
end

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
addpath(fullfile(projectRoot, 'src', 'matlab'));

if runMode == "paper"
    dataDir = fullfile(projectRoot, 'data', 'paper');
    outputStem = 'floquet_effective_validation_paper';
    figureStem = 'Figure_Floquet_Effective_Validation';
    params.Nc = 4 * 501;
    numTimes = 2401;
else
    dataDir = fullfile(projectRoot, 'data', 'development');
    outputStem = 'floquet_effective_validation_dev';
    figureStem = 'floquet_effective_validation_dev';
    params.Nc = 4 * 51;
    numTimes = 1201;
end

figureDir = fullfile(projectRoot, 'figures');
outputDir = fullfile(projectRoot, 'outputs', 'floquet');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

params.xi = 1.0;
params.wc = 0.0;
params.g = 0.1;
params.n1 = 6;
params.n2 = 6;
params.x1 = 0;
params.Dx = 2;
params.Omega0 = params.wc;

u0 = 0.9;
Tpassage = 20 / params.xi;
nu = 8 * params.xi;
cdScale = 1;

model = giant_atom_static_model(params);
controls = floquet_controls(u0, Tpassage, nu);
tEval = linspace(0, Tpassage, numTimes);

psi0 = zeros(model.dimension, 1);
psi0(1) = 1;

rhsOptions.counterdiabaticScale = cdScale;
effectiveSolverOptions = odeset( ...
    'RelTol', 1e-10, ...
    'AbsTol', 1e-12, ...
    'MaxStep', 0.02);
labSolverOptions = odeset( ...
    'RelTol', 1e-9, ...
    'AbsTol', 1e-11, ...
    'MaxStep', min(Tpassage / 4000, 2*pi / nu / 80));

fprintf('Running effective central-sideband passage with Nc=%d...\n', model.Nc);
timer = tic;
[tEff, psiEffRows] = ode113(@(t, psi) floquet_effective_rhs(t, psi, model, controls, rhsOptions), ...
    tEval, psi0, effectiveSolverOptions);
effectiveRuntimeSeconds = toc(timer);

fprintf('Running exact laboratory-frame passage with Nc=%d and nu/xi=%.2f...\n', model.Nc, nu / params.xi);
timer = tic;
[tLab, psiLabRows] = ode113(@(t, psi) floquet_lab_rhs(t, psi, model, controls, rhsOptions), ...
    tEval, psi0, labSolverOptions);
labRuntimeSeconds = toc(timer);

if max(abs(tEff - tLab)) > 1e-12
    error('Effective and laboratory integrations returned different time grids.');
end

psiEff = transpose(psiEffRows);
psiLab = transpose(psiLabRows);
psiLabSlow = micromotion_to_rotating_frame(tLab.', psiLab, controls);

obsEff = observables_single_excitation(psiEff);
obsLabRaw = observables_single_excitation(psiLab);
obsLabSlow = observables_single_excitation(psiLabSlow);

theta = controls.theta(tEff.');
targetConcurrence = sin(2 * theta);

stateOverlap = sum(conj(psiEff) .* psiLabSlow, 1);
stateInfidelity = max(0, 1 - abs(stateOverlap).^2);
deltaConcurrence = abs(obsEff.concurrence - obsLabSlow.concurrence);
deltaBellFidelity = abs(obsEff.bellPlusFidelity - obsLabSlow.bellPlusFidelity);
deltaAtomicPopulation = abs(obsEff.atomicPopulation - obsLabSlow.atomicPopulation);

sidebandReport = sideband_validity(model, controls, linspace(0, Tpassage, 201), 6);

summary = struct();
summary.runMode = runMode;
summary.model = "full laboratory-frame finite-mode dynamics versus central-sideband effective dynamics";
summary.validationStage = "Floquet Bell-BIC passage after local |eg> preparation";
summary.Nc = model.Nc;
summary.dimension = model.dimension;
summary.kConvention = model.kConvention;
summary.xi = params.xi;
summary.gOverXi = params.g / params.xi;
summary.nuOverXi = nu / params.xi;
summary.u0 = u0;
summary.TpassageXi = params.xi * Tpassage;
summary.counterdiabaticScale = cdScale;
summary.effectiveRuntimeSeconds = effectiveRuntimeSeconds;
summary.labRuntimeSeconds = labRuntimeSeconds;
summary.maxNormErrorEffective = max(abs(obsEff.norm - 1));
summary.maxNormErrorLab = max(abs(obsLabRaw.norm - 1));
summary.finalConcurrenceEffective = obsEff.concurrence(end);
summary.finalConcurrenceLabSlow = obsLabSlow.concurrence(end);
summary.finalBellFidelityEffective = obsEff.bellPlusFidelity(end);
summary.finalBellFidelityLabSlow = obsLabSlow.bellPlusFidelity(end);
summary.finalBellFidelityLabRaw = obsLabRaw.bellPlusFidelity(end);
summary.finalAtomicPopulationEffective = obsEff.atomicPopulation(end);
summary.finalAtomicPopulationLabSlow = obsLabSlow.atomicPopulation(end);
summary.maxAbsDeltaConcurrence = max(deltaConcurrence);
summary.maxAbsDeltaBellFidelity = max(deltaBellFidelity);
summary.maxAbsDeltaAtomicPopulation = max(deltaAtomicPopulation);
summary.maxStateInfidelitySlowFrame = max(stateInfidelity);
summary.finalStateInfidelitySlowFrame = stateInfidelity(end);
summary.sidebandReport = sidebandReport;

sourceData = make_source_table(tEff, params.xi, targetConcurrence, obsEff, ...
    obsLabRaw, obsLabSlow, deltaConcurrence, deltaBellFidelity, ...
    deltaAtomicPopulation, stateInfidelity);

fig = make_validation_figure(tEff, params.xi, targetConcurrence, obsEff, ...
    obsLabRaw, obsLabSlow, deltaConcurrence, deltaBellFidelity, ...
    stateInfidelity);

if runMode == "paper"
    pdfPath = fullfile(figureDir, [figureStem '.pdf']);
    pngPath = fullfile(figureDir, [figureStem '.png']);
    tiffPath = fullfile(figureDir, [figureStem '.tif']);
else
    pdfPath = fullfile(outputDir, [figureStem '.pdf']);
    pngPath = fullfile(outputDir, [figureStem '.png']);
    tiffPath = fullfile(outputDir, [figureStem '.tif']);
end

exportgraphics(fig, pdfPath, 'ContentType', 'vector');
exportgraphics(fig, pngPath, 'Resolution', 300);
exportgraphics(fig, tiffPath, 'Resolution', 600);

matPath = fullfile(dataDir, [outputStem '.mat']);
jsonPath = fullfile(dataDir, [outputStem '_summary.json']);
csvPath = fullfile(dataDir, [outputStem '_source_data.csv']);

if runMode == "paper"
    save(matPath, 'tEff', 'obsEff', 'obsLabRaw', 'obsLabSlow', ...
        'targetConcurrence', 'deltaConcurrence', ...
        'deltaBellFidelity', 'deltaAtomicPopulation', 'stateInfidelity', ...
        'summary', 'params', 'u0', 'Tpassage', 'nu');
else
    save(matPath, 'tEff', 'psiEff', 'psiLab', 'psiLabSlow', 'obsEff', ...
        'obsLabRaw', 'obsLabSlow', 'targetConcurrence', ...
        'deltaConcurrence', 'deltaBellFidelity', 'deltaAtomicPopulation', ...
        'stateInfidelity', 'summary', 'params', 'u0', 'Tpassage', 'nu');
end
write_json(jsonPath, summary);
writetable(sourceData, csvPath);

summary.outputs = struct( ...
    'pdf', string(pdfPath), ...
    'png', string(pngPath), ...
    'tiff', string(tiffPath), ...
    'mat', string(matPath), ...
    'json', string(jsonPath), ...
    'csv', string(csvPath));
write_json(jsonPath, summary);

fprintf(['Validation complete.\n' ...
    '  max |Delta C|       = %.3e\n' ...
    '  max |Delta F|       = %.3e\n' ...
    '  max state infid.    = %.3e\n' ...
    '  final F_eff/F_lab   = %.6f / %.6f\n' ...
    'Outputs:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n'], ...
    summary.maxAbsDeltaConcurrence, ...
    summary.maxAbsDeltaBellFidelity, ...
    summary.maxStateInfidelitySlowFrame, ...
    summary.finalBellFidelityEffective, ...
    summary.finalBellFidelityLabSlow, ...
    pdfPath, pngPath, tiffPath, matPath, jsonPath, csvPath);
end

function fig = make_validation_figure(t, xi, targetConcurrence, obsEff, ...
    obsLabRaw, obsLabSlow, deltaConcurrence, deltaBellFidelity, ...
    stateInfidelity)

tXi = xi * t(:).';
colors.effective = [0.0000 0.4470 0.7410];
colors.lab = [0.8500 0.3250 0.0980];
colors.raw = [0.55 0.55 0.55];
colors.target = [0.15 0.15 0.15];

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 16]);
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

axA = nexttile; hold(axA, 'on'); box(axA, 'on'); grid(axA, 'on');
plot(axA, tXi, targetConcurrence, ':', 'Color', colors.target, 'LineWidth', 2.0);
plot(axA, tXi, obsEff.concurrence, '-', 'Color', colors.effective, 'LineWidth', 2.1);
plot(axA, tXi, obsLabSlow.concurrence, '--', 'Color', colors.lab, 'LineWidth', 1.9);
xlabel(axA, '$\xi t$', 'Interpreter', 'latex');
ylabel(axA, '$\mathcal{C}(t)$', 'Interpreter', 'latex');
ylim(axA, [-0.03 1.03]);
title(axA, 'Concurrence', 'Interpreter', 'latex');
legend(axA, {'target', '$H_{\rm eff}$', 'lab, slow frame'}, ...
    'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 11);
format_axis(axA);
panel_label(axA, 'a');

axB = nexttile; hold(axB, 'on'); box(axB, 'on'); grid(axB, 'on');
plot(axB, tXi, obsLabRaw.bellPlusFidelity, '-', 'Color', colors.raw, ...
    'LineWidth', 1.0);
plot(axB, tXi, obsEff.bellPlusFidelity, '-', 'Color', colors.effective, ...
    'LineWidth', 2.1);
plot(axB, tXi, obsLabSlow.bellPlusFidelity, '--', 'Color', colors.lab, ...
    'LineWidth', 1.9);
xlabel(axB, '$\xi t$', 'Interpreter', 'latex');
ylabel(axB, '$F_{\Psi^+}(t)$', 'Interpreter', 'latex');
ylim(axB, [-0.03 1.03]);
title(axB, 'Bell fidelity', 'Interpreter', 'latex');
legend(axB, {'lab frame', '$H_{\rm eff}$', 'lab, slow frame'}, ...
    'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 11);
format_axis(axB);
panel_label(axB, 'b');

axC = nexttile; hold(axC, 'on'); box(axC, 'on'); grid(axC, 'on');
plot(axC, tXi, obsEff.atomicPopulation, '-', 'Color', colors.effective, ...
    'LineWidth', 2.1);
plot(axC, tXi, obsLabSlow.atomicPopulation, '--', 'Color', colors.lab, ...
    'LineWidth', 1.9);
plot(axC, tXi, obsEff.photonicPopulation, ':', 'Color', colors.effective, ...
    'LineWidth', 1.8);
plot(axC, tXi, obsLabSlow.photonicPopulation, '-.', 'Color', colors.lab, ...
    'LineWidth', 1.6);
xlabel(axC, '$\xi t$', 'Interpreter', 'latex');
ylabel(axC, 'population', 'Interpreter', 'latex');
ylim(axC, [-0.03 1.03]);
title(axC, 'Atomic and photonic population', 'Interpreter', 'latex');
legend(axC, {'$P_{\rm at}^{\rm eff}$', '$P_{\rm at}^{\rm lab}$', ...
    '$P_{\gamma}^{\rm eff}$', '$P_{\gamma}^{\rm lab}$'}, ...
    'Interpreter', 'latex', 'Location', 'east', 'FontSize', 10);
format_axis(axC);
panel_label(axC, 'c');

axD = nexttile; hold(axD, 'on'); box(axD, 'on'); grid(axD, 'on');
errorScale = 1e3;
plot(axD, tXi, errorScale * deltaConcurrence, '-', 'Color', colors.effective, ...
    'LineWidth', 1.9);
plot(axD, tXi, errorScale * deltaBellFidelity, '--', 'Color', colors.lab, ...
    'LineWidth', 1.9);
plot(axD, tXi, errorScale * stateInfidelity, ':', 'Color', colors.target, ...
    'LineWidth', 2.0);
xlabel(axD, '$\xi t$', 'Interpreter', 'latex');
ylabel(axD, '$10^3\times$ error', 'Interpreter', 'latex');
title(axD, 'Effective-model error', 'Interpreter', 'latex');
legend(axD, {'$|\Delta\mathcal{C}|$', '$|\Delta F_{\Psi^+}|$', ...
    '$1-|\langle\psi_{\rm eff}|\psi_{\rm lab}\rangle|^2$'}, ...
    'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 10);
maxScaledError = errorScale * max([deltaConcurrence, deltaBellFidelity, stateInfidelity], [], 'all');
ylim(axD, [0 1.20 * max(maxScaledError, 1e-3)]);
format_axis(axD);
panel_label(axD, 'd');

end

function format_axis(ax)
set(ax, 'FontSize', 13, 'TickLabelInterpreter', 'latex', ...
    'LineWidth', 0.8, 'Layer', 'top');
xlim(ax, [0 max(ax.XLim)]);
if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end
end

function panel_label(ax, label)
text(ax, 0.03, 0.93, ['\textbf{' label '}'], 'Units', 'normalized', ...
    'Interpreter', 'latex', 'FontSize', 13, 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');
end

function sourceData = make_source_table(t, xi, targetConcurrence, obsEff, ...
    obsLabRaw, obsLabSlow, deltaConcurrence, deltaBellFidelity, ...
    deltaAtomicPopulation, stateInfidelity)

sourceData = table();
sourceData.xi_t = xi * t(:);
sourceData.target_concurrence = targetConcurrence(:);
sourceData.C_effective = obsEff.concurrence(:);
sourceData.C_lab_slow_frame = obsLabSlow.concurrence(:);
sourceData.F_Bell_effective = obsEff.bellPlusFidelity(:);
sourceData.F_Bell_lab_slow_frame = obsLabSlow.bellPlusFidelity(:);
sourceData.F_Bell_lab_frame = obsLabRaw.bellPlusFidelity(:);
sourceData.P_atomic_effective = obsEff.atomicPopulation(:);
sourceData.P_atomic_lab_slow_frame = obsLabSlow.atomicPopulation(:);
sourceData.P_photon_effective = obsEff.photonicPopulation(:);
sourceData.P_photon_lab_slow_frame = obsLabSlow.photonicPopulation(:);
sourceData.abs_delta_concurrence = deltaConcurrence(:);
sourceData.abs_delta_bell_fidelity = deltaBellFidelity(:);
sourceData.abs_delta_atomic_population = deltaAtomicPopulation(:);
sourceData.state_infidelity_slow_frame = stateInfidelity(:);
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
