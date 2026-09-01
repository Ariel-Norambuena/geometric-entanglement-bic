function summary = run_two_stage_bell_bic_protocol(runMode)
%RUN_TWO_STAGE_BELL_BIC_PROTOCOL Full preparation-plus-Bell-BIC diagnostic.
%
% The state is propagated in one continuous finite-mode dynamics:
%
%   Stage I:  |gg,0> -> |eg,0> with a resonant local pi pulse on atom 1.
%   Stage II: |eg,0> -> Bell-BIC atomic component with Floquet controls and
%            the optional atomic counterdiabatic term.
%
% This is a development figure. It demonstrates the control principle in the
% effective central-sideband model; laboratory-frame and dressed-BIC
% validation remain separate checks.

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
Tprep = 2 / params.xi;
Tpassage = 20 / params.xi;
nu = 8 * params.xi;
drivePhase = pi / 2;

model = giant_atom_static_model(params);
controls = floquet_controls(u0, Tpassage, nu);
protocol = struct( ...
    'Tprep', Tprep, ...
    'Tpassage', Tpassage, ...
    'totalTime', Tprep + Tpassage, ...
    'drivePhase', drivePhase, ...
    'piPulseArea', pi);

tPrep = linspace(0, Tprep, 301);
tPassage = Tprep + linspace(0, Tpassage, 1501);
tEval = unique([tPrep tPassage]);

caseSpecs = struct( ...
    'name', {'without_cd', 'with_cd'}, ...
    'plotLabel', {'without CD', 'with CD'}, ...
    'counterdiabaticScale', {0, 1});

results = repmat(empty_result(), numel(caseSpecs), 1);
psiAll = cell(numel(caseSpecs), 1);

for idx = 1:numel(caseSpecs)
    psi0 = zeros(model.dimension + 1, 1);
    psi0(1) = 1;

    rhsOptions.counterdiabaticScale = caseSpecs(idx).counterdiabaticScale;
    solverOptions = odeset( ...
        'RelTol', 1e-9, ...
        'AbsTol', 1e-11, ...
        'MaxStep', min(Tprep/300, 0.02));

    timer = tic;
    [tOut, psiOutRows] = ode113(@(t, psi) two_stage_rhs(t, psi, model, controls, protocol, rhsOptions), ...
        tEval, psi0, solverOptions);
    runtimeSeconds = toc(timer);

    psi = transpose(psiOutRows);
    obs = two_stage_observables(psi);
    target = passage_target(tOut.', controls, protocol);
    targetOverlap = abs(conj(target(1, :)) .* psi(2, :) + conj(target(2, :)) .* psi(3, :)).^2;
    targetConditionalFidelity = targetOverlap ./ obs.atomicSinglePopulation;

    [~, prepIdx] = min(abs(tOut - Tprep));

    results(idx).name = caseSpecs(idx).name;
    results(idx).counterdiabaticScale = caseSpecs(idx).counterdiabaticScale;
    results(idx).runtimeSeconds = runtimeSeconds;
    results(idx).maxNormError = max(abs(obs.norm - 1));
    results(idx).afterPrepGroundPopulation = obs.groundPopulation(prepIdx);
    results(idx).afterPrepEgPopulation = obs.egPopulation(prepIdx);
    results(idx).afterPrepGePopulation = obs.gePopulation(prepIdx);
    results(idx).afterPrepPhotonicPopulation = obs.photonicPopulation(prepIdx);
    results(idx).finalGroundPopulation = obs.groundPopulation(end);
    results(idx).finalAtomicSinglePopulation = obs.atomicSinglePopulation(end);
    results(idx).finalPhotonicPopulation = obs.photonicPopulation(end);
    results(idx).finalConcurrence = obs.concurrence(end);
    results(idx).finalConditionalConcurrence = obs.conditionalConcurrence(end);
    results(idx).finalBellPlusFidelity = obs.bellPlusFidelity(end);
    results(idx).finalConditionalBellPlusFidelity = obs.conditionalBellPlusFidelity(end);
    results(idx).finalTargetConditionalFidelity = targetConditionalFidelity(end);
    results(idx).minimumPassageTargetConditionalFidelity = min( ...
        targetConditionalFidelity(tOut.' >= Tprep));
    results(idx).obs = obs;
    results(idx).targetConditionalFidelity = targetConditionalFidelity;

    psiAll{idx} = psi;

    fprintf(['%-10s | P_eg(after prep)=%.6f | C_final=%.6f | ' ...
        'F_Bell_final=%.6f | F_Bell_cond=%.6f | max norm err=%.2e\n'], ...
        caseSpecs(idx).plotLabel, ...
        results(idx).afterPrepEgPopulation, ...
        results(idx).finalConcurrence, ...
        results(idx).finalBellPlusFidelity, ...
        results(idx).finalConditionalBellPlusFidelity, ...
        results(idx).maxNormError);
end

tXi = params.xi * tOut;
totalTimeXi = params.xi * protocol.totalTime;
thetaGlobal = global_theta(tOut.', controls, protocol);
targetConcurrence = sin(2 * thetaGlobal);
targetConcurrence(tOut.' < Tprep) = NaN;
[u1Global, u2Global] = global_u(tOut.', controls, protocol, u0);
omegaCD = global_theta_dot(tOut.', controls, protocol);
omegaPi = pi_pulse_envelope(tOut.', protocol);

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 1 18 27]);
tiledlayout(fig, 5, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on; box on; grid on;
plot(tXi, results(2).obs.groundPopulation, 'LineWidth', 1.8);
plot(tXi, results(2).obs.egPopulation, 'LineWidth', 1.8);
plot(tXi, results(2).obs.gePopulation, 'LineWidth', 1.8);
plot(tXi, results(2).obs.photonicPopulation, 'LineWidth', 1.8);
xlim([0 totalTimeXi]);
mark_stage_boundary(Tprep, params.xi);
ylabel('population', 'Interpreter', 'latex');
ylim([-0.03 1.03]);
legend({'$P_{gg}$', '$P_{eg}$', '$P_{ge}$', '$P_{\gamma}$'}, ...
    'Interpreter', 'latex', 'Location', 'east');
set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');
hide_axes_toolbar(gca);

nexttile; hold on; box on; grid on;
plot(tXi, targetConcurrence, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 2.0);
plot(tXi, results(1).obs.concurrence, 'LineWidth', 1.8);
plot(tXi, results(2).obs.concurrence, 'LineWidth', 2.2);
xlim([0 totalTimeXi]);
mark_stage_boundary(Tprep, params.xi);
ylabel('$\mathcal{C}(t)$', 'Interpreter', 'latex');
ylim([-0.03 1.03]);
legend({'target passage', 'without CD', 'with CD'}, ...
    'Interpreter', 'latex', 'Location', 'southeast');
set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');
hide_axes_toolbar(gca);

nexttile; hold on; box on; grid on;
plot(tXi, results(1).obs.bellPlusFidelity, 'LineWidth', 1.8);
plot(tXi, results(2).obs.bellPlusFidelity, 'LineWidth', 2.2);
plot(tXi, results(2).obs.conditionalBellPlusFidelity, '--', 'LineWidth', 1.8);
xlim([0 totalTimeXi]);
mark_stage_boundary(Tprep, params.xi);
ylabel('$F_{\Psi^+}(t)$', 'Interpreter', 'latex');
ylim([-0.03 1.03]);
legend({'without CD', 'with CD', 'with CD, conditional'}, ...
    'Interpreter', 'latex', 'Location', 'southeast');
set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');
hide_axes_toolbar(gca);

nexttile; hold on; box on; grid on;
plot(tXi, u1Global, 'LineWidth', 2.0);
plot(tXi, u2Global, '--', 'LineWidth', 2.0);
xlim([0 totalTimeXi]);
mark_stage_boundary(Tprep, params.xi);
ylabel('$u_i(t)$', 'Interpreter', 'latex');
ylim([-0.03 1.03]);
legend({'$u_1(t)$', '$u_2(t)$'}, 'Interpreter', 'latex', 'Location', 'east');
set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');
hide_axes_toolbar(gca);

nexttile; hold on; box on; grid on;
yyaxis left;
plot(tXi, omegaPi / params.xi, 'LineWidth', 2.0);
ylabel('$\Omega_{\pi}(t)/\xi$', 'Interpreter', 'latex');
ylim([-0.05 1.15 * max(omegaPi / params.xi)]);
yyaxis right;
plot(tXi, omegaCD / params.xi, 'k--', 'LineWidth', 2.0);
ylabel('$\Omega_{\rm CD}(t)/\xi$', 'Interpreter', 'latex');
ylim([-0.005 1.15 * max(omegaCD / params.xi)]);
xlim([0 totalTimeXi]);
mark_stage_boundary(Tprep, params.xi);
xlabel('$\xi t$', 'Interpreter', 'latex');
legend({'$\Omega_{\pi}(t)$', '$\Omega_{\rm CD}(t)$'}, ...
    'Interpreter', 'latex', 'Location', 'northeast');
set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');
hide_axes_toolbar(gca);

pngPath = fullfile(outputDir, 'two_stage_bell_bic_protocol_dev.png');
pdfPath = fullfile(outputDir, 'two_stage_bell_bic_protocol_dev.pdf');
exportgraphics(fig, pngPath, 'Resolution', 300);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');

summaryCases = rmfield(results, {'obs', 'targetConditionalFidelity'});
summary = struct();
summary.runMode = runMode;
summary.model = "two-stage central-sideband effective finite-mode dynamics";
summary.stageI = "local resonant pi pulse on atom 1, with u1=0 and u2=u0";
summary.stageII = "Floquet dark-state passage with optional H_CD = thetaDot(t) sigma_y";
summary.Nc = model.Nc;
summary.dimension = model.dimension + 1;
summary.Tprep = Tprep;
summary.Tpassage = Tpassage;
summary.totalTime = protocol.totalTime;
summary.nu = nu;
summary.u0 = u0;
summary.drivePhase = drivePhase;
summary.maxOmegaPiOverXi = max(omegaPi / params.xi);
summary.maxOmegaCDOverXi = max(omegaCD / params.xi);
summary.rtol = 1e-9;
summary.atol = 1e-11;
summary.cases = summaryCases;
summary.outputs = struct('png', string(pngPath), 'pdf', string(pdfPath));

matPath = fullfile(dataDir, 'two_stage_bell_bic_protocol_dev.mat');
jsonPath = fullfile(dataDir, 'two_stage_bell_bic_protocol_dev_summary.json');
save(matPath, 'tOut', 'psiAll', 'results', 'caseSpecs', 'u1Global', ...
    'u2Global', 'omegaPi', 'omegaCD', 'targetConcurrence', 'summary', ...
    'params', 'u0', 'Tprep', 'Tpassage', 'nu');
write_json(jsonPath, summary);

fprintf('Two-stage protocol outputs written to:\n  %s\n  %s\n  %s\n  %s\n', ...
    matPath, jsonPath, pngPath, pdfPath);
end

function dpsi = two_stage_rhs(t, psi, model, controls, protocol, options)
c0 = psi(1);
c1 = psi(2);
c2 = psi(3);
phi = psi(4:end);

[u1, u2] = instantaneous_u(t, controls, protocol);
thetaDot = instantaneous_theta_dot(t, controls, protocol);
omegaPi = pi_pulse_envelope(t, protocol);
omegaCD = options.counterdiabaticScale * thetaDot;

f1k = u1 * model.g1k;
f2k = u2 * model.g2k;

dc0 = -1i * (omegaPi / 2) * exp(-1i * protocol.drivePhase) * c1;
dc1 = -1i * (model.Omega0 * c1 + transpose(f1k) * phi ...
    + (omegaPi / 2) * exp(1i * protocol.drivePhase) * c0);
dc2 = -1i * (model.Omega0 * c2 + transpose(f2k) * phi);
dphi = -1i * (model.wk .* phi + conj(f1k) * c1 + conj(f2k) * c2);

dc1 = dc1 - omegaCD * c2;
dc2 = dc2 + omegaCD * c1;

dpsi = [dc0; dc1; dc2; dphi];
end

function obs = two_stage_observables(psi)
c0 = psi(1, :);
c1 = psi(2, :);
c2 = psi(3, :);
phi = psi(4:end, :);

obs.groundPopulation = abs(c0).^2;
obs.egPopulation = abs(c1).^2;
obs.gePopulation = abs(c2).^2;
obs.atomicSinglePopulation = obs.egPopulation + obs.gePopulation;
obs.photonicPopulation = sum(abs(phi).^2, 1);
obs.norm = obs.groundPopulation + obs.atomicSinglePopulation + obs.photonicPopulation;
obs.concurrence = 2 * abs(c1 .* c2);
obs.conditionalConcurrence = obs.concurrence ./ obs.atomicSinglePopulation;
obs.bellPlusFidelity = abs(c1 + c2).^2 / 2;
obs.conditionalBellPlusFidelity = obs.bellPlusFidelity ./ obs.atomicSinglePopulation;
end

function [u1, u2] = instantaneous_u(t, controls, protocol)
if t < protocol.Tprep
    u1 = 0;
    u2 = controls.u0;
else
    [u1, u2] = controls.u(t - protocol.Tprep);
end
end

function thetaDot = instantaneous_theta_dot(t, controls, protocol)
if t < protocol.Tprep
    thetaDot = 0;
else
    thetaDot = controls.thetaDot(t - protocol.Tprep);
end
end

function omegaPi = pi_pulse_envelope(t, protocol)
inside = (t >= 0) & (t <= protocol.Tprep);
s = min(1, max(0, t ./ protocol.Tprep));
omegaPeak = 2 * protocol.piPulseArea / protocol.Tprep;
omegaPi = omegaPeak .* sin(pi * s).^2 .* inside;
end

function theta = global_theta(t, controls, protocol)
tau = min(protocol.Tpassage, max(0, t - protocol.Tprep));
theta = controls.theta(tau);
theta(t < protocol.Tprep) = 0;
end

function [u1, u2] = global_u(t, controls, protocol, u0)
tau = min(protocol.Tpassage, max(0, t - protocol.Tprep));
[u1, u2] = controls.u(tau);
u1(t < protocol.Tprep) = 0;
u2(t < protocol.Tprep) = u0;
end

function thetaDot = global_theta_dot(t, controls, protocol)
tau = min(protocol.Tpassage, max(0, t - protocol.Tprep));
thetaDot = controls.thetaDot(tau);
thetaDot(t < protocol.Tprep) = 0;
end

function target = passage_target(t, controls, protocol)
theta = global_theta(t, controls, protocol);
target = [cos(theta); sin(theta)];
target(:, t < protocol.Tprep) = NaN;
end

function result = empty_result()
result.name = "";
result.counterdiabaticScale = [];
result.runtimeSeconds = [];
result.maxNormError = [];
result.afterPrepGroundPopulation = [];
result.afterPrepEgPopulation = [];
result.afterPrepGePopulation = [];
result.afterPrepPhotonicPopulation = [];
result.finalGroundPopulation = [];
result.finalAtomicSinglePopulation = [];
result.finalPhotonicPopulation = [];
result.finalConcurrence = [];
result.finalConditionalConcurrence = [];
result.finalBellPlusFidelity = [];
result.finalConditionalBellPlusFidelity = [];
result.finalTargetConditionalFidelity = [];
result.minimumPassageTargetConditionalFidelity = [];
result.obs = [];
result.targetConditionalFidelity = [];
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

function mark_stage_boundary(Tprep, xi)
xline(xi * Tprep, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
yl = ylim;
text(0.18 * xi * Tprep, yl(2) - 0.12 * diff(yl), 'prep', ...
    'Interpreter', 'latex', 'FontSize', 11);
text(xi * Tprep + 0.04 * diff(xlim), yl(2) - 0.12 * diff(yl), 'Bell-BIC passage', ...
    'Interpreter', 'latex', 'FontSize', 11);
end

function hide_axes_toolbar(ax)
if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end
end
