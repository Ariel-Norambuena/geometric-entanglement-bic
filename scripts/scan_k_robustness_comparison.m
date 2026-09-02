function summary = scan_k_robustness_comparison(runMode)
%SCAN_K_ROBUSTNESS_COMPARISON Compare k-sensitivity benchmarks.
%
% The first benchmark reproduces the spirit of Effect_dk_Bell_states.m from
% the PRA code package: a passive Markovian Bell-BIC state is initialized at
% K = pi/2 + dk and then compared with the fixed Bell state |Psi+>.
%
% The second benchmark scans the two-stage Floquet-CD protocol under the
% same fractional resonant-wave-vector mismatch. For a detuned K, the atomic
% frequency is Omega(K)=wc-2 xi cos(K), and the dynamics are simulated in the
% resonant atomic rotating frame by shifting the photonic detunings.

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
    outputStem = 'k_robustness_comparison_paper';
    params.Nc = 4 * 501;
    deltaFractions = [-0.10 -0.075 -0.05 -0.025 -0.01 0 0.01 0.025 0.05 0.075 0.10];
    durationGridXi = [5 10 15 20 30 40];
    praTime = linspace(0, 500, 1601);
else
    dataDir = fullfile(projectRoot, 'data', 'development');
    outputStem = 'k_robustness_comparison_dev';
    params.Nc = 4 * 51;
    deltaFractions = [-0.10 -0.05 -0.01 0 0.01 0.05 0.10];
    durationGridXi = [5 10 20];
    praTime = linspace(0, 500, 801);
end

figureDir = fullfile(projectRoot, 'figures');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

params.xi = 1.0;
params.wc = 0.0;
params.g = 0.1;
params.n1 = 6;
params.n2 = 6;
params.x1 = 0;
params.Dx = 2;
params.Omega0 = 0;

k0 = pi / 2;
u0 = 0.9;
nu = 8 * params.xi;
Tprep = 2 / params.xi;
Treference = 20 / params.xi;
drivePhase = pi / 2;

praDeltas = [0.001 0.005 0.01 0.05 0.10];
pra = simulate_pra_markov_sensitivity(praDeltas, praTime);

numDelta = numel(deltaFractions);
numDurations = numel(durationGridXi);
withoutCD = repmat(empty_protocol_result(), numDelta, 1);
withCDReference = repmat(empty_protocol_result(), numDelta, 1);
fidelityGrid = zeros(numDurations, numDelta);
concurrenceGrid = zeros(numDurations, numDelta);
atomicPopulationGrid = zeros(numDurations, numDelta);

fprintf('Scanning Floquet-CD robustness with Nc=%d...\n', params.Nc);
for id = 1:numDelta
    delta = deltaFractions(id);
    K = k0 * (1 + delta);

    withoutCD(id) = simulate_two_stage_detuned(params, K, u0, nu, Tprep, ...
        Treference, drivePhase, 0);
    withCDReference(id) = simulate_two_stage_detuned(params, K, u0, nu, Tprep, ...
        Treference, drivePhase, 1);

    for it = 1:numDurations
        result = simulate_two_stage_detuned(params, K, u0, nu, Tprep, ...
            durationGridXi(it) / params.xi, drivePhase, 1);
        fidelityGrid(it, id) = result.finalBellPlusFidelity;
        concurrenceGrid(it, id) = result.finalConcurrence;
        atomicPopulationGrid(it, id) = result.finalAtomicSinglePopulation;
    end

    fprintf(['  deltaK/K0=% .3f | no CD F=%.5f | CD F=%.5f | ' ...
        'CD C=%.5f | P_at=%.5f\n'], ...
        delta, withoutCD(id).finalBellPlusFidelity, ...
        withCDReference(id).finalBellPlusFidelity, ...
        withCDReference(id).finalConcurrence, ...
        withCDReference(id).finalAtomicSinglePopulation);
end

summary = struct();
summary.runMode = runMode;
summary.description = "PRA-style K sensitivity versus two-stage Floquet-CD robustness";
summary.k0 = k0;
summary.deltaFractions = deltaFractions;
summary.durationGridXi = durationGridXi;
summary.TprepXi = params.xi * Tprep;
summary.TreferenceXi = params.xi * Treference;
summary.u0 = u0;
summary.nuOverXi = nu / params.xi;
summary.params = params;
summary.pra = pra;
summary.withoutCD = withoutCD;
summary.withCDReference = withCDReference;
summary.fidelityGrid = fidelityGrid;
summary.concurrenceGrid = concurrenceGrid;
summary.atomicPopulationGrid = atomicPopulationGrid;

sourceData = make_source_data(deltaFractions, durationGridXi, withoutCD, ...
    withCDReference, fidelityGrid, concurrenceGrid, atomicPopulationGrid);
praCurveData = make_pra_curve_data(pra);

matPath = fullfile(dataDir, [outputStem '.mat']);
jsonPath = fullfile(dataDir, [outputStem '_summary.json']);
csvPath = fullfile(dataDir, [outputStem '_source_data.csv']);
praCsvPath = fullfile(dataDir, [outputStem '_pra_fidelity_curves.csv']);
save(matPath, 'summary', 'sourceData', 'pra', 'withoutCD', 'withCDReference', ...
    'fidelityGrid', 'concurrenceGrid', 'atomicPopulationGrid', 'praCurveData');
write_json(jsonPath, compact_json_summary(summary));
writetable(sourceData, csvPath);
writetable(praCurveData, praCsvPath);

paths = make_robustness_figure(summary, figureDir);
summary.outputs = paths;
summary.sourceData = struct('summaryCsv', string(csvPath), 'praCurvesCsv', string(praCsvPath));
write_json(jsonPath, compact_json_summary(summary));

fprintf('K-robustness comparison written to:\n  %s\n  %s\n  %s\n  %s\n', ...
    matPath, jsonPath, csvPath, praCsvPath);
fprintf('Figure files:\n  %s\n  %s\n  %s\n', paths.pdf, paths.png, paths.tiff);
end

function pra = simulate_pra_markov_sensitivity(deltaFractions, t)
Natoms = 2;
dimT = 2^Natoms;
Is = eye(dimT);

excited = [1; 0];
ground = [0; 1];
sp = excited * ground';
Sp = cell(1, Natoms);
Sm = cell(1, Natoms);
for iAtom = 1:Natoms
    Sp{iAtom} = local_atomic_operator(sp, iAtom, Natoms);
    Sm{iAtom} = Sp{iAtom}';
end

PsiPlus = (kron(ground, excited) + kron(excited, ground)) / sqrt(2);

xi = 1;
g = 0.5 * xi;
Omega = 0;
n = 8;
dx = 2;
x = [0 dx];
k0 = pi / 2;

fidelityCurves = zeros(numel(deltaFractions), numel(t));
initialFidelity = zeros(size(deltaFractions));
finalFidelity = zeros(size(deltaFractions));

for idx = 1:numel(deltaFractions)
    K = k0 * (1 + deltaFractions(idx));
    initialState = kron(excited, ground) - exp(-1i * K * dx) * kron(ground, excited);
    initialState = initialState / norm(initialState);

    term1 = 2 * exp(1i * K * abs(x - x.'));
    term2 = exp(1i * K * abs(x + n - x.'));
    term3 = exp(1i * K * abs(x - x.' - n));
    collectiveA = g^2 * (term1 + term2 + term3) / (2 * xi);

    J = imag(collectiveA);
    heff = atomic_markov_hamiltonian(Natoms, Omega, J, Sp, Sm);
    liouvillianHamiltonian = -1i * kron(Is, heff) + 1i * kron(heff.', Is);

    gamma = real(collectiveA);
    dissipator = zeros(dimT^2);
    for ii = 1:Natoms
        for jj = 1:Natoms
            dissipator = dissipator + mixed_lindblad_superoperator(gamma(ii, jj), ...
                Sp{ii}, Sm{jj}, dimT);
        end
    end

    liouvillian = dissipator + liouvillianHamiltonian;
    rho0 = initialState * initialState';
    fidelityCurves(idx, :) = pure_state_fidelity_curve(liouvillian, rho0, PsiPlus, t);
    initialFidelity(idx) = fidelityCurves(idx, 1);
    finalFidelity(idx) = fidelityCurves(idx, end);
end

pra = struct();
pra.description = "Markovian passive Bell-BIC benchmark following Effect_dk_Bell_states.m";
pra.deltaFractions = deltaFractions;
pra.timeXi = t;
pra.initialFidelity = initialFidelity;
pra.finalFidelity = finalFidelity;
pra.fidelityCurves = fidelityCurves;
pra.parameters = struct('gOverXi', g / xi, 'n', n, 'dx', dx, ...
    'k0', k0, 'tfXi', t(end));
end

function result = simulate_two_stage_detuned(params, K, u0, nu, Tprep, Tpassage, drivePhase, cdScale)
k0 = pi / 2;
omegaK = params.wc - 2 * params.xi * cos(K);

model = giant_atom_static_model(params);
model.wk = model.wk - omegaK;
model.Omega0 = 0;
model.detunedResonantK = K;
model.detunedOmega = omegaK;

controls = floquet_controls(u0, Tpassage, nu);
protocol = struct( ...
    'Tprep', Tprep, ...
    'Tpassage', Tpassage, ...
    'totalTime', Tprep + Tpassage, ...
    'drivePhase', drivePhase, ...
    'piPulseArea', pi);

tPrep = linspace(0, Tprep, 51);
tPassage = Tprep + linspace(0, Tpassage, max(151, ceil(60 * Tpassage)));
tEval = unique([tPrep tPassage]);

psi0 = zeros(model.dimension + 1, 1);
psi0(1) = 1;
rhsOptions.counterdiabaticScale = cdScale;
solverOptions = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, ...
    'MaxStep', min(Tprep/250, 0.03));

[tOut, psiOutRows] = ode113(@(t, psi) two_stage_rhs(t, psi, model, controls, protocol, rhsOptions), ...
    tEval, psi0, solverOptions);
psi = transpose(psiOutRows);
obs = two_stage_observables(psi);

[~, prepIdx] = min(abs(tOut - Tprep));

result.deltaKOverK0 = K / k0 - 1;
result.K = K;
result.detunedOmegaOverXi = omegaK / params.xi;
result.counterdiabaticScale = cdScale;
result.TpassageXi = params.xi * Tpassage;
result.afterPrepEgPopulation = obs.egPopulation(prepIdx);
result.finalAtomicSinglePopulation = obs.atomicSinglePopulation(end);
result.finalPhotonicPopulation = obs.photonicPopulation(end);
result.finalConcurrence = obs.concurrence(end);
result.finalConditionalConcurrence = obs.conditionalConcurrence(end);
result.finalBellPlusFidelity = obs.bellPlusFidelity(end);
result.finalConditionalBellPlusFidelity = obs.conditionalBellPlusFidelity(end);
result.maxNormError = max(abs(obs.norm - 1));
end

function paths = make_robustness_figure(summary, figureDir)
colors.pra = [0.60 0.60 0.60; 0.80 0.47 0.65; 0.90 0.62 0.00; 0.00 0.62 0.45; 0.82 0.33 0.00];
colors.without = [0.82 0.33 0.00];
colors.with = [0.00 0.45 0.70];

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18.5 17.5]);
layout = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

axA = nexttile(layout); hold(axA, 'on'); box(axA, 'on');
for idx = 1:numel(summary.pra.deltaFractions)
    plot(axA, summary.pra.timeXi, summary.pra.fidelityCurves(idx, :), ...
        'LineWidth', 1.6, 'Color', colors.pra(idx, :));
end
xlabel(axA, '$\xi t$', 'Interpreter', 'latex');
ylabel(axA, '$F_{\Psi^+}$', 'Interpreter', 'latex');
title(axA, 'PRA-style passive benchmark', 'Interpreter', 'latex');
legend(axA, compose('$\\delta K/K_0=%.3g$', summary.pra.deltaFractions), ...
    'Interpreter', 'latex', 'Location', 'southeast', 'Box', 'on', 'Color', 'white');
format_axis(axA);
xlim(axA, [0 500]);
ylim(axA, [-0.03 1.03]);
panel_label(axA, 'a');

axB = nexttile(layout); hold(axB, 'on'); box(axB, 'on');
deltaPercent = 100 * summary.deltaFractions;
fNoCD = [summary.withoutCD.finalBellPlusFidelity];
fCD = [summary.withCDReference.finalBellPlusFidelity];
plot(axB, deltaPercent, fNoCD, 'o-', 'Color', colors.without, ...
    'MarkerFaceColor', colors.without, 'LineWidth', 1.7);
plot(axB, deltaPercent, fCD, 's-', 'Color', colors.with, ...
    'MarkerFaceColor', colors.with, 'LineWidth', 1.9);
xlabel(axB, '$100\,\delta K/K_0$', 'Interpreter', 'latex');
ylabel(axB, '$F_{\Psi^+}(t_f)$', 'Interpreter', 'latex');
title(axB, 'Two-stage protocol, $\xi T=20$', 'Interpreter', 'latex');
legend(axB, {'without CD', 'with CD'}, 'Interpreter', 'latex', ...
    'Location', 'east', 'Box', 'on', 'Color', 'white');
format_axis(axB);
ylim(axB, [0.45 1.01]);
panel_label(axB, 'b');

axC = nexttile(layout); box(axC, 'on');
imagesc(axC, deltaPercent, summary.durationGridXi, summary.fidelityGrid);
set(axC, 'YDir', 'normal');
colormap(axC, parula);
clim(axC, [0.985 1.00]);
cb = colorbar(axC);
cb.Label.String = '';
xlabel(axC, '$100\,\delta K/K_0$', 'Interpreter', 'latex');
ylabel(axC, '$\xi T$', 'Interpreter', 'latex');
title(axC, 'CD fidelity versus time and mismatch', 'Interpreter', 'latex');
format_axis(axC);
panel_label(axC, 'c');

axD = nexttile(layout); hold(axD, 'on'); box(axD, 'on');
durationGrid = summary.durationGridXi;
selectedDeltas = [0 0.05 0.10];
lineStyles = {'-', '--', ':'};
for idx = 1:numel(selectedDeltas)
    [~, deltaIdx] = min(abs(summary.deltaFractions - selectedDeltas(idx)));
    plot(axD, durationGrid, summary.fidelityGrid(:, deltaIdx), lineStyles{idx}, ...
        'Color', colors.with, 'LineWidth', 1.9);
end
yline(axD, 0.99, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xlabel(axD, '$\xi T$', 'Interpreter', 'latex');
ylabel(axD, '$F_{\Psi^+}(t_f)$', 'Interpreter', 'latex');
title(axD, 'Time-to-fidelity tradeoff with CD', 'Interpreter', 'latex');
legend(axD, {'$0\%$', '$5\%$', '$10\%$'}, ...
    'Interpreter', 'latex', 'Location', 'southeast', 'Box', 'on', 'Color', 'white');
format_axis(axD);
ylim(axD, [0.90 1.01]);
panel_label(axD, 'd');

hide_all_axes_toolbars(fig);

paths.pdf = fullfile(figureDir, 'Figure_K_Robustness_Comparison.pdf');
paths.png = fullfile(figureDir, 'Figure_K_Robustness_Comparison.png');
paths.tiff = fullfile(figureDir, 'Figure_K_Robustness_Comparison.tif');
exportgraphics(fig, paths.pdf, 'ContentType', 'vector');
exportgraphics(fig, paths.png, 'Resolution', 600);
exportgraphics(fig, paths.tiff, 'Resolution', 600);
end

function jsonSummary = compact_json_summary(summary)
jsonSummary = summary;
jsonSummary.pra = compact_pra_summary(summary.pra);
end

function praSummary = compact_pra_summary(pra)
praSummary = rmfield(pra, {'timeXi', 'fidelityCurves'});
end

function sourceData = make_source_data(deltaFractions, durationGridXi, withoutCD, withCDReference, ...
    fidelityGrid, concurrenceGrid, atomicPopulationGrid)
numDelta = numel(deltaFractions);
sourceData = table();
sourceData.deltaK_over_K0 = deltaFractions(:);
sourceData.reference_T_xi = repmat(20, numDelta, 1);
sourceData.without_CD_final_F = [withoutCD.finalBellPlusFidelity].';
sourceData.without_CD_final_C = [withoutCD.finalConcurrence].';
sourceData.with_CD_final_F = [withCDReference.finalBellPlusFidelity].';
sourceData.with_CD_final_C = [withCDReference.finalConcurrence].';
sourceData.with_CD_final_atomic_population = [withCDReference.finalAtomicSinglePopulation].';

for idx = 1:numel(durationGridXi)
    tag = matlab.lang.makeValidName(sprintf('Txi_%g', durationGridXi(idx)));
    sourceData.(['CD_F_' tag]) = fidelityGrid(idx, :).';
    sourceData.(['CD_C_' tag]) = concurrenceGrid(idx, :).';
    sourceData.(['CD_Patom_' tag]) = atomicPopulationGrid(idx, :).';
end
end

function praCurveData = make_pra_curve_data(pra)
praCurveData = table();
praCurveData.xi_t = pra.timeXi(:);
for idx = 1:numel(pra.deltaFractions)
    tag = matlab.lang.makeValidName(sprintf('deltaK_over_K0_%g', pra.deltaFractions(idx)));
    praCurveData.(tag) = pra.fidelityCurves(idx, :).';
end
end

function result = empty_protocol_result()
result.deltaKOverK0 = [];
result.K = [];
result.detunedOmegaOverXi = [];
result.counterdiabaticScale = [];
result.TpassageXi = [];
result.afterPrepEgPopulation = [];
result.finalAtomicSinglePopulation = [];
result.finalPhotonicPopulation = [];
result.finalConcurrence = [];
result.finalConditionalConcurrence = [];
result.finalBellPlusFidelity = [];
result.finalConditionalBellPlusFidelity = [];
result.maxNormError = [];
end

function fidelityCurve = pure_state_fidelity_curve(L, rho0, target, t)
[V, D] = eig(L, 'vector');
coefficients = V \ rho0(:);
fidelityCurve = zeros(size(t));
for idx = 1:numel(t)
    rhoVec = V * (coefficients .* exp(D * t(idx)));
    rho = reshape(rhoVec, size(rho0));
    rho = (rho + rho') / 2;
    fidelityCurve(idx) = real(target' * rho * target);
end
fidelityCurve = max(0, min(1, fidelityCurve));
end

function Heff = atomic_markov_hamiltonian(Natoms, Omega, J, Sp, Sm)
Ha = 0;
Hi = 0;
for ii = 1:Natoms
    Ha = Ha + Omega * Sp{ii} * Sm{ii};
    for jj = 1:Natoms
        Hi = Hi + J(ii, jj) * (Sp{ii} * Sm{jj} + Sm{ii} * Sp{jj});
    end
end
Heff = Ha + Hi;
end

function L = mixed_lindblad_superoperator(gamma, A, B, dim)
It = eye(dim);
L = 2 * gamma * (kron(B, A.') - 0.5 * kron(A * B, It) - 0.5 * kron(It, (A * B).'));
end

function Sci = local_atomic_operator(sc, targetAtom, numAtoms)
Is = eye(2);
operators = cell(1, numAtoms);
for site = 1:numAtoms
    if site == targetAtom
        operators{site} = sc;
    else
        operators{site} = Is;
    end
end
Sci = operators{1};
for site = 2:numAtoms
    Sci = kron(Sci, operators{site});
end
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

function format_axis(ax)
grid(ax, 'on');
ax.GridAlpha = 0.12;
ax.LineWidth = 0.8;
ax.FontSize = 10;
ax.TickLabelInterpreter = 'latex';
end

function panel_label(ax, label)
text(ax, 0.02, 0.95, ['\textbf{' label '}'], ...
    'Units', 'normalized', 'Interpreter', 'latex', ...
    'FontSize', 11, 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', 'Margin', 1);
end

function hide_all_axes_toolbars(fig)
axesHandles = findall(fig, 'Type', 'axes');
for idx = 1:numel(axesHandles)
    ax = axesHandles(idx);
    if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
        ax.Toolbar.Visible = 'off';
    end
end
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
