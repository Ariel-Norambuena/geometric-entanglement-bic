function summary = scan_protocol_parameter_robustness(runMode)
%SCAN_PROTOCOL_PARAMETER_ROBUSTNESS Map final Bell fidelity versus errors.
%
% This script probes the robustness of the complete two-stage protocol in
% the validated central-sideband effective model. Each panel reports the
% final unconditional Bell fidelity F_{Psi+}(t_f) after the loading pulse and
% the Floquet-CD passage.
%
% Scanned parameter pairs:
%   1. local loading pulse area and phase,
%   2. independent Floquet coupling calibration errors,
%   3. passage time and CD-amplitude calibration,
%   4. common and relative atomic detunings.

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
    params.Nc = 4 * 201;
    numGrid = 25;
    dataDir = fullfile(projectRoot, 'data', 'paper');
    figureDir = fullfile(projectRoot, 'figures');
    outputStem = 'protocol_parameter_robustness_paper';
    figureStem = 'Figure_Protocol_Parameter_Robustness';
else
    params.Nc = 4 * 51;
    numGrid = 11;
    dataDir = fullfile(projectRoot, 'data', 'development');
    figureDir = fullfile(projectRoot, 'outputs', 'floquet');
    outputStem = 'protocol_parameter_robustness_dev';
    figureStem = 'protocol_parameter_robustness_dev';
end

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
params.Omega0 = params.wc;

base = struct();
base.u0 = 0.9;
base.nu = 8 * params.xi;
base.Tprep = 2 / params.xi;
base.Tpassage = 20 / params.xi;
base.drivePhase = pi / 2;
base.piPulseArea = pi;
base.counterdiabaticScale = 1;
base.u1Scale = 1;
base.u2Scale = 1;
base.deltaCommon = 0;
base.deltaRelative = 0;

model = giant_atom_static_model(params);

scan = define_scans(numGrid);
numScans = numel(scan);
maps = repmat(empty_map(), numScans, 1);

for scanIdx = 1:numScans
    maps(scanIdx).name = scan(scanIdx).name;
    maps(scanIdx).x = scan(scanIdx).x;
    maps(scanIdx).y = scan(scanIdx).y;
    maps(scanIdx).xLabel = scan(scanIdx).xLabel;
    maps(scanIdx).yLabel = scan(scanIdx).yLabel;
    maps(scanIdx).xPlotScale = scan(scanIdx).xPlotScale;
    maps(scanIdx).yPlotScale = scan(scanIdx).yPlotScale;
    maps(scanIdx).title = scan(scanIdx).title;
    maps(scanIdx).fidelity = zeros(numel(scan(scanIdx).y), numel(scan(scanIdx).x));
    maps(scanIdx).concurrence = zeros(numel(scan(scanIdx).y), numel(scan(scanIdx).x));
    maps(scanIdx).atomicPopulation = zeros(numel(scan(scanIdx).y), numel(scan(scanIdx).x));
    maps(scanIdx).maxNormError = zeros(numel(scan(scanIdx).y), numel(scan(scanIdx).x));

    fprintf('Scanning %s (%d x %d)...\n', scan(scanIdx).name, ...
        numel(scan(scanIdx).x), numel(scan(scanIdx).y));

    for row = 1:numel(scan(scanIdx).y)
        for col = 1:numel(scan(scanIdx).x)
            protocol = apply_scan_point(base, scan(scanIdx).name, ...
                scan(scanIdx).x(col), scan(scanIdx).y(row));
            result = simulate_final_protocol(model, params, protocol);

            maps(scanIdx).fidelity(row, col) = result.finalBellPlusFidelity;
            maps(scanIdx).concurrence(row, col) = result.finalConcurrence;
            maps(scanIdx).atomicPopulation(row, col) = result.finalAtomicSinglePopulation;
            maps(scanIdx).maxNormError(row, col) = result.maxNormError;
        end

        fprintf('  row %2d/%2d complete\n', row, numel(scan(scanIdx).y));
    end

    maps(scanIdx).minimumFidelity = min(maps(scanIdx).fidelity, [], 'all');
    maps(scanIdx).maximumFidelity = max(maps(scanIdx).fidelity, [], 'all');
    maps(scanIdx).nominalFidelity = maps(scanIdx).fidelity( ...
        nearest_index(scan(scanIdx).y, 0), nearest_index(scan(scanIdx).x, 0));
end

summary = struct();
summary.runMode = runMode;
summary.model = "two-stage central-sideband effective finite-mode dynamics";
summary.description = "Final Bell fidelity maps under control and detuning errors";
summary.Nc = model.Nc;
summary.dimension = model.dimension + 1;
summary.kConvention = model.kConvention;
summary.xi = params.xi;
summary.gOverXi = params.g / params.xi;
summary.u0 = base.u0;
summary.nuOverXi = base.nu / params.xi;
summary.TprepXi = params.xi * base.Tprep;
summary.TpassageXi = params.xi * base.Tpassage;
summary.counterdiabaticScale = base.counterdiabaticScale;
summary.numGrid = numGrid;
summary.maps = maps;

paths = make_robustness_figure(summary, figureDir, figureStem);
summary.outputs = paths;

matPath = fullfile(dataDir, [outputStem '.mat']);
jsonPath = fullfile(dataDir, [outputStem '_summary.json']);
csvPath = fullfile(dataDir, [outputStem '_source_data.csv']);

sourceData = make_source_table(summary);
save(matPath, 'summary', 'params', 'base');
writetable(sourceData, csvPath);
summary.outputs.mat = string(matPath);
summary.outputs.json = string(jsonPath);
summary.outputs.csv = string(csvPath);
write_json(jsonPath, compact_json_summary(summary));

fprintf(['Protocol-parameter robustness complete.\n' ...
    '  global min F = %.6f\n' ...
    '  nominal F    = %.6f\n' ...
    'Outputs:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n'], ...
    min(arrayfun(@(m) m.minimumFidelity, maps)), ...
    maps(1).nominalFidelity, ...
    paths.pdf, paths.png, paths.tiff, matPath, jsonPath, csvPath);
end

function scan = define_scans(numGrid)
scan = repmat(struct('name', "", 'x', [], 'y', [], ...
    'xLabel', "", 'yLabel', "", 'title', "", ...
    'xPlotScale', 1, 'yPlotScale', 1), 4, 1);

scan(1).name = "loading";
scan(1).x = linspace(-0.10, 0.10, numGrid);
scan(1).y = linspace(-0.10, 0.10, numGrid);
scan(1).xLabel = '$100\,\delta A_{\pi}/A_{\pi}$';
scan(1).yLabel = '$\delta\varphi_{\pi}/\pi$';
scan(1).title = 'Loading pulse';
scan(1).xPlotScale = 100;

scan(2).name = "floquet_calibration";
scan(2).x = linspace(-0.10, 0.10, numGrid);
scan(2).y = linspace(-0.10, 0.10, numGrid);
scan(2).xLabel = '$100\,\delta u_1/u_1$';
scan(2).yLabel = '$100\,\delta u_2/u_2$';
scan(2).title = 'Floquet-coupling calibration';
scan(2).xPlotScale = 100;
scan(2).yPlotScale = 100;

scan(3).name = "timing_cd";
scan(3).x = linspace(-0.50, 0.50, numGrid);
scan(3).y = linspace(-0.30, 0.30, numGrid);
scan(3).xLabel = '$100\,\delta T/T_0$';
scan(3).yLabel = '$100\,\delta\alpha_{\rm CD}$';
scan(3).title = 'Timing and CD strength';
scan(3).xPlotScale = 100;
scan(3).yPlotScale = 100;

scan(4).name = "detuning";
scan(4).x = linspace(-0.05, 0.05, numGrid);
scan(4).y = linspace(-0.05, 0.05, numGrid);
scan(4).xLabel = '$\Delta_{\rm c}/\xi$';
scan(4).yLabel = '$\Delta_{\rm r}/\xi$';
scan(4).title = 'Atomic detunings';
end

function protocol = apply_scan_point(base, scanName, x, y)
protocol = base;

switch scanName
    case "loading"
        protocol.piPulseArea = base.piPulseArea * (1 + x);
        protocol.drivePhase = base.drivePhase + pi * y;

    case "floquet_calibration"
        protocol.u1Scale = 1 + x;
        protocol.u2Scale = 1 + y;

    case "timing_cd"
        protocol.Tpassage = base.Tpassage * (1 + x);
        protocol.counterdiabaticScale = base.counterdiabaticScale * (1 + y);

    case "detuning"
        protocol.deltaCommon = x;
        protocol.deltaRelative = y;

    otherwise
        error('Unknown scan name: %s.', scanName);
end

protocol.totalTime = protocol.Tprep + protocol.Tpassage;
end

function result = simulate_final_protocol(model, params, protocol)
controls = floquet_controls(protocol.u0, protocol.Tpassage, protocol.nu);

psi0 = zeros(model.dimension + 1, 1);
psi0(1) = 1;

solverOptions = odeset( ...
    'RelTol', 1e-8, ...
    'AbsTol', 1e-10, ...
    'MaxStep', min([protocol.Tprep / 200, protocol.Tpassage / 700, 0.04]));

tSpan = [0 protocol.Tprep protocol.totalTime];
rhsOptions.counterdiabaticScale = protocol.counterdiabaticScale;

[~, psiRows] = ode113(@(t, psi) two_stage_rhs(t, psi, model, controls, protocol, rhsOptions), ...
    tSpan, psi0, solverOptions);

psi = transpose(psiRows);
obs = two_stage_observables(psi);

result.finalAtomicSinglePopulation = obs.atomicSinglePopulation(end);
result.finalPhotonicPopulation = obs.photonicPopulation(end);
result.finalConcurrence = obs.concurrence(end);
result.finalConditionalConcurrence = obs.conditionalConcurrence(end);
result.finalBellPlusFidelity = obs.bellPlusFidelity(end);
result.finalConditionalBellPlusFidelity = obs.conditionalBellPlusFidelity(end);
result.maxNormError = max(abs(obs.norm - 1));
result.params = params;
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
[delta1, delta2] = atomic_detunings(protocol);

dc0 = -1i * (omegaPi / 2) * exp(-1i * protocol.drivePhase) * c1;
dc1 = -1i * ((model.Omega0 + delta1) * c1 + transpose(f1k) * phi ...
    + (omegaPi / 2) * exp(1i * protocol.drivePhase) * c0);
dc2 = -1i * ((model.Omega0 + delta2) * c2 + transpose(f2k) * phi);
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
    u2 = controls.u0 * protocol.u2Scale;
else
    [u1, u2] = controls.u(t - protocol.Tprep);
    u1 = protocol.u1Scale * u1;
    u2 = protocol.u2Scale * u2;
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

function [delta1, delta2] = atomic_detunings(protocol)
delta1 = protocol.deltaCommon + protocol.deltaRelative / 2;
delta2 = protocol.deltaCommon - protocol.deltaRelative / 2;
end

function paths = make_robustness_figure(summary, figureDir, figureStem)
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18.5 16.5]);
layout = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for idx = 1:numel(summary.maps)
    ax = nexttile(layout);
    draw_map(ax, summary.maps(idx), idx);
end

paths.pdf = string(fullfile(figureDir, [figureStem '.pdf']));
paths.png = string(fullfile(figureDir, [figureStem '.png']));
paths.tiff = string(fullfile(figureDir, [figureStem '.tif']));

exportgraphics(fig, paths.pdf, 'ContentType', 'vector');
exportgraphics(fig, paths.png, 'Resolution', 600);
exportgraphics(fig, paths.tiff, 'Resolution', 600);
end

function draw_map(ax, map, panelIndex)
xPlot = map.x * map.xPlotScale;
yPlot = map.y * map.yPlotScale;
imagesc(ax, xPlot, yPlot, map.fidelity);
set(ax, 'YDir', 'normal');
colormap(ax, turbo);
clim(ax, [0.90 1.00]);
hold(ax, 'on');

contour(ax, xPlot, yPlot, map.fidelity, [0.99 0.99], ...
    'LineColor', 'w', 'LineWidth', 1.0);
contour(ax, xPlot, yPlot, map.fidelity, [0.995 0.995], ...
    'LineColor', 'w', 'LineStyle', '--', 'LineWidth', 1.0);
plot(ax, 0, 0, 'wo', 'MarkerSize', 5, 'MarkerFaceColor', 'k', ...
    'LineWidth', 0.9);

cb = colorbar(ax);
cb.Label.String = '$F_{\Psi^+}(t_f)$';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';

xlabel(ax, map.xLabel, 'Interpreter', 'latex');
ylabel(ax, map.yLabel, 'Interpreter', 'latex');
title(ax, map.title, 'Interpreter', 'latex');
format_axis(ax);
panel_label(ax, char('a' + panelIndex - 1));
end

function format_axis(ax)
box(ax, 'on');
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 0.8;
if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end
end

function panel_label(ax, label)
text(ax, 0.03, 0.94, ['\textbf{' label '}'], 'Units', 'normalized', ...
    'Interpreter', 'latex', 'FontSize', 13, 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'Color', 'w');
end

function sourceData = make_source_table(summary)
sourceData = table();
for mapIdx = 1:numel(summary.maps)
    map = summary.maps(mapIdx);
    [xGrid, yGrid] = meshgrid(map.x, map.y);
    numPoints = numel(xGrid);
    mapTable = table();
    mapTable.map = repmat(string(map.name), numPoints, 1);
    mapTable.x = xGrid(:);
    mapTable.y = yGrid(:);
    mapTable.final_Bell_fidelity = map.fidelity(:);
    mapTable.final_concurrence = map.concurrence(:);
    mapTable.final_atomic_population = map.atomicPopulation(:);
    mapTable.max_norm_error = map.maxNormError(:);
    sourceData = [sourceData; mapTable]; %#ok<AGROW>
end
end

function jsonSummary = compact_json_summary(summary)
jsonSummary = summary;
for idx = 1:numel(jsonSummary.maps)
    jsonSummary.maps(idx).fidelity = [];
    jsonSummary.maps(idx).concurrence = [];
    jsonSummary.maps(idx).atomicPopulation = [];
    jsonSummary.maps(idx).maxNormError = [];
end
end

function map = empty_map()
map.name = "";
map.x = [];
map.y = [];
map.xLabel = "";
map.yLabel = "";
map.title = "";
map.xPlotScale = 1;
map.yPlotScale = 1;
map.fidelity = [];
map.concurrence = [];
map.atomicPopulation = [];
map.maxNormError = [];
map.minimumFidelity = [];
map.maximumFidelity = [];
map.nominalFidelity = [];
end

function idx = nearest_index(values, target)
[~, idx] = min(abs(values - target));
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
