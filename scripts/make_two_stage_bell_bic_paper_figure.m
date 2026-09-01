function paths = make_two_stage_bell_bic_paper_figure()
%MAKE_TWO_STAGE_BELL_BIC_PAPER_FIGURE Export the paper-ready protocol figure.
%
% The script uses the production-grid source data produced by
% run_two_stage_bell_bic_protocol("paper"). It exports vector PDF, high
% resolution PNG/TIFF, and the source-data CSV.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);

dataPath = fullfile(projectRoot, 'data', 'paper', 'two_stage_bell_bic_protocol_paper.mat');
if ~exist(dataPath, 'file')
    error(['Missing paper data. Run addpath(''scripts''); ' ...
        'run_two_stage_bell_bic_protocol("paper") first.']);
end

loaded = load(dataPath, 'tOut', 'results', 'u1Global', 'u2Global', ...
    'omegaPi', 'omegaCD', 'targetConcurrence', 'summary', 'params', ...
    'Tprep');

outputDir = fullfile(projectRoot, 'figures');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

tXi = loaded.params.xi * loaded.tOut;
totalTimeXi = loaded.params.xi * loaded.summary.totalTime;
prepTimeXi = loaded.params.xi * loaded.Tprep;

colors.target = [0.15 0.15 0.15];
colors.without = [0.82 0.33 0.00];
colors.with = [0.00 0.45 0.70];
colors.green = [0.00 0.62 0.45];
colors.yellow = [0.90 0.62 0.00];
colors.purple = [0.58 0.40 0.74];

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 19 21]);
layout = tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

axA = nexttile(layout, [1 2]); hold(axA, 'on'); box(axA, 'on');
plot(axA, tXi, loaded.results(2).obs.groundPopulation, 'Color', [0.35 0.35 0.35], 'LineWidth', 1.4);
plot(axA, tXi, loaded.results(2).obs.egPopulation, 'Color', colors.without, 'LineWidth', 1.6);
plot(axA, tXi, loaded.results(2).obs.gePopulation, 'Color', colors.yellow, 'LineWidth', 1.6);
plot(axA, tXi, loaded.results(2).obs.photonicPopulation, 'Color', colors.purple, 'LineWidth', 1.4);
format_axis(axA, totalTimeXi, prepTimeXi);
ylabel(axA, 'population', 'Interpreter', 'latex');
ylim(axA, [-0.03 1.03]);
lgA = legend(axA, {'$P_{gg}$', '$P_{eg}$', '$P_{ge}$', '$P_{\gamma}$'}, ...
    'Interpreter', 'latex', 'Location', 'eastoutside', 'Box', 'off');
lgA.FontSize = 9;
panel_label(axA, 'a');
stage_labels(axA, prepTimeXi, totalTimeXi);

axB = nexttile(layout); hold(axB, 'on'); box(axB, 'on');
plot(axB, tXi, loaded.targetConcurrence, ':', 'Color', colors.target, 'LineWidth', 1.7);
plot(axB, tXi, loaded.results(1).obs.concurrence, 'Color', colors.without, 'LineWidth', 1.6);
plot(axB, tXi, loaded.results(2).obs.concurrence, 'Color', colors.with, 'LineWidth', 1.9);
format_axis(axB, totalTimeXi, prepTimeXi);
ylabel(axB, '$\mathcal{C}(t)$', 'Interpreter', 'latex');
ylim(axB, [-0.03 1.03]);
lgB = legend(axB, {'target', 'without CD', 'with CD'}, ...
    'Interpreter', 'latex', 'Location', 'southeast', 'Box', 'off');
lgB.FontSize = 9;
panel_label(axB, 'b');

axC = nexttile(layout); hold(axC, 'on'); box(axC, 'on');
plot(axC, tXi, loaded.results(1).obs.bellPlusFidelity, 'Color', colors.without, 'LineWidth', 1.6);
plot(axC, tXi, loaded.results(2).obs.bellPlusFidelity, 'Color', colors.with, 'LineWidth', 1.9);
plot(axC, tXi, loaded.results(2).obs.conditionalBellPlusFidelity, '--', ...
    'Color', colors.with, 'LineWidth', 1.5);
format_axis(axC, totalTimeXi, prepTimeXi);
ylabel(axC, '$F_{\Psi^+}(t)$', 'Interpreter', 'latex');
ylim(axC, [-0.03 1.03]);
lgC = legend(axC, {'without CD', 'with CD', 'with CD, cond.'}, ...
    'Interpreter', 'latex', 'Location', 'southeast', 'Box', 'off');
lgC.FontSize = 9;
panel_label(axC, 'c');

axD = nexttile(layout); hold(axD, 'on'); box(axD, 'on');
plot(axD, tXi, loaded.u1Global, 'Color', colors.with, 'LineWidth', 1.8);
plot(axD, tXi, loaded.u2Global, '--', 'Color', colors.without, 'LineWidth', 1.8);
format_axis(axD, totalTimeXi, prepTimeXi);
xlabel(axD, '$\xi t$', 'Interpreter', 'latex');
ylabel(axD, '$u_i(t)$', 'Interpreter', 'latex');
ylim(axD, [-0.03 1.0]);
lgD = legend(axD, {'$u_1$', '$u_2$'}, 'Interpreter', 'latex', ...
    'Location', 'northeast', 'Box', 'off');
lgD.FontSize = 9;
panel_label(axD, 'd');

axE = nexttile(layout); hold(axE, 'on'); box(axE, 'on');
yyaxis(axE, 'left');
plot(axE, tXi, loaded.omegaPi / loaded.params.xi, 'Color', colors.with, 'LineWidth', 1.8);
ylabel(axE, '$\Omega_{\pi}/\xi$', 'Interpreter', 'latex');
ylim(axE, [-0.08 3.5]);
axE.YAxis(1).Color = colors.with;
yyaxis(axE, 'right');
plot(axE, tXi, loaded.omegaCD / loaded.params.xi, 'k--', 'LineWidth', 1.8);
ylabel(axE, '$\Omega_{\rm CD}/\xi$', 'Interpreter', 'latex');
ylim(axE, [-0.003 0.085]);
axE.YAxis(2).Color = [0 0 0];
format_axis(axE, totalTimeXi, prepTimeXi);
xlabel(axE, '$\xi t$', 'Interpreter', 'latex');
lgE = legend(axE, {'$\Omega_{\pi}$', '$\Omega_{\rm CD}$'}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'Box', 'off');
lgE.FontSize = 9;
panel_label(axE, 'e');

hide_all_axes_toolbars(fig);

paths.pdf = fullfile(outputDir, 'Figure_TwoStage_Bell_BIC.pdf');
paths.png = fullfile(outputDir, 'Figure_TwoStage_Bell_BIC.png');
paths.tiff = fullfile(outputDir, 'Figure_TwoStage_Bell_BIC.tif');
exportgraphics(fig, paths.pdf, 'ContentType', 'vector');
exportgraphics(fig, paths.png, 'Resolution', 600);
exportgraphics(fig, paths.tiff, 'Resolution', 600);

fprintf('Paper-ready two-stage figure exported to:\n  %s\n  %s\n  %s\n', ...
    paths.pdf, paths.png, paths.tiff);
end

function format_axis(ax, totalTimeXi, prepTimeXi)
xlim(ax, [0 totalTimeXi]);
xline(ax, prepTimeXi, 'k--', 'LineWidth', 0.8, 'Alpha', 0.7, 'HandleVisibility', 'off');
grid(ax, 'on');
ax.GridAlpha = 0.12;
ax.LineWidth = 0.8;
ax.FontSize = 10;
ax.TickLabelInterpreter = 'latex';
ax.XTick = 0:5:totalTimeXi;
end

function panel_label(ax, label)
text(ax, 0.015, 0.95, ['\textbf{' label '}'], ...
    'Units', 'normalized', 'Interpreter', 'latex', ...
    'FontSize', 11, 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', 'Margin', 1);
end

function stage_labels(ax, prepTimeXi, totalTimeXi)
yl = ylim(ax);
text(ax, 0.50 * prepTimeXi, yl(1) + 0.18 * diff(yl), 'loading', ...
    'Interpreter', 'latex', 'FontSize', 9, 'HorizontalAlignment', 'center');
text(ax, prepTimeXi + 0.42 * (totalTimeXi - prepTimeXi), ...
    yl(2) - 0.12 * diff(yl), 'Floquet--CD Bell-BIC conversion', ...
    'Interpreter', 'latex', 'FontSize', 9, 'HorizontalAlignment', 'center');
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
