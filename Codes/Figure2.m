%% Figure2.m
% Reproduce Fig. 2 of the manuscript.
%
% The figure summarizes the analytic connection between the geometry of the
% two giant atoms and the entanglement stored in the atomic BIC:
%   (left)  concurrence C(lambda) = 2 lambda/(1 + lambda^2);
%   (right) fidelity with a target Bell-like state as a function of the
%           polar coordinates x = lambda cos(theta), y = lambda sin(theta).
%
% Running this script creates a two-panel figure and exports PDF/PNG copies
% to ../outputs.

clear; close all; clc;

%% Output folder
scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    projectRoot = pwd;
else
    projectRoot = fileparts(fileparts(scriptPath));
end

outputDir = fullfile(projectRoot, 'outputs');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Analytic concurrence curve
lambdaCurve = logspace(-2, 2, 10000);
concurrenceCurve = 2 * abs(lambdaCurve) ./ (1 + lambdaCurve.^2);

%% Fidelity map in polar geometry coordinates
phiTarget = pi/4;                    % target phase in |Phi(phiTarget)>
lambdaMax = 5;
numLambda = 700;
numTheta = 700;

theta = linspace(0, 2*pi, numTheta);
lambda = linspace(0, lambdaMax, numLambda);
[THETA, LAMBDA] = meshgrid(theta, lambda);

concurrenceMap = 2 * abs(LAMBDA) ./ (1 + LAMBDA.^2);
fidelityMap = (1 - concurrenceMap .* cos(phiTarget - THETA)) / 2;

X = LAMBDA .* cos(THETA);
Y = LAMBDA .* sin(THETA);

isoFidelityLevels = [0.01 0.30 0.50 0.80 0.99];

%% Plot
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 25 10]);
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Left panel: concurrence.
nexttile;
semilogx(lambdaCurve, concurrenceCurve, 'LineWidth', 2.2, 'Color', [0.05 0.25 0.75]);
box on;
xlabel('$\lambda=n_1/n_2$', 'Interpreter', 'latex');
ylabel('$\mathcal{C}(\lambda)$', 'Interpreter', 'latex');
ylim([0 1.05]);
xlim([min(lambdaCurve) max(lambdaCurve)]);
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex');

% Right panel: fidelity landscape.
nexttile;
contourf(X, Y, fidelityMap, 40, 'LineColor', 'none');
hold on;
contour(X, Y, fidelityMap, isoFidelityLevels, '-k', 'LineWidth', 1.2);
hold off;
axis equal tight;
xlim([-lambdaMax lambdaMax]);
ylim([-lambdaMax lambdaMax]);
xticks(-lambdaMax:1:lambdaMax);
yticks(-lambdaMax:1:lambdaMax);
box on;
clim([0 1]);
cb = colorbar;
cb.Ticks = 0:0.2:1;
cb.TickLabelInterpreter = 'latex';
xlabel('$x=\lambda\cos\theta$', 'Interpreter', 'latex');
ylabel('$y=\lambda\sin\theta$', 'Interpreter', 'latex');
title('$\mathcal{F}$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex');

%% Export
exportgraphics(fig, fullfile(outputDir, 'Figure2.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outputDir, 'Figure2.png'), 'Resolution', 300);

fprintf('Figure 2 exported to:\n  %s\n  %s\n', ...
    fullfile(outputDir, 'Figure2.pdf'), ...
    fullfile(outputDir, 'Figure2.png'));
