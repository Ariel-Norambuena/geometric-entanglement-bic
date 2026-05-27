%% Figure3_Right_Panel.m
% Reproduce the decay-rate maps shown in the right part of Fig. 3.
%
% The maps evaluate the on-shell Markovian rates
%   Gamma_pm(Omega) = [2 g^2/(xi |sin k*|)] f_pm(k*)
% over the geometry plane x = lambda cos(theta), y = lambda sin(theta).
% The comparison between even and odd n2 highlights the parity dependence
% of the giant-atom form factor at k* = pi/2.
%
% Running this script exports PDF/PNG copies to ../outputs.

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

%% Physical and geometric parameters
xi = 1.0;
g = 0.1 * xi;
kstar = pi/2;

denom = xi * abs(sin(kstar));
if denom < 1e-12
    error('k* is too close to a band edge; the on-shell rate is singular.');
end

n2Even = 6;
n2Odd = 3;

lambdaMax = 3;
numLambda = 420;
numTheta = 720;

lambda = linspace(0, lambdaMax, numLambda);
theta = linspace(0, 2*pi, numTheta);
[LAMBDA, THETA] = meshgrid(lambda, theta);

X = LAMBDA .* cos(THETA);
Y = LAMBDA .* sin(THETA);

%% Compute rate maps
[gammaPlusEven, gammaMinusEven] = local_gamma_maps(LAMBDA, THETA, n2Even, kstar, g, denom, xi);
[gammaPlusOdd, gammaMinusOdd] = local_gamma_maps(LAMBDA, THETA, n2Odd, kstar, g, denom, xi);

allRates = [gammaPlusEven(:); gammaMinusEven(:); gammaPlusOdd(:); gammaMinusOdd(:)];
colorLimits = [min(allRates) max(allRates)];

%% Plot
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 24 18]);
tiledlayout(fig, 2, 2, 'Padding', 'loose', 'TileSpacing', 'compact');

plot_rate_panel(X, Y, gammaPlusEven, colorLimits, ...
    sprintf('$\\Gamma_{+}/\\xi$, $n_2=%d$ (even)', n2Even));
plot_rate_panel(X, Y, gammaPlusOdd, colorLimits, ...
    sprintf('$\\Gamma_{+}/\\xi$, $n_2=%d$ (odd)', n2Odd));
plot_rate_panel(X, Y, gammaMinusEven, colorLimits, ...
    sprintf('$\\Gamma_{-}/\\xi$, $n_2=%d$ (even)', n2Even));
plot_rate_panel(X, Y, gammaMinusOdd, colorLimits, ...
    sprintf('$\\Gamma_{-}/\\xi$, $n_2=%d$ (odd)', n2Odd));

colormap(turbo);

%% Export
exportgraphics(fig, fullfile(outputDir, 'Figure3_Right_Panel.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outputDir, 'Figure3_Right_Panel.png'), 'Resolution', 300);

fprintf('Figure 3 right panel exported to:\n  %s\n  %s\n', ...
    fullfile(outputDir, 'Figure3_Right_Panel.pdf'), ...
    fullfile(outputDir, 'Figure3_Right_Panel.png'));

%% Local functions
function [gammaPlus, gammaMinus] = local_gamma_maps(LAMBDA, THETA, n2, kstar, g, denom, xi)
% Evaluate normalized on-shell rates for a fixed n2.
%
% The radial coordinate fixes n1 through n1 = lambda n2. We keep n1
% continuous here to expose the full geometry landscape rather than only
% the experimentally selected integer separations.

    n1 = LAMBDA * n2;
    A1 = cos(kstar * n1 / 2);
    A2 = cos(kstar * n2 / 2);
    phase = THETA + (kstar / 2) .* (n2 - n1);

    fplus = abs(A1 + A2 .* exp(1i * phase)).^2;
    fminus = abs(A1 - A2 .* exp(1i * phase)).^2;

    gammaPlus = (2 * g^2 / denom) .* fplus / xi;
    gammaMinus = (2 * g^2 / denom) .* fminus / xi;
end

function plot_rate_panel(X, Y, rates, colorLimits, panelTitle)
% Draw one polar heat map using the same color scale as the other panels.

    nexttile;
    pcolor(X, Y, rates);
    shading flat;
    axis equal tight;
    box on;
    set(gca, 'CLim', colorLimits, 'FontSize', 13, 'TickLabelInterpreter', 'latex');
    title(panelTitle, 'Interpreter', 'latex');
    xlabel('$x=\lambda\cos\theta$', 'Interpreter', 'latex');
    ylabel('$y=\lambda\sin\theta$', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex');
end
