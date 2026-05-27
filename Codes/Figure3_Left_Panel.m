%% Figure3_Left_Panel.m
% Reproduce the exact concurrence dynamics shown in the left part of Fig. 3.
%
% The script benchmarks the memory-kernel equation in the manuscript by
% propagating the full finite Hamiltonian in the single-excitation sector.
% The continuum is represented by k in [0, pi) plus explicit right/left
% propagation directions, which is equivalent to retaining the on-shell
% pair +/-k in the Brillouin zone.
%
% Solid curves: exact Hamiltonian propagation.
% Dashed curves: short-time Markov prediction C(t) = C(0) exp[-Gamma t].
%
% Running this script exports PDF/PNG copies to ../outputs. The exact
% diagonalizations are intentionally publication-grade and may take time.

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

%% Global parameters
xi = 1.0;
wc = 0.0 * xi;
g = 0.1 * xi;

fastTest = strcmpi(getenv('FIG3_FAST_TEST'), '1');
if fastTest
    fprintf('FIG3_FAST_TEST=1: using a reduced grid for a quick smoke test.\n');
    Nc = 4 * 51;
    outputSuffix = '_fast_test';
else
    Nc = 4 * 501;
    outputSuffix = '';
end

% A multiple of four includes k* = pi/2 exactly.
k = make_k_grid_0_pi(Nc);      % positive wave vectors; left movers are added explicitly
wk = wc - 2 * xi * cos(k);

kstar0 = pi/2;
Omega0 = wc - 2 * xi * cos(kstar0);

n2Ideal = 6;
n1Ideal = 6;
lambdaIdeal = n1Ideal / n2Ideal;

x1 = 0;
pct = 0.10;
n1LambdaDetuned = (1 + pct) * n2Ideal;
deltaTheta = pct * (2*pi);

kstarDetuned = kstar0 * (1 + pct);
OmegaKDetuned = wc - 2 * xi * cos(kstarDetuned);

%% Time grid
if fastTest
    tmax = 100 / xi;
    numTimes = 500;
else
    tmax = 1000 / xi;
    numTimes = 100000;
end
tvec = linspace(0, tmax, numTimes);

% Chunking avoids storing the full dim x numTimes wavefunction in memory.
timeChunkSize = min(2000, numTimes);

%% Panel (a): initial state aligned with |Psi+>
% The BIC atomic state is |eg> - exp(-i k* Dx)|ge>. Choosing k*Dx = pi
% turns this state into the symmetric Bell state |Psi+>.
phasePlusTarget = pi;
branchPlus = 0;
DxPlus = (phasePlusTarget + 2*pi*branchPlus) / kstar0;
phasePlus = exp(-1i * kstar0 * DxPlus);

psi0Plus = zeros(2, 1);
psi0Plus(1) = 1 / sqrt(1 + lambdaIdeal^2);
psi0Plus(2) = -lambdaIdeal * phasePlus / sqrt(1 + lambdaIdeal^2);

DxPlusThetaDetuned = DxPlus + deltaTheta / kstar0;

A_ideal = simulate_case_exact(k, wk, Omega0, g, x1, DxPlus, n1Ideal, n2Ideal, psi0Plus, tvec, timeChunkSize);
A_lam = simulate_case_exact(k, wk, Omega0, g, x1, DxPlus, n1LambdaDetuned, n2Ideal, psi0Plus, tvec, timeChunkSize);
A_theta = simulate_case_exact(k, wk, Omega0, g, x1, DxPlusThetaDetuned, n1Ideal, n2Ideal, psi0Plus, tvec, timeChunkSize);
A_kdet = simulate_case_exact(k, wk, OmegaKDetuned, g, x1, DxPlus, n1Ideal, n2Ideal, psi0Plus, tvec, timeChunkSize);

[gammaPlusIdeal, ~, kOnShellIdeal] = gamma_analytic(wc, xi, g, Omega0, DxPlus, n1Ideal, n2Ideal);
[gammaPlusLambda, ~, ~] = gamma_analytic(wc, xi, g, Omega0, DxPlus, n1LambdaDetuned, n2Ideal);
[gammaPlusTheta, ~, ~] = gamma_analytic(wc, xi, g, Omega0, DxPlusThetaDetuned, n1Ideal, n2Ideal);
[gammaPlusKdet, ~, kOnShellDetuned] = gamma_analytic(wc, xi, g, OmegaKDetuned, DxPlus, n1Ideal, n2Ideal);

A_ideal.Cexp = exp_prediction(A_ideal.t, A_ideal.C(1), gammaPlusIdeal);
A_lam.Cexp = exp_prediction(A_lam.t, A_lam.C(1), gammaPlusLambda);
A_theta.Cexp = exp_prediction(A_theta.t, A_theta.C(1), gammaPlusTheta);
A_kdet.Cexp = exp_prediction(A_kdet.t, A_kdet.C(1), gammaPlusKdet);

%% Panel (b): initial state aligned with |Psi->
% Choosing k*Dx = 0 turns the same BIC construction into the antisymmetric
% Bell state |Psi->.
phaseMinusTarget = 0;
branchMinus = 0;
DxMinus = (phaseMinusTarget + 2*pi*branchMinus) / kstar0;
phaseMinus = exp(-1i * kstar0 * DxMinus);

psi0Minus = zeros(2, 1);
psi0Minus(1) = 1 / sqrt(1 + lambdaIdeal^2);
psi0Minus(2) = -lambdaIdeal * phaseMinus / sqrt(1 + lambdaIdeal^2);

DxMinusThetaDetuned = DxMinus + deltaTheta / kstar0;

B_ideal = simulate_case_exact(k, wk, Omega0, g, x1, DxMinus, n1Ideal, n2Ideal, psi0Minus, tvec, timeChunkSize);
B_lam = simulate_case_exact(k, wk, Omega0, g, x1, DxMinus, n1LambdaDetuned, n2Ideal, psi0Minus, tvec, timeChunkSize);
B_theta = simulate_case_exact(k, wk, Omega0, g, x1, DxMinusThetaDetuned, n1Ideal, n2Ideal, psi0Minus, tvec, timeChunkSize);
B_kdet = simulate_case_exact(k, wk, OmegaKDetuned, g, x1, DxMinus, n1Ideal, n2Ideal, psi0Minus, tvec, timeChunkSize);

[~, gammaMinusIdeal, ~] = gamma_analytic(wc, xi, g, Omega0, DxMinus, n1Ideal, n2Ideal);
[~, gammaMinusLambda, ~] = gamma_analytic(wc, xi, g, Omega0, DxMinus, n1LambdaDetuned, n2Ideal);
[~, gammaMinusTheta, ~] = gamma_analytic(wc, xi, g, Omega0, DxMinusThetaDetuned, n1Ideal, n2Ideal);
[~, gammaMinusKdet, ~] = gamma_analytic(wc, xi, g, OmegaKDetuned, DxMinus, n1Ideal, n2Ideal);

B_ideal.Cexp = exp_prediction(B_ideal.t, B_ideal.C(1), gammaMinusIdeal);
B_lam.Cexp = exp_prediction(B_lam.t, B_lam.C(1), gammaMinusLambda);
B_theta.Cexp = exp_prediction(B_theta.t, B_theta.C(1), gammaMinusTheta);
B_kdet.Cexp = exp_prediction(B_kdet.t, B_kdet.C(1), gammaMinusKdet);

%% Plot exact dynamics and Markov overlays
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 24 11]);
tiledlayout(fig, 1, 2, 'Padding', 'loose', 'TileSpacing', 'compact');

nexttile; hold on; box on; grid on;
h1 = plot(xi * A_ideal.t, A_ideal.C, 'LineWidth', 2.2);
h2 = plot(xi * A_lam.t, A_lam.C, 'LineWidth', 1.6);
h3 = plot(xi * A_theta.t, A_theta.C, 'LineWidth', 1.6);
h4 = plot(xi * A_kdet.t, A_kdet.C, 'LineWidth', 1.6);
plot(xi * A_ideal.t, A_ideal.Cexp, '--', 'LineWidth', 1.6, 'Color', h1.Color);
plot(xi * A_lam.t, A_lam.Cexp, '--', 'LineWidth', 1.6, 'Color', h2.Color);
plot(xi * A_theta.t, A_theta.Cexp, '--', 'LineWidth', 1.6, 'Color', h3.Color);
plot(xi * A_kdet.t, A_kdet.Cexp, '--', 'LineWidth', 1.6, 'Color', h4.Color);
xlabel('$\xi t$', 'Interpreter', 'latex');
ylabel('$\mathcal{C}(t)$', 'Interpreter', 'latex');
title('$|\Psi(0)\rangle = |\Psi^{+}\rangle$', 'Interpreter', 'latex', 'FontSize', 11);
ylim([-0.05 1.05]);
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
hide_axes_toolbar(gca);

nexttile; hold on; box on; grid on;
j1 = plot(xi * B_ideal.t, B_ideal.C, 'LineWidth', 2.2);
j2 = plot(xi * B_lam.t, B_lam.C, 'LineWidth', 1.6);
j3 = plot(xi * B_theta.t, B_theta.C, 'LineWidth', 1.6);
j4 = plot(xi * B_kdet.t, B_kdet.C, 'LineWidth', 1.6);
plot(xi * B_ideal.t, B_ideal.Cexp, '--', 'LineWidth', 1.6, 'Color', j1.Color);
plot(xi * B_lam.t, B_lam.Cexp, '--', 'LineWidth', 1.6, 'Color', j2.Color);
plot(xi * B_theta.t, B_theta.Cexp, '--', 'LineWidth', 1.6, 'Color', j3.Color);
plot(xi * B_kdet.t, B_kdet.Cexp, '--', 'LineWidth', 1.6, 'Color', j4.Color);
xlabel('$\xi t$', 'Interpreter', 'latex');
ylabel('$\mathcal{C}(t)$', 'Interpreter', 'latex');
title('$|\Psi(0)\rangle = |\Psi^{-}\rangle$', 'Interpreter', 'latex', 'FontSize', 11);
ylim([-0.05 1.05]);
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
hide_axes_toolbar(gca);

lg = legend({ ...
    'ideal (num)', ...
    sprintf('$\\lambda=1(1+%.0f\\%%)$ (num)', 100*pct), ...
    '$\delta\theta=0.1\times 2\pi$ (num)', ...
    sprintf('$k^{\\star}=k^{\\star}_0(1+%.0f\\%%)$ (num)', 100*pct), ...
    'ideal (exp)', ...
    '$\lambda$ detune (exp)', ...
    '$\theta$ detune (exp)', ...
    '$k^{\star}$ detune (exp)'}, ...
    'Interpreter', 'latex');
lg.Layout.Tile = 'south';
lg.NumColumns = 4;
lg.FontSize = 12;

%% Export
pdfPath = fullfile(outputDir, ['Figure3_Left_Panel' outputSuffix '.pdf']);
pngPath = fullfile(outputDir, ['Figure3_Left_Panel' outputSuffix '.png']);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');
exportgraphics(fig, pngPath, 'Resolution', 300);

fprintf('\nAnalytic on-shell k* (ideal Omega):     %.6f\n', kOnShellIdeal);
fprintf('Analytic on-shell k* (detuned Omega):   %.6f\n', kOnShellDetuned);

fprintf('\nPanel (a), plus channel rates\n');
fprintf('  ideal : Gamma+ = %.6e\n', gammaPlusIdeal);
fprintf('  lambda: Gamma+ = %.6e\n', gammaPlusLambda);
fprintf('  theta : Gamma+ = %.6e\n', gammaPlusTheta);
fprintf('  kdet  : Gamma+ = %.6e\n', gammaPlusKdet);

fprintf('\nPanel (b), minus channel rates\n');
fprintf('  ideal : Gamma- = %.6e\n', gammaMinusIdeal);
fprintf('  lambda: Gamma- = %.6e\n', gammaMinusLambda);
fprintf('  theta : Gamma- = %.6e\n', gammaMinusTheta);
fprintf('  kdet  : Gamma- = %.6e\n', gammaMinusKdet);

fprintf('\nFigure 3 left panel exported to:\n  %s\n  %s\n', ...
    pdfPath, ...
    pngPath);

%% Local functions
function k = make_k_grid_0_pi(Nc)
% Uniform positive-momentum grid in [0, pi).

    m = (0:Nc-1).';
    k = pi * m / Nc;
end

function out = simulate_case_exact(k, wk, Omega, g, x1, Dx, n1, n2, psi0Atomic, tvec, chunkSize)
% Propagate the finite Hamiltonian exactly in the one-excitation manifold.
%
% Basis ordering:
%   1              |e,g,0>
%   2              |g,e,0>
%   3:2+Nc         |g,g,1_{k,R}>
%   3+Nc:2+2*Nc    |g,g,1_{k,L}>

    Nc = numel(k);
    dim = 2 + 2*Nc;
    x2 = x1 + Dx;

    A1 = cos(k * n1 / 2);
    A2 = cos(k * n2 / 2);

    phase1 = x1 + n1/2;
    phase2 = x2 + n2/2;

    % This normalization matches a full Brillouin-zone discretization with
    % 2*Nc waveguide modes while keeping right and left movers explicit.
    prefactor = 2 * g / sqrt(2*Nc);

    g1R = prefactor .* A1 .* exp(-1i * k * phase1);
    g2R = prefactor .* A2 .* exp(-1i * k * phase2);
    g1L = prefactor .* A1 .* exp(+1i * k * phase1);
    g2L = prefactor .* A2 .* exp(+1i * k * phase2);

    H = zeros(dim, dim);
    H(1, 1) = Omega;
    H(2, 2) = Omega;
    H(3:2+Nc, 3:2+Nc) = diag(wk);
    H(3+Nc:end, 3+Nc:end) = diag(wk);

    H(1, 3:2+Nc) = g1R.';
    H(2, 3:2+Nc) = g2R.';
    H(3:2+Nc, 1) = conj(g1R);
    H(3:2+Nc, 2) = conj(g2R);

    H(1, 3+Nc:end) = g1L.';
    H(2, 3+Nc:end) = g2L.';
    H(3+Nc:end, 1) = conj(g1L);
    H(3+Nc:end, 2) = conj(g2L);

    psi0 = zeros(dim, 1);
    psi0(1:2) = psi0Atomic(1:2);

    [V, D] = eig(H, 'vector');
    amplitudes = V' * psi0;
    atomicRows = V(1:2, :);

    numTimes = numel(tvec);
    c1 = zeros(numTimes, 1);
    c2 = zeros(numTimes, 1);
    chunkSize = max(1, min(chunkSize, numTimes));

    for first = 1:chunkSize:numTimes
        last = min(first + chunkSize - 1, numTimes);
        idx = first:last;
        phase = exp(-1i * D(:) * tvec(idx));
        atomicState = atomicRows * (amplitudes .* phase);
        c1(idx) = atomicState(1, :).';
        c2(idx) = atomicState(2, :).';
    end

    concurrence = 2 * abs(c1 .* c2);

    % The waveguide population is reconstructed from unitarity for a compact
    % diagnostic without storing all photonic amplitudes at every time.
    photonPopulation = max(0, 1 - abs(c1).^2 - abs(c2).^2);
    normTotal = abs(c1).^2 + abs(c2).^2 + photonPopulation;

    out.t = tvec(:);
    out.C = concurrence;
    out.norm = normTotal;
end

function [gammaPlus, gammaMinus, kstar] = gamma_analytic(wc, xi, g, Omega, Dx, n1, n2)
% On-shell Markov decay rates for the symmetric and antisymmetric channels.

    arg = (wc - Omega) / (2*xi);
    arg = min(1, max(-1, arg));
    kstar = acos(arg);

    denom = xi * abs(sin(kstar));
    if denom < 1e-12
        gammaPlus = NaN;
        gammaMinus = NaN;
        return;
    end

    A1 = cos(kstar * n1 / 2);
    A2 = cos(kstar * n2 / 2);
    phase = kstar * (Dx + (n2 - n1)/2);

    fplus = abs(A1 + A2 * exp(1i * phase))^2;
    fminus = abs(A1 - A2 * exp(1i * phase))^2;

    gammaPlus = (2 * g^2 / denom) * fplus;
    gammaMinus = (2 * g^2 / denom) * fminus;
end

function Cexp = exp_prediction(t, C0, gamma)
% Build the short-time exponential prediction robustly.

    if ~isfinite(gamma) || gamma < 0
        Cexp = NaN(size(t));
        return;
    end
    Cexp = C0 * exp(-gamma * t);
end

function hide_axes_toolbar(ax)
% Prevent MATLAB's interactive axes toolbar from appearing in exported files.

    if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
        ax.Toolbar.Visible = 'off';
    end
end
