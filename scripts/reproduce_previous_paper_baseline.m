function baseline = reproduce_previous_paper_baseline(mode)
%REPRODUCE_PREVIOUS_PAPER_BASELINE Freeze selected observables from old Fig. 3.
%
% Usage:
%   reproduce_previous_paper_baseline()
%   reproduce_previous_paper_baseline("write")
%   reproduce_previous_paper_baseline("verify")
%
% The script records a compact machine-readable baseline for the exact
% finite-mode dynamics used in the previous giant-atom BIC paper. It uses
% the same [0, pi) + explicit R/L convention as Codes/Figure3_Left_Panel.m.
% This convention is equivalent to a full Brillouin-zone discretization but
% makes the two on-shell branches explicit in the basis.

if nargin < 1
    mode = "write";
else
    mode = string(mode);
end

if ~ismember(mode, ["write", "verify"])
    error('Unknown mode "%s". Use "write" or "verify".', mode);
end

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
dataDir = fullfile(projectRoot, 'data');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

baselinePath = fullfile(dataDir, 'baseline_previous_paper.json');

baseline = compute_baseline(projectRoot);

if mode == "write"
    write_json(baselinePath, baseline);
    fprintf('Baseline written to %s\n', baselinePath);
else
    if ~exist(baselinePath, 'file')
        error('Baseline file does not exist: %s', baselinePath);
    end

    reference = jsondecode(fileread(baselinePath));
    verify_against_reference(baseline, reference);
    fprintf('Baseline regression passed against %s\n', baselinePath);
end
end

function baseline = compute_baseline(projectRoot)
timer = tic;

params.xi = 1.0;
params.wc = 0.0 * params.xi;
params.g = 0.1 * params.xi;
params.NcPositive = 4 * 501;
params.totalDimension = 2 + 2 * params.NcPositive;
params.kConvention = "[0, pi) with explicit right/left propagation directions";
params.kstar0 = pi/2;
params.Omega0 = params.wc - 2 * params.xi * cos(params.kstar0);
params.n1Ideal = 6;
params.n2Ideal = 6;
params.lambdaIdeal = params.n1Ideal / params.n2Ideal;
params.x1 = 0;
params.detuningFraction = 0.10;
params.n1LambdaDetuned = (1 + params.detuningFraction) * params.n2Ideal;
params.deltaTheta = params.detuningFraction * (2*pi);
params.kstarDetuned = params.kstar0 * (1 + params.detuningFraction);
params.OmegaKDetuned = params.wc - 2 * params.xi * cos(params.kstarDetuned);
params.sampleTimesXi = [0 100 200 400 600 800 1000];
params.referenceTolerance = 5e-8;

k = make_k_grid_0_pi(params.NcPositive);
wk = params.wc - 2 * params.xi * cos(k);
tvec = params.sampleTimesXi / params.xi;

DxPlus = pi / params.kstar0;
phasePlus = exp(-1i * params.kstar0 * DxPlus);
psi0Plus = [
    1 / sqrt(1 + params.lambdaIdeal^2);
    -params.lambdaIdeal * phasePlus / sqrt(1 + params.lambdaIdeal^2)
];

DxMinus = 0 / params.kstar0;
phaseMinus = exp(-1i * params.kstar0 * DxMinus);
psi0Minus = [
    1 / sqrt(1 + params.lambdaIdeal^2);
    -params.lambdaIdeal * phaseMinus / sqrt(1 + params.lambdaIdeal^2)
];

cases = [
    make_case("plus_ideal", params.Omega0, DxPlus, params.n1Ideal, params.n2Ideal, psi0Plus)
    make_case("plus_lambda_detuned", params.Omega0, DxPlus, params.n1LambdaDetuned, params.n2Ideal, psi0Plus)
    make_case("minus_ideal", params.Omega0, DxMinus, params.n1Ideal, params.n2Ideal, psi0Minus)
    make_case("minus_lambda_detuned", params.Omega0, DxMinus, params.n1LambdaDetuned, params.n2Ideal, psi0Minus)
];

for idx = 1:numel(cases)
    caseTimer = tic;
    out = simulate_case_exact(k, wk, cases(idx).Omega, params.g, params.x1, ...
        cases(idx).Dx, cases(idx).n1, cases(idx).n2, cases(idx).psi0Atomic, tvec);
    cases(idx).concurrence = out.C(:).';
    cases(idx).norm = out.norm(:).';
    cases(idx).maxNormError = max(abs(out.norm - 1));
    cases(idx).runtimeSeconds = toc(caseTimer);
    cases(idx).psi0Atomic = encode_complex_vector(cases(idx).psi0Atomic);
end

[~, gitCommit] = system(sprintf('git -C "%s" rev-parse HEAD', projectRoot));

baseline.schemaVersion = "1.0";
baseline.generatedAt = char(datetime('now', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
baseline.gitCommit = strtrim(gitCommit);
baseline.matlabVersion = version;
baseline.computer = computer;
baseline.productionGrid = true;
baseline.parameters = params;
baseline.cases = cases;
baseline.totalRuntimeSeconds = toc(timer);
end

function item = make_case(name, Omega, Dx, n1, n2, psi0Atomic)
item.name = name;
item.Omega = Omega;
item.Dx = Dx;
item.n1 = n1;
item.n2 = n2;
item.psi0Atomic = psi0Atomic;
item.concurrence = [];
item.norm = [];
item.maxNormError = [];
item.runtimeSeconds = [];
end

function k = make_k_grid_0_pi(Nc)
m = (0:Nc-1).';
k = pi * m / Nc;
end

function out = simulate_case_exact(k, wk, Omega, g, x1, Dx, n1, n2, psi0Atomic, tvec)
Nc = numel(k);
dim = 2 + 2*Nc;
x2 = x1 + Dx;

A1 = cos(k * n1 / 2);
A2 = cos(k * n2 / 2);

phase1 = x1 + n1/2;
phase2 = x2 + n2/2;
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
amplitudes = V \ psi0;
psiT = V * (amplitudes .* exp(-1i * D(:) * tvec));

c1 = psiT(1, :).';
c2 = psiT(2, :).';

out.C = 2 * abs(c1 .* c2);
out.norm = sum(abs(psiT).^2, 1).';
end

function encoded = encode_complex_vector(v)
encoded = struct('real', real(v(:)).', 'imag', imag(v(:)).');
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

function verify_against_reference(current, reference)
tolerance = reference.parameters.referenceTolerance;
if current.parameters.NcPositive ~= reference.parameters.NcPositive
    error('NcPositive changed: current %d, reference %d.', ...
        current.parameters.NcPositive, reference.parameters.NcPositive);
end

if current.parameters.totalDimension ~= reference.parameters.totalDimension
    error('Total dimension changed: current %d, reference %d.', ...
        current.parameters.totalDimension, reference.parameters.totalDimension);
end

for idx = 1:numel(reference.cases)
    currentCase = current.cases(idx);
    referenceCase = reference.cases(idx);
    if string(currentCase.name) ~= string(referenceCase.name)
        error('Case order changed at index %d: current %s, reference %s.', ...
            idx, currentCase.name, referenceCase.name);
    end

    maxConcurrenceError = max(abs(currentCase.concurrence(:) - referenceCase.concurrence(:)));
    if maxConcurrenceError > tolerance
        error('Case %s changed by %.3e, above tolerance %.3e.', ...
            char(currentCase.name), maxConcurrenceError, tolerance);
    end

    if currentCase.maxNormError > 1e-8
        error('Case %s has norm error %.3e, above 1e-8.', ...
            char(currentCase.name), currentCase.maxNormError);
    end
end
end
