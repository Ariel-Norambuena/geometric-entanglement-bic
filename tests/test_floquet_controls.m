%% test_floquet_controls.m
% Unit tests for the Floquet dark-state-passage control functions.

clear; clc;

testDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(testDir);
addpath(fullfile(projectRoot, 'src', 'matlab'));

u0 = 0.9;
T = 50;
nu = 8;
controls = floquet_controls(u0, T, nu);

t = linspace(0, T, 101);
[u1, u2] = controls.u(t);
[beta1, beta2] = controls.beta(t);

assert(max(abs(besselj(0, beta1) - u1), [], 'all') < 1e-11);
assert(max(abs(besselj(0, beta2) - u2), [], 'all') < 1e-11);

assert(abs(u1(1)) < 1e-14);
assert(abs(u2(1) - u0) < 1e-14);
assert(abs(u1(end) - u0/sqrt(2)) < 1e-12);
assert(abs(u2(end) - u0/sqrt(2)) < 1e-12);

[beta1Dot, beta2Dot] = controls.betaDot(t);
assert(all(isfinite(beta1Dot), 'all'));
assert(all(isfinite(beta2Dot), 'all'));
assert(abs(beta1Dot(1)) < 1e-12);
assert(abs(beta2Dot(1)) < 1e-12);
assert(abs(beta1Dot(end)) < 1e-12);
assert(abs(beta2Dot(end)) < 1e-12);

% Two-quadrature identity: integrating Omega_i(t)-Omega0 must recover
% beta_i(t) sin(nu t) up to numerical quadrature error.
tFine = linspace(0, T, 20001);
omegaAtoms = controls.atomicFrequencies(tFine, 0);
[beta1Fine, beta2Fine] = controls.beta(tFine);
targetPhase = [
    beta1Fine .* sin(nu * tFine);
    beta2Fine .* sin(nu * tFine)
];
integratedPhase = cumtrapz(tFine, omegaAtoms, 2);
assert(max(abs(integratedPhase - targetPhase), [], 'all') < 1e-4);

fprintf('Floquet control tests passed.\n');
