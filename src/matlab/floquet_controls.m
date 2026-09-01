function controls = floquet_controls(u0, T, nu)
%FLOQUET_CONTROLS Constant-norm dark-state-passage Floquet controls.
%
% The effective coupling path is
%   u1(t) = u0 sin(theta(t)),  u2(t) = u0 cos(theta(t)),
% with theta(t) moving from 0 to pi/4 through a smooth quintic schedule.
% The physical modulation depths beta_i(t) satisfy J0[beta_i(t)] = u_i(t)
% on the first monotonic Bessel branch.

arguments
    u0 (1,1) double {mustBePositive, mustBeLessThanOrEqual(u0, 1)} = 0.9
    T (1,1) double {mustBePositive} = 50
    nu (1,1) double {mustBePositive} = 8
end

controls.u0 = u0;
controls.T = T;
controls.nu = nu;
controls.theta = @(t) theta_schedule(t, T);
controls.thetaDot = @(t) theta_dot_schedule(t, T);
controls.u = @(t) coupling_controls(t, u0, T);
controls.uDot = @(t) coupling_control_derivatives(t, u0, T);
controls.beta = @(t) beta_controls(t, u0, T);
controls.betaDot = @(t) beta_control_derivatives(t, u0, T);
controls.atomicFrequencies = @(t, Omega0) driven_atomic_frequencies(t, Omega0, controls);
end

function theta = theta_schedule(t, T)
s = clipped_time(t, T);
h = 10*s.^3 - 15*s.^4 + 6*s.^5;
theta = (pi/4) * h;
end

function thetaDot = theta_dot_schedule(t, T)
s = clipped_time(t, T);
inside = (t >= 0) & (t <= T);
hDot = (30*s.^2 - 60*s.^3 + 30*s.^4) / T;
thetaDot = (pi/4) * hDot .* inside;
end

function [u1, u2] = coupling_controls(t, u0, T)
theta = theta_schedule(t, T);
u1 = u0 * sin(theta);
u2 = u0 * cos(theta);
end

function [u1Dot, u2Dot] = coupling_control_derivatives(t, u0, T)
theta = theta_schedule(t, T);
thetaDot = theta_dot_schedule(t, T);
u1Dot = u0 * cos(theta) .* thetaDot;
u2Dot = -u0 * sin(theta) .* thetaDot;
end

function [beta1, beta2] = beta_controls(t, u0, T)
[u1, u2] = coupling_controls(t, u0, T);
beta1 = besselj0_inverse_first_branch(u1);
beta2 = besselj0_inverse_first_branch(u2);
end

function [beta1Dot, beta2Dot] = beta_control_derivatives(t, u0, T)
[u1Dot, u2Dot] = coupling_control_derivatives(t, u0, T);
[beta1, beta2] = beta_controls(t, u0, T);
beta1Dot = beta_dot_from_u_dot(beta1, u1Dot);
beta2Dot = beta_dot_from_u_dot(beta2, u2Dot);
end

function betaDot = beta_dot_from_u_dot(beta, uDot)
denom = besselj(1, beta);
betaDot = zeros(size(beta));
regular = abs(denom) > 1e-12;
betaDot(regular) = -uDot(regular) ./ denom(regular);

if any(~regular & abs(uDot) > 1e-10, 'all')
    error('Encountered singular beta derivative away from a stationary endpoint.');
end
end

function omega = driven_atomic_frequencies(t, Omega0, controls)
[beta1, beta2] = controls.beta(t);
[beta1Dot, beta2Dot] = controls.betaDot(t);

omega = [
    Omega0 + controls.nu * beta1 .* cos(controls.nu * t) + beta1Dot .* sin(controls.nu * t);
    Omega0 + controls.nu * beta2 .* cos(controls.nu * t) + beta2Dot .* sin(controls.nu * t)
];
end

function s = clipped_time(t, T)
s = min(1, max(0, t ./ T));
end
