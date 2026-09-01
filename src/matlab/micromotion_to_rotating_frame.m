function psiRot = micromotion_to_rotating_frame(t, psiLab, controls)
%MICROMOTION_TO_ROTATING_FRAME Remove local atomic micromotion phases.
%
% The photonic amplitudes are unchanged. Atomic amplitudes transform as
% c_i^rot(t) = exp(+i beta_i(t) sin(nu t)) c_i^lab(t).

psiRot = psiLab;
[beta1, beta2] = controls.beta(t);
phase1 = exp(1i * beta1 .* sin(controls.nu * t));
phase2 = exp(1i * beta2 .* sin(controls.nu * t));

psiRot(1, :) = phase1 .* psiLab(1, :);
psiRot(2, :) = phase2 .* psiLab(2, :);
end
