function [B, dB, info] = compact_bic_branch(model, u0, theta)
%COMPACT_BIC_BRANCH Normalized dressed BIC and its theta derivative.
% Valid for connections (0,n),(2,n+2), integer n>2, at the band centre.
% The Fourier convention agrees with giant_atom_static_model: a photon at
% site x has momentum amplitudes exp(+i*k*x)/sqrt(Nc).

assert(model.x1 == 0 && model.Dx == 2 && model.n1 == model.n2);
assert(model.n1 > 2 && model.n1 == round(model.n1));
assert(abs(model.Omega0-model.wc) < 1e-12);
q = model.g*u0/model.xi;
p = [0; 0; (exp(1i*model.k) + exp(1i*model.k*(model.n1+1))) ...
    / sqrt(2*model.Nc)];
D = [cos(theta); sin(theta); zeros(model.Nc,1)];
dD = [-sin(theta); cos(theta); zeros(model.Nc,1)];
eta = atan(q*sin(2*theta)/sqrt(2));
etaPrime = sqrt(2)*q*cos(2*theta)/(1+q^2*sin(2*theta)^2/2);
B = cos(eta)*D + sin(eta)*p;
dB = cos(eta)*dD + etaPrime*(-sin(eta)*D + cos(eta)*p);
info = struct('D',D,'p',p,'eta',eta,'etaPrime',etaPrime, ...
    'atomicWeight',cos(eta)^2);
end
