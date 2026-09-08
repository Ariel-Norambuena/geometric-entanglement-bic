function [B, info] = geometric_bic_state(model, u)
%GEOMETRIC_BIC_STATE Dressed state selected by the geometric amplitude rule.
% Band-centre family: n_i=4*l_i+2 and even separation between left connections.
% The amplitude rule is the original analytical construction of A. R. Legon.
% It cancels the first momentum derivative at both resonant roots, selecting
% a superposition within the BIC subspace. It does not remove all photons.

assert(abs(model.Omega0-model.wc)<1e-12);
assert(all(mod([model.n1 model.n2],4)==2));
assert(mod(model.Dx,2)==0 && mod(model.Nc,4)==0);
assert(numel(u)==2 && isreal(u) && u(2)~=0);
lambda=model.n1*u(1)/(model.n2*u(2));
a=[1;-lambda*exp(-1i*pi*model.Dx/2)]/sqrt(1+lambda^2);
n=[model.n1 model.n2]; x=[model.x1 model.x2];
phi=zeros(model.Nc,1);
for j=1:2
    sites=x(j)+(1:n(j)-1);
    weights=(model.g*u(j)/model.xi)*sin(pi*(1:n(j)-1)/2);
    phi=phi+a(j)*exp(1i*model.k*sites)*weights.'/sqrt(model.Nc);
end
Z=1/(1+norm(phi)^2);
B=sqrt(Z)*[a;phi];
info=struct('lambda',lambda,'a',a,'Z',Z,'photonNormBeforeNormalization',norm(phi)^2);
end
