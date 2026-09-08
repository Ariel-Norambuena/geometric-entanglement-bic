function dy = frozen_envelope_rhs(t, y, model, controls, cdScale)
%FROZEN_ENVELOPE_RHS Auxiliary ODE for u_j(t-tau) -> u_j(t) in the memory.
% y=[c1;c2;p1_k;p2_k], with p_jk=-i*g_jk^* int exp(-i*w_k*tau)c_j(t-tau)d tau.
% Only c1,c2 are physical amplitudes in this approximate equation. The p_jk
% are memory variables, not independent photon modes or normalized states.
% For constant u, phi_k=u1*p1_k+u2*p2_k recovers the exact vacuum-input model.
if nargin<5, cdScale=0; end
nc=model.Nc;
c=y(1:2); p1=y(3:2+nc); p2=y(3+nc:end);
[u1,u2]=controls.u(t);
phi=u1*p1+u2*p2;
omegaCD=cdScale*controls.thetaDot(t);
dc=-1i*[model.Omega0*c(1)+u1*transpose(model.g1k)*phi; ...
    model.Omega0*c(2)+u2*transpose(model.g2k)*phi];
dc=dc+omegaCD*[-c(2);c(1)];
dp1=-1i*(model.wk.*p1+conj(model.g1k)*c(1));
dp2=-1i*(model.wk.*p2+conj(model.g2k)*c(2));
dy=[dc;dp1;dp2];
end
