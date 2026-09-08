% Verify the original amplitude rule using the complete finite Hamiltonian.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src','matlab'));
for n2=[6 10]
    m=giant_atom_static_model(struct('Nc',804,'n1',6,'n2',n2));
    for lambda=[0 0.25 1 2]
        u=0.9*[lambda*n2 m.n1]/hypot(lambda*n2,m.n1);
        [B,info]=geometric_bic_state(m,u);
        ctrl=floquet_controls(); ctrl.u=@(t) fixed_u(t,u);
        residual=norm(floquet_effective_rhs(0,B,m,ctrl)+1i*m.wc*B);
        assert(residual<1e-12 && abs(norm(B)-1)<1e-12);
        numerator=u(1)*conj(m.g1k)*info.a(1)+u(2)*conj(m.g2k)*info.a(2);
        denom=m.wc-m.wk; regular=abs(denom)>1e-12;
        phi=zeros(m.Nc,1); phi(regular)=numerator(regular)./denom(regular);
        reconstructed=sqrt(info.Z)*[info.a;phi];
        assert(norm(reconstructed-B)<1e-12);
        atomic=B(1:2)/norm(B(1:2));
        assert(abs(2*abs(prod(atomic))-2*lambda/(1+lambda^2))<1e-12);
    end
end

% Coincident identical connection sets admit an exactly photon-free state.
m=giant_atom_static_model(struct('Nc',804,'Dx',0,'g',0.4));
[B,info]=geometric_bic_state(m,[0.3 0.7]);
assert(norm(B(3:end))<1e-13 && abs(info.Z-1)<1e-13);

% Equal lengths and equal frozen couplings admit exact unitary Bell channels.
m=giant_atom_static_model(struct('Nc',804));
f=0.9/sqrt(2)*[transpose(m.g1k);transpose(m.g2k)];
U=[1 1;1 -1]/sqrt(2);
for z=[0.31+0.07i -0.42+0.03i]
    sigma=(f./transpose(z-m.wk))*f';
    transformed=U'*sigma*U;
    assert(norm(transformed-diag(diag(transformed)))<1e-12);
end

% At fixed controls, the memory-variable realization is exact, including
% its reconstructed photon amplitudes. No renormalization is applied.
m=giant_atom_static_model(struct('Nc',204));
ctrl=floquet_controls(); ctrl.u=@(t) fixed_u(t,[0.3 0.7]);
p0=zeros(m.dimension,1); p0(1)=1;
y0=zeros(2+2*m.Nc,1); y0(1)=1;
t=linspace(0,20,101); opts=odeset('RelTol',1e-11,'AbsTol',1e-13,'MaxStep',0.04);
[~,p]=ode113(@(t,p) floquet_effective_rhs(t,p,m,ctrl),t,p0,opts);
[~,y]=ode113(@(t,y) frozen_envelope_rhs(t,y,m,ctrl),t,y0,opts);
reconstructed=[y(:,1:2) 0.3*y(:,3:2+m.Nc)+0.7*y(:,3+m.Nc:end)];
assert(max(vecnorm(p-reconstructed,2,2))<1e-9);
fprintf('Geometric BIC identities and constant-envelope memory limit passed.\n');

function [u1,u2]=fixed_u(t,u)
u1=u(1)+0*t;u2=u(2)+0*t;
end
