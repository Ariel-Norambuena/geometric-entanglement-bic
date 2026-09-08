% Exact eigenstate, derivative, restricted-CD residual, and degeneracy tests.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src','matlab'));
for n=[3 6 8]
    m=giant_atom_static_model(struct('Nc',204,'n1',n,'n2',n));
    c=floquet_controls(0.9,20,8);
    for t=[0 3 10 17 20]
        theta=c.theta(t);
        [b,db,info]=compact_bic_branch(m,c.u0,theta);
        assert(abs(norm(b)-1)<1e-12);
        assert(norm(floquet_effective_rhs(t,b,m,c))<1e-12);
        step=1e-6;
        finite=(compact_bic_branch(m,c.u0,theta+step) ...
            -compact_bic_branch(m,c.u0,theta-step))/(2*step);
        assert(norm(db-finite)<1e-9);
        cd=1i*[0 -1;1 0];
        h= [cd*b(1:2); zeros(m.Nc,1)];
        rr=1i*db-h; rr=rr-b*(b'*rr);
        assert(abs(norm(rr)-abs(info.etaPrime))<1e-12);
        ha=1i*info.etaPrime*(info.p*(info.D'*b)-info.D*(info.p'*b));
        assert(norm(1i*db-h-ha)<1e-12);
    end
end
m=giant_atom_static_model(struct('Nc',204));
theta=0.3; u0=0.9; u1=u0*sin(theta); u2=u0*cos(theta);
b1=[1;0;(m.g*u1/m.xi)/sqrt(m.Nc)*(exp(1i*m.k) ...
    -exp(3i*m.k)+exp(5i*m.k))];
b2=[0;1;(m.g*u2/m.xi)/sqrt(m.Nc)*(exp(3i*m.k) ...
    -exp(5i*m.k)+exp(7i*m.k))];
c.u=@(t) deal(u1,u2); c.thetaDot=@(t) 0;
assert(norm(floquet_effective_rhs(0,b1,m,c))<1e-12);
assert(norm(floquet_effective_rhs(0,b2,m,c))<1e-12);
assert(min(eig([b1 b2]'*[b1 b2]))>0.9);
b=cos(theta)*b1+sin(theta)*b2; b=b/norm(b);
assert(norm(b-compact_bic_branch(m,u0,theta))<1e-12);
fprintf('Compact BIC, derivative, and exact dressed-CD identities passed.\n');
