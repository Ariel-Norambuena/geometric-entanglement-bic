function [L,H,Gamma,A]=giant_atom_markov_generator(g,xi,K,n,dx)
%GIANT_ATOM_MARKOV_GENERATOR Trace-preserving static tight-binding benchmark.
% Basis: |ee>, |eg>, |ge>, |gg>; rho is vectorized column by column.
% Gamma is the coefficient of 2*s_i*rho*s_j' - {s_j'*s_i,rho}.
% The exact on-shell tight-binding density of states supplies 1/sin(K).
assert(xi>0 && K>0 && K<pi);
x=[0 dx];
A=(g^2/(2*xi*sin(K)))*(2*exp(1i*K*abs(x-x.')) ...
    +exp(1i*K*abs(x+n-x.'))+exp(1i*K*abs(x-x.'-n)));
Gamma=real(A); J=imag(A);
sm=[0 0;1 0];
S={kron(sm,eye(2)),kron(eye(2),sm)};
H=zeros(4);
for i=1:2
    for j=1:2
        H=H+J(i,j)*(S{i}'*S{j});
    end
end
L=-1i*(kron(eye(4),H)-kron(H.',eye(4)));
for i=1:2
    for j=1:2
        loss=S{j}'*S{i};
        L=L+Gamma(i,j)*(2*kron(conj(S{j}),S{i}) ...
            -kron(eye(4),loss)-kron(loss.',eye(4)));
    end
end
end
