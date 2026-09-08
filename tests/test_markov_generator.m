root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src','matlab'));
identity=eye(4); ground=[0;0;0;1]; rho=ground*ground';
for K=[0.7 pi/2 1.7 2.4]
    [L,H,Gamma]=giant_atom_markov_generator(0.5,1,K,8,2);
    assert(norm(identity(:)'*L)<1e-12);
    assert(norm(L*rho(:))<1e-12);
    assert(norm(H-H')<1e-12);
    assert(min(eig(Gamma))>-1e-12);
    bell=[0;1;1;0]/sqrt(2); rr=bell*bell';
    evolved=reshape(expm(3*L)*rr(:),4,4);
    assert(abs(trace(evolved)-1)<1e-12);
    assert(norm(evolved-evolved')<1e-12);
    assert(min(eig((evolved+evolved')/2))>-1e-12);
end
fprintf('Markov generator: trace, vacuum, Hermiticity, and positivity passed.\n');
