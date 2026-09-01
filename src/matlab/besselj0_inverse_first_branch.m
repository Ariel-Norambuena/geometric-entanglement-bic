function beta = besselj0_inverse_first_branch(u)
%BESSELJ0_INVERSE_FIRST_BRANCH Invert J0(beta)=u on 0 <= beta <= z01.
%
% The first branch is monotonic from J0(0)=1 to J0(z01)=0. A bracketed
% solver is used for every nontrivial value to avoid Newton jumps to another
% branch.

z01 = 2.404825557695773;
tol = 1e-13;

if any(u(:) < -tol | u(:) > 1 + tol)
    error('u must lie in the first-branch range 0 <= u <= 1.');
end

u = min(1, max(0, u));
beta = zeros(size(u));

for idx = 1:numel(u)
    target = u(idx);
    if abs(target - 1) < tol
        beta(idx) = 0;
    elseif abs(target) < tol
        beta(idx) = z01;
    else
        beta(idx) = fzero(@(x) besselj(0, x) - target, [0, z01]);
    end
end
end
