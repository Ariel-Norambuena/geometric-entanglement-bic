function obs = observables_single_excitation(psi)
%OBSERVABLES_SINGLE_EXCITATION Basic observables for two atoms plus waveguide.
%
% Input psi may be a column vector or a matrix whose columns are states.

c1 = psi(1, :);
c2 = psi(2, :);
phi = psi(3:end, :);

atomicPopulation = abs(c1).^2 + abs(c2).^2;
photonicPopulation = sum(abs(phi).^2, 1);
normValue = atomicPopulation + photonicPopulation;

concurrence = 2 * abs(c1 .* c2);
conditionalConcurrence = nan(size(concurrence));
active = atomicPopulation > 0;
conditionalConcurrence(active) = concurrence(active) ./ atomicPopulation(active);

bellPlusFidelity = abs(c1 + c2).^2 / 2;
conditionalBellPlusFidelity = nan(size(bellPlusFidelity));
conditionalBellPlusFidelity(active) = bellPlusFidelity(active) ./ atomicPopulation(active);

obs.norm = normValue;
obs.atomicPopulation = atomicPopulation;
obs.photonicPopulation = photonicPopulation;
obs.concurrence = concurrence;
obs.conditionalConcurrence = conditionalConcurrence;
obs.bellPlusFidelity = bellPlusFidelity;
obs.conditionalBellPlusFidelity = conditionalBellPlusFidelity;
end
