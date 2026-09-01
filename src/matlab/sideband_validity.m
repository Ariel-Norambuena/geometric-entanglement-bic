function report = sideband_validity(model, controls, tSamples, mMax)
%SIDEBAND_VALIDITY Quantify separation of closed Floquet sidebands.
%
% The central-sideband model is safest when the off-resonant sideband gap
% Delta_F is large compared with both the sideband coupling strengths and
% the envelope-rate scale max|betaDot|.

arguments
    model struct
    controls struct
    tSamples double
    mMax (1,1) double {mustBeInteger, mustBePositive} = 6
end

mValues = [-mMax:-1 1:mMax];
gapValues = inf(numel(mValues), 1);
for idx = 1:numel(mValues)
    m = mValues(idx);
    gapValues(idx) = min(abs(model.Omega0 - model.wk + m * controls.nu));
end
deltaF = min(gapValues);

maxBareCoupling = max([abs(model.g1k); abs(model.g2k)]);
maxSidebandCoupling = 0;
maxBetaDot = 0;

for t = tSamples(:).'
    [beta1, beta2] = controls.beta(t);
    [beta1Dot, beta2Dot] = controls.betaDot(t);
    maxBetaDot = max(maxBetaDot, max(abs([beta1Dot beta2Dot])));

    for m = mValues
        harmonicWeight = max(abs([besselj(m, beta1), besselj(m, beta2)]));
        maxSidebandCoupling = max(maxSidebandCoupling, maxBareCoupling * harmonicWeight);
    end
end

report.deltaF = deltaF;
report.maxBareCoupling = maxBareCoupling;
report.maxSidebandCoupling = maxSidebandCoupling;
report.maxBetaDot = maxBetaDot;
report.deltaFOverSidebandCoupling = deltaF / maxSidebandCoupling;
report.deltaFOverBetaDot = deltaF / maxBetaDot;
report.mMax = mMax;
report.numTimeSamples = numel(tSamples);
end
