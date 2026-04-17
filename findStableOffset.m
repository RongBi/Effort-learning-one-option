function [stableOffset_idx] = findStableOffset(forceTrace,stableOnset_idx,smoothKernel,zeroStableThres,threshStableForce)

% define offset of stable period to calculate motor noise
% Rong Bi,2026

[nTrials,~,nSubj] = size(forceTrace);
stableOffset_idx = nan(nSubj,nTrials);

for s = 1:nSubj
    for tr = 1:nTrials
        trace_to_smooth = forceTrace(tr,:,s);
        d_trace = diff(smooth(trace_to_smooth,smoothKernel));

        if ~isnan(stableOnset_idx(s,tr))
            stableOffset_idx(s,tr) = max(intersect(find(d_trace > zeroStableThres), find(trace_to_smooth > threshStableForce)));
        end

    end

end