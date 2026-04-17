function [stableOnset_idx] = findStableOnset(forceTrace,prior_idx,smoothKernel,zeroSumDerivaDiffThres,sampleWindow,sumDerivaThres,threshStableForce)
% define stable onset
% Rong Bi, 2026

[nTrials,~,nSubj] = size(forceTrace);
stableOnset_idx = nan(nSubj,nTrials);

for s = 1:nSubj

    for tr = 1:nTrials
        trace_to_smooth = forceTrace(tr,:,s);
        d_trace = abs(smooth(diff(trace_to_smooth),smoothKernel));

        for p = 1:length(d_trace)-sampleWindow(1)
            sumDeriva(tr,p) = sum(d_trace(p:p+sampleWindow(1)));
        end
             
        %find the time point when sumDerivatives become flat
        d_sumD_trace = abs(smooth(diff(sumDeriva(tr,:)),smoothKernel));

        stable_all_idx = find(d_sumD_trace'<zeroSumDerivaDiffThres & sumDeriva(tr,1:end-1)<=sumDerivaThres(1));

        if ~isnan(prior_idx(s,tr))
            stableIdxAfterPlateau = min(find(stable_all_idx' >prior_idx(s,tr) & trace_to_smooth(stable_all_idx)>threshStableForce));
        else
            stableIdxAfterPlateau = NaN;
        end

        if ~isnan(stableIdxAfterPlateau)
            stable_all_idx = stable_all_idx(stableIdxAfterPlateau:end);
        else
            stable_all_idx =[];
        end

        if ~isempty(stable_all_idx)
            stableOnset_idx(s,tr) = stable_all_idx(1);
        else
            stableOnset_idx(s,tr) = NaN;
        end          

    end

end