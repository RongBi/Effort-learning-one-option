function [prior_idx,prior_force] = findPrior(forceTrace,startIdx,smoothKernel,zeroPriorThres,forceWindow,threshPlateauForce)
% define prior idx and prior force
% Rong Bi, 2026

[nTrials,~,nSubj] = size(forceTrace);
prior_idx = nan(nSubj,nTrials);
prior_force = nan(nSubj,nTrials);

for s = 1:nSubj
    for tr = 1:nTrials
        trace_to_smooth = forceTrace(tr,:,s);
        d_trace = diff(trace_to_smooth);
        plateau_all_idx = find(abs(smooth(d_trace,smoothKernel)) < zeroPriorThres);

        % remove the first (maybe second) idx if it's before start time or force < 0.15
        if ~isnan(startIdx(s,tr))
            plateauIdxAfterStart = min(find(plateau_all_idx >startIdx(s,tr) & plateau_all_idx >forceWindow & trace_to_smooth(plateau_all_idx)'>threshPlateauForce));

        else
            plateauIdxAfterStart = NaN;
        end

        if ~isnan(plateauIdxAfterStart)
            plateau_all_idx = plateau_all_idx(plateauIdxAfterStart:end);
        else
            plateau_all_idx = [];
        end

        if ~isempty(plateau_all_idx)
            prior_idx(s,tr) = plateau_all_idx(1);
        else
            prior_idx(s,tr) = NaN;
        end

        if ~isnan(prior_idx(s,tr))
            prior_force(s,tr) = mean(trace_to_smooth(prior_idx(s,tr)-forceWindow/2:prior_idx(s,tr)+forceWindow/2));
        else
            prior_force(s,tr) = NaN;
        end
    end

end
