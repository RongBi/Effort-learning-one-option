function [nanCutoffAll, badSubIDs, badTrialsAll] = findBadForceTraceTrials(forceTrace, subjs, nSamplesThresh)
% FINDBADFORCETRIALS
%   forceTrace: [nTrials x nSamples x nSubj]
%   subjs:      cellstr of subject IDs (e.g., 'S01', 'S12', ...)
%   minSamples: threshold for a trial to be considered "good" (default 1300)
%
% Returns
%   cutoffPerSubj:   per-subject cutoff (#samples) = min length among GOOD trials
%   badSubIDs:       per-subject numeric ID if any bad trials exist, else NaN
%   badTrialsPerSubj: cell array; each cell has a row vector of bad trial indices

if nargin < 3 || isempty(nSamplesThresh)
    nSamplesThresh = 1300;
end

[~, ~, nSubj] = size(forceTrace);

nanCutoffAll     = nan(nSubj,1);
badSubIDs        = nan(nSubj,1);
badTrialsAll  = cell(nSubj,1);

for s = 1:nSubj
    subID = str2double(subjs{s}(end-1:end));
    forceTraceSub = forceTrace(:,:,s);   % [nTrials x nSamples]

    % Trial length up to the FIRST NaN:
    % cumprod(~isnan)) stays 1 until first NaN, then 0 afterwards; summing gives length.
    nanCutOff = sum(cumprod(~isnan(forceTraceSub), 2), 2);  % [nTrials x 1]

    % Bad trials are those with length < nSamplesThresh
    badTrialsSub = find(nanCutOff < nSamplesThresh);
    badTrialsAll{s} = badTrialsSub;

    if ~isempty(badTrialsSub)
        badSubIDs(s) = subID;
    end

    % Per-subject cutoff = minimum length among GOOD trials (>= minSamples).
    goodCutOff = nanCutOff >= nSamplesThresh;
    if any(goodCutOff)
        nanCutoffAll(s) = min(nanCutOff(goodCutOff));
    else
        nanCutoffAll(s) = NaN;  % or 0, if you prefer
    end
end

end
