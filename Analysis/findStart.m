function [startIdx, nTrials_smallerThresh] = findStart(forceTrace, nanCutoffAll, badTrAll,threshStartFit)
% This script is to define the start time
% Rong Bi, 2026

[nTrials,~,nSubj] = size(forceTrace);

startIdx = nan(nSubj,nTrials);
nTrials_smallerThresh = zeros(nSubj,1);

for s = 1:nSubj
    badTr = badTrAll{s,1};
    nanCutoff = nanCutoffAll(s);
    forceTraceSub = forceTrace(:,:,s);
    idx_slow = 0;

    for tr = 1:nTrials

        if ismember(tr,badTr)
            continue
        else

            yData = forceTraceSub(tr,1:nanCutoff-1);

            % find start idx
            if isempty(find(yData > threshStartFit))
                startIdx(s,tr) = NaN;
            else
                sm_trace = smooth(yData,40);
                sm_trace = sm_trace - mean(sm_trace(1:10));
                d_trace = diff(sm_trace);
                [~,idx_max] = max(d_trace);

                if idx_max > 400   % trials where the maximum derivative occurred later than 400 samples            
                    threshold = 0.0001;
                    idx_slow = idx_slow + 1;
                else 
                    threshold = 0.001;
                end

                if ~isnan(min(intersect(find(d_trace > threshold), find(sm_trace > threshStartFit))))
                    startIdx(s,tr) = min(intersect(find(d_trace > threshold), find(sm_trace > threshStartFit)));
                else
                    startIdx(s,tr) = min(intersect(find(d_trace > 0.0001), find(sm_trace > threshStartFit)));
                end
                
            end

        end

    end

    nTrials_smallerThresh(s,:) = idx_slow;

end
