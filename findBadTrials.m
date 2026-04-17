function [nanCutoffAll,badSubAll,badTrAll] = findBadTrials(forceTrace,subjs)
% for one option task with 96*4 trials
% detect bad trials with samples less than 1300

% Rong Bi 
% Sep 2025

[nTrials,nSamples,nSubj] = size(forceTrace);
trials = 1:nTrials;

cutoff1 = nan(nSubj,2);
cutoff2 = nan(nSubj,2);
cutoff3 = nan(nSubj,2);
cutoff4 = nan(nSubj,2);

nanCutoffAll = nan(nSubj,1);
badSubAll = nan(nSubj,1);

for s = 1:nSubj

    IDSub = str2double(subjs{s}(end-1:end));
    forceTraceSub = forceTrace(:,:,s);

    % detect bad trials
    nanCutoff = min(find(sum(isnan(forceTraceSub))>0));
    trCutoff = min(find(isnan(forceTraceSub(:,nanCutoff))));
    cutoff1(s,1) = trCutoff;
    cutoff1(s,2) = nanCutoff;

    if nanCutoff < 1300 % check whether there is enough sample for cutoff trial
        badSub = IDSub;
        badTr = trCutoff;% just suitable for 1 bad trial in one block
        goodTr = find(~ismember(trials,badTr));
        nanCutoff = min(find(sum(isnan(forceTraceSub(goodTr,:)))>0));
        trCutoff = min(find(isnan(forceTraceSub(goodTr,nanCutoff)))) + sum(min(find(isnan(forceTraceSub(goodTr,nanCutoff))))>badTr);
        cutoff2(s,1) = trCutoff;
        cutoff2(s,2) = nanCutoff;
        if nanCutoff < 1300
            badTr = [badTr trCutoff];
            goodTr = find(~ismember(trials,badTr));%contain the indices of elements in trials that are not equal to any element in badTr.
            nanCutoff = min(find(sum(isnan(forceTraceSub(goodTr,:)))>0));
            trCutoff = min(find(isnan(forceTraceSub(goodTr,nanCutoff)))) + sum(min(find(isnan(forceTraceSub(goodTr,nanCutoff))))>badTr);
            cutoff3(s,1) = trCutoff;
            cutoff3(s,2) = nanCutoff;

            if nanCutoff < 1300
                badTr = [badTr trCutoff];
                goodTr = find(~ismember(trials,badTr));%contain the indices of elements in trials that are not equal to any element in badTr.
                nanCutoff = min(find(sum(isnan(forceTraceSub(goodTr,:)))>0));
                trCutoff = min(find(isnan(forceTraceSub(goodTr,nanCutoff)))) + sum(min(find(isnan(forceTraceSub(goodTr,nanCutoff))))>badTr);
                cutoff4(s,1) = trCutoff;
                cutoff4(s,2) = nanCutoff;

            else
                cutoff4(s,1) = NaN;
                cutoff4(s,2) = NaN;
            end
        end
    else

        badSub = NaN;
        badTr = NaN;
        cutoff2(s,1) = NaN;
        cutoff2(s,2) = NaN;
    end

    % save bad subs and bad trials and nancutoff
    nanCutoffAll(s,:) = nanCutoff;
    badSubAll(s,:) = badSub;
    badTrAll{s,1} = badTr;

    % add bad trials based on inspection
    if IDSub==4
        unusedTr = [3*96+49 3*96+50 3*96+51]; % b4, tr 49,50,51

        if ~isnan(badSub)
            badTrAll{s,1} = [badTrAll{s,1} unusedTr];
        else
            badTrAll{s,1} = [unusedTr];
            badSubAll(s,:) = IDSub;
        end
    elseif IDSub==18
        unusedTr = [1 1*96+40 1*96+80]; % b2, tr 40,80

        if ~isnan(badSub)
            badTrAll{s,1} = [badTrAll{s,1} unusedTr];
        else
            badTrAll{s,1} = [unusedTr];
            badSubAll(s,:) = IDSub;
        end
    elseif IDSub==21
        unusedTr = [2*96+1];% b3, tr 1

        if ~isnan(badSub)
            badTrAll{s,1} = [badTrAll{s,1} unusedTr];
        else
            badTrAll{s,1} = [unusedTr];
            badSubAll(s,:) = IDSub;
        end
    elseif IDSub==31
        unusedTr = [1*96+36]; % b2, tr 36

        if ~isnan(badSub)
            badTrAll{s,1} = [badTrAll{s,1} unusedTr];
        else
            badTrAll{s,1} = [unusedTr];
            badSubAll(s,:) = IDSub;
        end

    end

end