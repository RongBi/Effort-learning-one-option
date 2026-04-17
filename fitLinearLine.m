function [slope] = fitLinearLine(data)
% This function is to get slope for each subject

% Rong Bi 2026

nSubj = size(data,1); % subj x effort/reward level
nLevels = size(data,2);

slope = NaN(nSubj,1);

for s = 1:nSubj
    lm = fitlm(1:nLevels, data(s,:));

    slope(s) = lm.Coefficients.Estimate(2); % Slope is the second coefficient
end