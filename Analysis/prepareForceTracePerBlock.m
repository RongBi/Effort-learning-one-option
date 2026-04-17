function [nTrials, requiredEffort, requiredEffortNoN, forceNorm] = prepareForceTracePerBlock(data,doHand,nSamples)
% load force data per block, get force trace and required effort
% Rong Bi,2026

nTrials = size(data.result.data,2); % 2=> returns the number of columns
requiredEffort = zeros(nTrials,1);
requiredEffortNoN = zeros(nTrials,1);
force = zeros(nTrials,nSamples);
forceNorm = zeros(nTrials,nSamples);


for tr = 1:nTrials 
    
    MVC = data.result.MVC(doHand);
    requiredEffort(tr) = data.result.params.allEfforts(tr); % required effort in Schedule
    requiredEffortNoN(tr) = data.result.params.allEffortsNoN(tr); % required effort in Schedule without noise
    
    dataStr = ['data',num2str(doHand)];
    force(tr,:) = data.result.data(1,tr).(dataStr); % force data recorded as a absolute griping value    
    forceNorm(tr,:) = force(tr,:)./MVC; % normalize between 0-1
   
end
