function [forceNorm,requEffort,requEffortNoN,volatility,noiseBL] = prepareForceTrace(behapath,subjs,s,subjStaF,subjVolF)

% Rong Bi 2026

nSamples = 2000;
S.nTrials = 96;
S.nTotalTrials = 96 * 4;

warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');

% load force data
if ismember(subjs{s},subjStaF)
    block(1) = load([behapath,subjs{s},'_h1_lowN_staF_data.mat']); % left hand, low noise
    block(2) = load([behapath,subjs{s},'_h2_lowN_staF_data.mat']); % right hand, low noise
    block(3) = load([behapath,subjs{s},'_h1_highN_staF_data.mat']); % left hand, high noise
    block(4) = load([behapath,subjs{s},'_h2_highN_staF_data.mat']); % right hand, high noise

elseif ismember(subjs{s},subjVolF)
    block(1) = load([behapath,subjs{s},'_h1_lowN_volF_data.mat']); % left hand, low noise
    block(2) = load([behapath,subjs{s},'_h2_lowN_volF_data.mat']); % right hand, low noise
    block(3) = load([behapath,subjs{s},'_h1_highN_volF_data.mat']); % left hand, high noise
    block(4) = load([behapath,subjs{s},'_h2_highN_volF_data.mat']); % right hand, high noise
end

nBlock = size(block,2);

for b = 1:nBlock
    % prepare data
    data = block(b);
    doHand = data.result.data.hand; % 1=left,2=right

    doVolFirst = data.result.params.doVolFirst;

    % function: load force data, compute force and derivative
    [nTrials, requiredEffortB, requiredEffortNoNB, forceNormB] = prepareForceTracePerBlock(data,doHand,nSamples);

    requEffort((b-1)*S.nTrials+1:b*S.nTrials,1) = requiredEffortB;
    requEffortNoN((b-1)*S.nTrials+1:b*S.nTrials,1) = requiredEffortNoNB;

    % scale forceNorm by required effort
    for tr = 1:nTrials
        forceNorm((b-1)*S.nTrials+tr,:) = forceNormB(tr,:);
    end

end

if doVolFirst
    volatility([1:48, 96+1:96+48, 96*2+1:96*2+48, 96*3+1:96*3+48],1) = 2; % high vol,1-48 trials in each block
    volatility([49:96, 96+49:96*2, 96*2+49:96*3, 96*3+49:96*4],1) = 1; % low vol,49-96 trials in each block
else
    volatility([1:48, 96+1:96+48, 96*2+1:96*2+48, 96*3+1:96*3+48],1) = 1; % low vol,1-48 trials in each block
    volatility([49:96, 96+49:96*2, 96*2+49:96*3, 96*3+49:96*4],1) = 2; % high vol,49-96 trials in each block
end

noiseBL(1:S.nTotalTrials/2,1) = 1;% low noise, block1 & block2
noiseBL(S.nTotalTrials/2+1:S.nTotalTrials,1) = 2; % high noise, block3 & block4