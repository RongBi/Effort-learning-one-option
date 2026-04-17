% This is the main script to analyse the force trace for the Effort_learning_oneOption study
% Rong Bi, April 2026

clear; clc;close all;

% Relative path
basepath = fileparts(mfilename('fullpath'));

% Define subfolders
analypath   = fullfile(basepath, 'Analysis');
functionPath = fullfile(basepath, 'Analysis', 'functions');
behapath    = fullfile(basepath, 'Data');
figurePath  = fullfile(basepath, 'Plots');

% Output file
savepath = fullfile(figurePath, 'S1 Data.xlsx');

% Add paths
addpath(genpath(analypath));
addpath(behapath);

cd(analypath);

% Subjects
subjs = {'s02','s05','s06','s07','s09','s10','s11',...
    's12','s13','s14','s15','s17','s19','s20','s21','s22',...
    's23','s24','s25','s26','s27',...
    's28','s30','s32',...
    's33','s34','s35','s36'};
% exclude subjs = {'s03','s04','s08','s16','s18','s29','s31'};

% do low volatility or high volatility block first
subjStaF = {'s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s23','s24','s25','s26','s27','s33','s36'};
subjVolF = {'s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s28','s29','s30','s31','s32','s34','s35'};

% settings
S.nSamples = 2000;
S.samplRate = 500;

nSubj = length(subjs);
S.nTrials = 96;
S.nBlockGroup = 4;
S.nTotalTrials = S.nTrials * S.nBlockGroup;
S.nChanges = 40;

forceTrace = zeros(S.nTotalTrials,S.nSamples,nSubj);
requEff = zeros(S.nTotalTrials,nSubj);
requEffNoN = zeros(S.nTotalTrials,nSubj);
vol = zeros(S.nTotalTrials,nSubj);
noise = zeros(S.nTotalTrials,nSubj);

jumpSize = nan(S.nTotalTrials,nSubj);
jumpSizeNoN = nan(S.nTotalTrials,nSubj);
jumpSign = nan(S.nTotalTrials,nSubj);

requEff_t1 = nan(S.nTotalTrials,nSubj);
requEff_t2 = nan(S.nTotalTrials,nSubj);
requEff_t3 = nan(S.nTotalTrials,nSubj);

nCond = 4;
startRT_M = zeros(nSubj,nCond);
adjustRT_M = zeros(nSubj,nCond);
startRT_log_M = zeros(nSubj,nCond);
adjustRT_log_M = zeros(nSubj,nCond);

adjustRT_NoChags_M = zeros(nSubj,nCond);
adjustRT_log_NoChags_M = zeros(nSubj,nCond);

firstTrialsBlock = 1:(S.nTrials/2):S.nTotalTrials;
idxAllTrials = 1:S.nTotalTrials;
idxAllTrialsNofirst = setdiff(idxAllTrials, firstTrialsBlock);

F.grey = [169 169 169 100]/255;
F.red = [120, 0, 0]/255;

myColors = [0, 0.4470, 0.7410; 0.4940, 0.1840, 0.5560]; % Deep blue, Muted purple

% regression models
glmLables_start = {'interc','vol','noise','last effort','vol*noise','vol*last effort','noise*last effort','vol*noise*last effort'};
regressors_start = {'vol','noise','last effort','vol*noise','vol*last effort','noise*last effort','vol*noise*last effort'};

glmLables_start_fatigue = {'interc','vol','session1-2','session3-4','noiseDiff','trialInSession','last effort','vol*noiseDiff','vol*last effort','noiseDiff*last effort','vol*noiseDiff*last effort'};
regressors_start_fatigue = {'vol','block34-12','block78-56','noiseDiff','trial','last effort','vol*noiseDiff','vol*last effort','noiseDiff*last effort','vol*noiseDiff*last effort'};

glmLables_adjust = {'interc','vol','noise','jumpSign','absJumpSize','vol*noise'};
regressors_adjust = {'vol','noise','jumpSign','absJumpSize','vol*noise'};

condName = {'lowNoise_lowVol','lowNoise_highVol','highNoise_lowVol','highNoise_highVol'};

%% Force trace preprocessing
% load data and save as trials x samples x subjects matrix
for s = 1:nSubj

    [forceNorm,requEffort,requEffortNoN,volatility,noiseBL] = prepareForceTrace(behapath,subjs,s,subjStaF,subjVolF);

    forceTrace(:,:,s) = forceNorm;
    requEff(:,s) = requEffort;
    requEffNoN(:,s) = requEffortNoN;
    vol(:,s) = volatility;
    noise(:,s) = noiseBL;

end

% delete trials had less than 1300 samples, due to recording issues or brief squeezes
[nanCutoffAll,badSubAll,badTrAll] = findBadTrials(forceTrace,subjs);

% detect trials where force began from non-zero values due to accidental squeezing before the thermometer onset 
% set the baseline values as zero to aid data processing
beasline_window_mean = 1:20;
baseline_window = 1:160;
threthHighInitForce = 0.2;
highInitForce_idx = zeros(nSubj,1);
for s = 1:nSubj
    forceTraceSub = forceTrace(:,:,s);
    idx = 0;
    for tr = 1:S.nTotalTrials
        if mean(forceTraceSub(tr,beasline_window_mean),'omitnan') > threthHighInitForce
            idx = idx + 1;
            forceTraceSub(tr,baseline_window)=0;
        end
    end

    forceTrace(:,:,s) = forceTraceSub;

    highInitForce_idx(s,1) = idx; % only affect 5 trials across all subjects

end


%% ========================= Define Force initiation (Start RT) ===========================
% define 'Start' idx
threshStartFit = 0.01;
[startIdx, nTrials_smallerThresh] = findStart(forceTrace, nanCutoffAll, badTrAll,threshStartFit);

% plot to check start RT -> fix start by hand
[startIdx] = fixStartByHand(startIdx);

% delete outlier trials where the start RT was more than 'mean ± three SD" 
for s = 1:nSubj
    mean_start = mean(startIdx(s,:),'omitnan');
    sd_start = std(startIdx(s,:),'omitnan');
    b1 = mean_start-sd_start.*3;
    b2 = mean_start+sd_start.*3;

    for tr = 1:S.nTotalTrials
        if startIdx(s,tr) > b2 || startIdx(s,tr) < b1
            startIdx(s,tr)  = NaN;
        end
    end
end

for s = 1:nSubj
 sumTrials = sum(length(find(~isnan(startIdx(s,:)))));
 % disp([subjs{s} ': ' num2str(sumTrials./S.nTotalTrials)]); % over 97% trials left
end
disp('Preprocessing done!')

%% ========================= Define prior in the first plateau ===========================
% find prior idx and prior force
smoothKernel = 40;
zeroPriorThres = 0.001;
forceWindow = 10; % first plateau
threshPlateauForce = 0.15;
[prior_idx,prior_force] = findPrior(forceTrace,startIdx,smoothKernel,zeroPriorThres,forceWindow,threshPlateauForce);

% fix prior idx by hand
[prior_idx] = fixPriorByHand(prior_idx);

% prior duration: how quickly people arrive at the first plateau
priordura_idx = prior_idx - startIdx;

% delete outlier trials where the prior duration was more than 'mean ± three SD"
for s = 1:nSubj
    mean_prior = mean(priordura_idx(s,:),'omitnan');
    sd_prior = std(priordura_idx(s,:),'omitnan');
    b1 = mean_prior-sd_prior.*3;
    b2 = mean_prior+sd_prior.*3;

    for tr = 1:S.nTotalTrials
        if priordura_idx(s,tr) > b2 || priordura_idx(s,tr) < b1
            prior_idx(s,tr)  = NaN;
        end
    end
end

% update prior duration
priordura_idx = prior_idx - startIdx;

% check number of trials where prior was defined
for s = 1:nSubj
 sumTrials(s,1) = sum(length(find(~isnan(prior_idx(s,:)))));
 % disp('After excluding extreme prior duration')
 % disp([subjs{s} ': ' num2str(sumTrials(s,1)./S.nTotalTrials)]); % over 95.8 % trials left
end

p_sumTrials = sumTrials./S.nTotalTrials;

disp('Prior defined!')

%% ========================= Define onset of Stable period ===========================
smoothKernel = 40;
zeroSumDerivaDiffThres = 0.0001;
sumDerivaThres = 0.1;
sampleWindow = 300;
threshStableForce = 0.1;

% define stable onset 
[stableOnset_idx] = findStableOnset(forceTrace,prior_idx,smoothKernel,zeroSumDerivaDiffThres,sampleWindow,sumDerivaThres,threshStableForce);

% fix stable onset by hand
[stableOnset_idx] = fixStableOnsetByHand(stableOnset_idx);


% compute winth-trial adjustment starts from the prior to stable onset
adjustment_samp = stableOnset_idx - prior_idx;% sample points as unit

disp('Stable onset defined!')

%% percentage of trials in which the prior and stable period could be detected for each subject
% could use this information to exclude the bad subjects
nPlateauTrials = zeros(nSubj,1);
nStableTrials = zeros(nSubj,1);

for s = 1:nSubj
    nPlateauTrials(s) = sum(~isnan(prior_idx(s,:)));
    nStableTrials(s) = sum(~isnan(stableOnset_idx(s,:)));
end

plateauTrials_rate = (nPlateauTrials./S.nTotalTrials).*100;
stableTrials_rate = (nStableTrials./S.nTotalTrials).*100; % use % as uni

%% for measures related to effort learning, convert samples to time (ms)
% start RT
startRT = startIdx.*(1000/S.samplRate);% convert samples to time
startRT_log = log(startRT);

% adjustment RT
adjustRT = adjustment_samp.*(1000/S.samplRate);
adjustRT_log = log(adjustRT);

% prior RT 
priorRT = priordura_idx.*(1000/S.samplRate);
priorRT_mean = nan(nSubj,1);
% mean of prior duration
for s = 1:nSubj
    trPriorRT = intersect(idxAllTrialsNofirst,find(~isnan(priorRT(s,:))));

    priorRT_mean(s) = mean(priorRT(s,trPriorRT),'omitnan');
end

disp('Mean prior duration is')
disp(mean(priorRT_mean));
disp('sd prior duration is')
disp(std(priorRT_mean));


%% Get the change trials
% define idx for break
blockBreak = (S.nTrials/2+1):S.nTrials/2:S.nTotalTrials;
trChanges = zeros(nSubj,S.nChanges);
trChanges_t1 = zeros(nSubj,S.nChanges);

% find trial idx for change trial and following trials
for s = 1:nSubj
    % find change trials t
        idxChanges = find(diff(requEffNoN(:,s))~=0) + 1;
        unvalidTr = intersect(idxChanges,blockBreak); % sometimes there is no change before and after break
        idxChanges(ismember(idxChanges, unvalidTr)) = [];

        idxChanges_t1 = idxChanges + 1;

        trChanges(s,:) = idxChanges; % change trial t
        trChanges_t1(s,:) = idxChanges_t1; % trial t+1
end


%% prep for manuscript Fig.3 & Fig.5 - fit glm (linear regression) for start RT and with-trial adjustment RT
diffBlockGroup12 = [-ones(S.nTrials,1);  ones(S.nTrials,1); zeros(S.nTrials*2,1)];
diffBlockGroup34 = [zeros(S.nTrials*2,1); -ones(S.nTrials,1); ones(S.nTrials,1)];
noiseBlockDiff = [-ones(S.nTrials,1);  -ones(S.nTrials,1); ones(S.nTrials,1); ones(S.nTrials,1)];
trialNumber = [1:S.nTrials 1:S.nTrials 1:S.nTrials 1:S.nTrials]';

betaStartGlm = nan(length(glmLables_start),nSubj);
betaStartGlm_fatigue = nan(length(glmLables_start_fatigue),nSubj);
betaAdjustGlm = nan(length(glmLables_adjust),nSubj);
betaAdjustGlm_noChags = nan(length(glmLables_adjust),nSubj);

for s  = 1:nSubj

    % required effort in trial t-1, t-2, t-3
    requEff_t1(2:end,s) = requEff(1:end-1,s);
    requEff_t2(3:end,s) = requEff(1:end-2,s);
    requEff_t3(4:end,s) = requEff(1:end-3,s);

    jumpSize(2:end,s) = requEff(2:end,s) - requEff(1:end-1,s);% (t)-(t-1)
    jumpSizeAbs = abs(jumpSize);% absolute jump size, no sign

    jumpSizeNoN(2:end,s) = requEffNoN(2:end,s)-requEffNoN(1:end-1,s); % jump size for change trials ignoring noise

    % separate down and up for each trial
    jumpSign(jumpSize(:,s)>0,s) = 1; % jump up, effort change from low to high
    jumpSign(jumpSize(:,s)<=0,s) = -1; % jump down, effort change from high to low

    trFromChanges = sort([trChanges(s,:) trChanges_t1(s,:)]);
    trNochanges = setdiff(idxAllTrials,trFromChanges);

%--------------------------------start RT------------------------------------
% start RT for all trials (prep Fig.3A)
trialsStart = intersect(find(~isnan(startRT_log(s,:))), idxAllTrialsNofirst); % all nan start would be nan plateau
designStartTime = zscore([zscore(vol(trialsStart,s)) zscore(noise(trialsStart,s)) zscore(requEff_t1(trialsStart,s)) ...
                   zscore(vol(trialsStart,s)).*zscore(noise(trialsStart,s)) zscore(vol(trialsStart,s)).*zscore(requEff_t1(trialsStart,s)) zscore(noise(trialsStart,s)).*zscore(requEff_t1(trialsStart,s)) ...
                   zscore(vol(trialsStart,s)).*zscore(noise(trialsStart,s)).*zscore(requEff_t1(trialsStart,s))]);
betaStartGlm(:,s) = glmfit(designStartTime,startRT_log(s,trialsStart),'normal');

%--------------------------------start RT: add fatigue regressor into regression (prep Fig.S4A-B)-------------------
designStartGlm_fatigue = zscore([zscore(vol(trialsStart,s)) diffBlockGroup12(trialsStart) diffBlockGroup34(trialsStart) noiseBlockDiff(trialsStart) zscore(trialNumber(trialsStart)) zscore(requEff_t1(trialsStart,s)) ...
                   zscore(vol(trialsStart,s)).*zscore(noiseBlockDiff(trialsStart)) zscore(vol(trialsStart,s)).*zscore(requEff_t1(trialsStart,s)) zscore(noiseBlockDiff(trialsStart)).*zscore(requEff_t1(trialsStart,s)) ...
                   zscore(vol(trialsStart,s)).*zscore(noiseBlockDiff(trialsStart)).*zscore(requEff_t1(trialsStart,s))]);
betaStartGlm_fatigue(:,s) = glmfit(designStartGlm_fatigue,startRT_log(s,trialsStart),'normal');


%-------------------------------adjustment RT-----------------------------------
% adjustment RT for all trials (prep Fig.5A)
idxAdjust = intersect(find(~isnan(adjustRT_log(s,:))), find(~isnan(jumpSign(:,s)))');
trialsAdjust = intersect(idxAdjust, idxAllTrialsNofirst);
designAdjust = zscore([zscore(vol(trialsAdjust,s)) zscore(noise(trialsAdjust,s)) zscore(jumpSign(trialsAdjust,s)) zscore(jumpSizeAbs(trialsAdjust,s))...
                   zscore(vol(trialsAdjust,s)).*zscore(noise(trialsAdjust,s))]);
betaAdjustGlm(:,s) = glmfit(designAdjust,adjustRT_log(s,trialsAdjust),'normal');

% adjustment RT without changes (prep Fig.S6A)
trialsAdjust_noChags = intersect(trNochanges,trialsAdjust);
designAdjust_noChags = zscore([zscore(vol(trialsAdjust_noChags,s)) zscore(noise(trialsAdjust_noChags,s)) zscore(jumpSign(trialsAdjust_noChags,s)) zscore(jumpSizeAbs(trialsAdjust_noChags,s))...
                   zscore(vol(trialsAdjust_noChags,s)).*zscore(noise(trialsAdjust_noChags,s))]);
betaAdjustGlm_noChags(:,s) = glmfit(designAdjust_noChags,adjustRT_log(s,trialsAdjust_noChags),'normal');

end

%% prep: mean start RT, mean prior duration, mean adjustment RT acorss conditions -> to visualize the interaction in regression
% mean of start RT and adjustment RT in low/high volatility across low/high noise environments
for s = 1:nSubj
    c = 0;
    for i = 1:2 % noise
        for j = 1:2 % vol
            c = c+1;

            trCon = intersect(idxAllTrialsNofirst,idxAllTrials(noise(:,s) == i & vol(:,s) == j));
            trChagsCon = intersect(trChanges(s,:), trCon);
            trNoChagsCon = intersect(setdiff(idxAllTrials,trChanges(s,:)), trCon);

            % all trials
            startRT_M(s,c) = mean(startRT(s,trCon),'omitnan');
            startRT_log_M(s,c) = mean(startRT_log(s,trCon),'omitnan');

            adjustRT_M(s,c) = mean(adjustRT(s,trCon),'omitnan');
            adjustRT_log_M(s,c) = mean(adjustRT_log(s,trCon),'omitnan');
            
            % without changes
            adjustRT_NoChags_M(s,c) = mean(adjustRT(s,trNoChagsCon),'omitnan');
            adjustRT_log_NoChags_M(s,c) = mean(adjustRT_log(s,trNoChagsCon),'omitnan');

        end
    end
end


%% manuscript Fig.2A-C, force trace show evidence of learned effort
% manuscript Fig.2A - find start point
s = 26; % example participant
fsize = [0.1 0.1 0.2 0.26];
f2 = figure('color','w');hold on;
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); 
set(gcf,'units','normalized','position',fsize);
ylim([0 0.6]);yticks(0:0.2:0.6);
xticks(0:500:1500);xticklabels(xticks * 2);% convert to ms
xlabel('Time from trial onset (ms)');ylabel('Force (% MVC)');
title('Example trial');

tr = 86; % example change down trial
force_data = forceTrace(tr,:,s); 
plot(force_data,'linewidth',2,'color',myColors(1,:));hold on;
plot(startIdx(s,tr),force_data(startIdx(s,tr)),'+k','MarkerSize',14,'LineWidth', 2);
legend('Actual force','Start');legend boxoff;

picName = ('fig2A');
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f2,output_graph,'Resolution',300);


% manuscript Fig.2B - prior & stable in change down trials
s=26;
fsize = [0.1 0.1 .2 .26];
f2 = figure('color','w');hold on;
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); 
set(gcf,'units','normalized','position',fsize);
ylim([0 0.7]);yticks(0:0.2:0.8);
xticks(0:500:1500);xticklabels(xticks * 2);% convert to ms
xlabel('Time from trial onset (ms)');ylabel('Force (% MVC)');
title('Example change down trial');

tr = 86; % change down 
force_data = forceTrace(tr,:,s); 
plot(force_data,'linewidth',2,'color',myColors(1,:));hold on;
plot(prior_idx(s,tr),force_data(prior_idx(s,tr)),'*','MarkerSize',10,'LineWidth', 2,'color',myColors(1,:));
plot(stableOnset_idx(s,tr),force_data(stableOnset_idx(s,tr)),'o','MarkerSize',10,'LineWidth', 2,'color',myColors(1,:));

yline(requEff(tr,s),'--k','color',myColors(1,:),'linewidth',1.5);
yline(requEff(tr-1,s),'--k','linewidth',1.5);
yline(requEff(tr+1,s),'--k','color',myColors(2,:),'linewidth',1.5);

force_data_t1 = forceTrace(tr+1,:,s); 
plot(force_data_t1,'linewidth',2,'color',myColors(2,:));hold on;
plot(prior_idx(s,tr+1),force_data_t1(prior_idx(s,tr+1)),'*','MarkerSize',10,'color',myColors(2,:),'LineWidth', 2);
plot(stableOnset_idx(s,tr+1),force_data_t1(stableOnset_idx(s,tr+1)),'o','MarkerSize',10,'color',myColors(2,:),'LineWidth', 2);

x1 = startIdx(s,tr); x2 = prior_idx(s,tr); y1 = 0; y2 = force_data(prior_idx(s,tr)); % prior RT
x3 = prior_idx(s,tr); x4 = stableOnset_idx(s,tr); y3 = force_data(prior_idx(s,tr)); y4 = force_data(stableOnset_idx(s,tr)); % prior RT
fill([x1 x2 x2 x1], [y1 y1 y2 y2], ...
      [0.75 0.75 0.75], ...   
      'FaceAlpha', 0.3, ...  % Transparency
      'EdgeColor', 'none'); 

fill([x3 x4 x4 x3], [y3 y3 y4 y4], ...
      [0.5 0.5 0.5], ...      
      'FaceAlpha', 0.3, ...  
      'EdgeColor', 'none'); 
legend('','Prior','Steady state onset','location','best');legend boxoff;

picName = 'fig2B';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f2,output_graph,'Resolution',300);

% manuscript Fig.2C - prior & stable in change up trials
s = 15;% change up
fsize = [0.1 0.1 .2 .26];
f2 = figure('color','w');hold on;
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); 
set(gcf,'units','normalized','position',fsize);
ylim([0 0.7]);yticks(0:0.2:0.8);
xticks(0:500:1500);xticklabels(xticks * 2);% convert to ms
xlabel('Time from trial onset (ms)');ylabel('Force (% MVC)');
title('Example change up trial');

tr = 358; % change up % tr = 14
force_data = forceTrace(tr,:,s); 
plot(force_data,'linewidth',2,'color',myColors(1,:));hold on;
plot(prior_idx(s,tr),force_data(prior_idx(s,tr)),'*','MarkerSize',10,'LineWidth', 2,'color',myColors(1,:));
plot(stableOnset_idx(s,tr),force_data(stableOnset_idx(s,tr)),'o','MarkerSize',10,'LineWidth', 2,'color',myColors(1,:));

yline(requEff(tr,s),'--k','color',myColors(1,:),'linewidth',1.5);
yline(requEff(tr-1,s),'--k','linewidth',1.5);
yline(requEff(tr+1,s),'--k','color',myColors(2,:),'linewidth',1.5);

force_data_t1 = forceTrace(tr+1,:,s); 
plot(force_data_t1,'linewidth',2,'color',myColors(2,:));hold on;
plot(prior_idx(s,tr+1),force_data_t1(prior_idx(s,tr+1)),'*','MarkerSize',10,'color',myColors(2,:),'LineWidth', 2);
plot(stableOnset_idx(s,tr+1),force_data_t1(stableOnset_idx(s,tr+1)),'o','MarkerSize',10,'color',myColors(2,:),'LineWidth', 2);

x1 = startIdx(s,tr); x2 = prior_idx(s,tr); y1 = 0; y2 = force_data(prior_idx(s,tr)); % prior RT
x3 = prior_idx(s,tr); x4 = stableOnset_idx(s,tr); y3 = force_data(prior_idx(s,tr)); y4 = force_data(stableOnset_idx(s,tr)); % prior RT
fill([x1 x2 x2 x1], [y1 y1 y2 y2], ...
      [0.75 0.75 0.75], ...   
      'FaceAlpha', 0.3, ...  % Transparency
      'EdgeColor', 'none'); 

fill([x3 x4 x4 x3], [y3 y3 y4 y4], ...
      [0.5 0.5 0.5], ...      
      'FaceAlpha', 0.3, ...  
      'EdgeColor', 'none'); 
legend('','','','Trial t-1','Trial t','Trial t+1','location','best');legend boxoff;

picName = 'fig2C';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f2,output_graph,'Resolution',300);


%% manuscript Fig.2D - histogram of Success rate for each participant
success_time = nan(nSubj,S.nTotalTrials);
success_rate = nan(nSubj,S.nTotalTrials);
% define successful trial which force stays above required effort > 1000 ms
for s = 1:nSubj
    for tr = 1:S.nTotalTrials
        trace_to_analy = forceTrace(tr,:,s);
        success_time(s,tr) = numel(find(trace_to_analy>=requEff(tr,s))).*1000/S.samplRate;

        if success_time(s,tr) >= 1000
            success_rate(s,tr) = 1;
        else
            success_rate(s,tr) = 0;
        end

    end
end

success_rate_mean = mean(success_rate,2).*100;

% plot histogram for success rate (for manuscript)
fsize = [0.1 0.1 .20 .26];
f2 = figure('color','w');
histogram(success_rate_mean, 'BinWidth', 2,'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
hold on;

xlabel('Success rate (%)');
ylabel('No.participants');
title('Success rates (n = 28)');
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); 
set(gcf,'units','normalized','position',fsize);
xlim([80, 102]);
xticks(80:5:100);

picName = 'fig2D';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f2,output_graph,'Resolution',300);

%% prep manuscript Fig.2E & 2F & 4A - prior accuracy in all trials, in low/high noise, in change down/up trials
prior_changes = cell(nSubj,2);
jumpSize_changes = cell(nSubj,2);
prior_accu_allTrial = nan(nSubj,length(idxAllTrialsNofirst));
prior_accu_LN = nan(nSubj,length(idxAllTrialsNofirst)./2);
prior_accu_HN = nan(nSubj,length(idxAllTrialsNofirst)./2);

for s = 1:nSubj
    % find change trials in down and up
        idxChanges_up = find(diff(requEffNoN(:,s))>0) + 1;
        idxChanges_down = find(diff(requEffNoN(:,s))<0) + 1;

        trChanges_up = intersect(idxAllTrialsNofirst,idxChanges_up);
        trChanges_down = intersect(idxAllTrialsNofirst,idxChanges_down);

        prior_up_subj = nan(length(trChanges_up),1);
        prior_down_subj = nan(length(trChanges_down),1);
        jumpSize_up_sub = nan(length(trChanges_up),1);
        jumpSize_down_sub = nan(length(trChanges_down),1);

        for i = 1:length(trChanges_up)
            prior_up_subj(i) = prior_force(s,trChanges_up(i)) - requEff(trChanges_up(i),s);
            jumpSize_up_sub(i) = jumpSizeNoN(trChanges_up(i),s);
        end

        for j = 1:length(trChanges_down)
            prior_down_subj(j) = prior_force(s,trChanges_down(j)) - requEff(trChanges_down(j),s);
            jumpSize_down_sub(j) = jumpSizeNoN(trChanges_down(j),s);
        end

        prior_changes{s,1} = prior_up_subj;% <0
        prior_changes{s,2} = prior_down_subj;% >0

        jumpSize_changes{s,1} = jumpSize_up_sub;
        jumpSize_changes{s,2} = jumpSize_down_sub;

        prior_accu_allTrial(s,:) = prior_force(s,idxAllTrialsNofirst) - requEff_t1(idxAllTrialsNofirst,s)';
        prior_accu_LN(s,:) = prior_accu_allTrial(s,1:376/2);
        prior_accu_HN(s,:) = prior_accu_allTrial(s,376/2+1:end);
end

% combine all subjects' prior data into one long column vector
prior_up = [];prior_down = [];
prior_accuracy = []; prior_accuracy_LN = []; prior_accuracy_HN = [];
mean_prior_up = nan(nSubj,1);mean_prior_down = nan(nSubj,1);
mean_priorAcc = nan(nSubj,1);mean_priorAcc_LN = nan(nSubj,1);mean_priorAcc_HN = nan(nSubj,1);
sd_priorAcc = nan(nSubj,1);sd_priorAcc_LN = nan(nSubj,1);sd_priorAcc_HN = nan(nSubj,1);

for s = 1:nSubj
    prior_up = [prior_up; prior_changes{s,1}];  % vertically concatenate
    prior_down = [prior_down; prior_changes{s,2}];

    mean_prior_up(s,1) = mean(prior_changes{s,1},'omitnan');
    mean_prior_down(s,1) = mean(prior_changes{s,2},'omitnan');

    prior_accuracy = [prior_accuracy;prior_accu_allTrial(s,:)];
    prior_accuracy_LN = [prior_accuracy_LN; prior_accu_LN(s,:)];
    prior_accuracy_HN = [prior_accuracy_HN; prior_accu_HN(s,:)];

    mean_priorAcc(s,1) = mean(prior_accu_allTrial(s,:),'omitnan');
    mean_priorAcc_LN(s,1) = mean(prior_accu_LN(s,:),'omitnan');
    mean_priorAcc_HN(s,1) = mean(prior_accu_HN(s,:),'omitnan');
    sd_priorAcc_LN(s,1) = std(prior_accu_LN(s,:),'omitnan');
    sd_priorAcc_HN(s,1) = std(prior_accu_HN(s,:),'omitnan');
end

% compare mean or variance for prior accuracy in low or high noise
[h,p,ci,stats] = ttest(mean_priorAcc);
disp('p & t value of prior accuracy:');
disp(p);
disp(stats.tstat);
[h,p,ci,stats] = ttest(mean_priorAcc_LN, mean_priorAcc_HN, 'Tail','both');
disp('p & t value of prior accuracy - mean LN vs HN:');
disp(p);
disp(stats.tstat);
[h,p,ci,stats] = ttest(sd_priorAcc_LN, sd_priorAcc_HN, 'Tail','both');
disp('p & t value of prior accuracy - sd LN vs HN:');
disp(stats.df);
disp(p);
disp(stats.tstat);

% t-test against zero for priors on change down or up trials
[h,p,ci,stats] = ttest(mean_prior_up);
disp('p & t value of Prior bias up:');
disp(p);
disp(stats.tstat);
[h,p,ci,stats] = ttest(mean_prior_down);
disp('p & t value of Prior bias down:');
disp(p);
disp(stats.tstat);


%% manuscript Fig.2E - histogram of learning accuracy (prior - last effort) across all trials
fsize = [0.1 0.1 0.2 0.26];
f2 = figure('color','w');
histogram(prior_accuracy, 'BinWidth', 0.05,'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none');hold on;
xline(0, '--k', 'LineWidth', 1); 
xlabel('Prior(t) - required effort(t-1)');
ylabel('Frequency');
title('Learning accuracy');
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); 
set(gcf,'units','normalized','position',fsize);
xlim([-0.5, 0.5]);
xticks(-0.4:0.2:0.4);

picName = 'fig2E';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f2,output_graph,'Resolution',300);


%% manuscript Fig.2F - priors on change down and change up trials
% Plot histogram of (prior - required effort) for change up and change down trials
fsize = [0.1 0.1 0.2 0.26];
f2 = figure('color','w');
histogram(prior_up, 'BinWidth', 0.05, 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'DisplayName', 'Change up');hold on;
histogram(prior_down, 'BinWidth', 0.05, 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'DisplayName', 'Change down');

xline(0, '--k', 'LineWidth', 1); % Add vertical line at x=0
% legend('Change up','Change down');
xlabel('Prior(t) - required effort(t)');
ylabel('Frequency');
title('Priors on change trials');
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); 
set(gcf,'units','normalized','position',fsize);

picName = 'fig2F1';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f2,output_graph,'Resolution',300);

% bar plot of mean (prior-required effort) in change up and change down
mean_priorBias = [mean_prior_up mean_prior_down];
f2 = figure('color','w');hold on;
% fsize = [0.1 0.1 0.15 0.25]; 
fsize = [0.1 0.1 0.08 0.15];
set(gcf,'units','normalized','position',fsize);
set(gca, 'XTick',1:2,'XTickLabels',{'Up','Down'},'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); 
bargraph(mean(mean_priorBias),std(mean_priorBias,[],1)./sqrt(nSubj),mean_priorBias);

picName = 'fig2F2';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f2,output_graph,'Resolution',300);


%% manuscript Fig.3 - start RT regression & start RT
% Fig 3A - start RT regression for all trials
f3 = figure('color','w');
fsize = [0.1 0.1 0.32 0.32]; 
set(gcf,'units','normalized','position',fsize);
set(gca,'XTick',1:length(regressors_start),'XTickLabels',regressors_start, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); hold on;
title('Effects on log(start RT)');
ylabel('Beta weight');
bargraph(mean(betaStartGlm(2:end,:),2),std(betaStartGlm(2:end,:),[],2)/sqrt(size(betaStartGlm(2:end,:),2)),betaStartGlm(2:end,:)');
clear pAll statsAll;
for i = 1:length(glmLables_start)
    [h,p,ci,stats] = ttest(betaStartGlm(i,:));
    pAll(i) = p;
    statsAll(i) = stats.tstat;
end
disp('p & t vale of start RT:');
disp(pAll);
disp(statsAll);

picName = 'fig3A';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f3,output_graph,'Resolution',300);

% 3B - start RT for all trials in four conditions
fsize = [0.1 0.1 .20 .26];
f3 = figure('color','w');hold on;
set(gca,'XTick',1:2,'XTickLabels',{'Low noise','High noise'},'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k');
set(gcf,'units','normalized','position',fsize);
ylabel('Start RT (ms)');
ylim([200 600]);
data_startRT_M(:,1,:) = [startRT_M(:,[1,3])];% there are a 3D matrix for individual data plot: data(ib,s,1:4)
data_startRT_M(:,2,:) = [startRT_M(:,[2,4])];
xpos = bargraph2([mean(startRT_M(:,1:2));mean(startRT_M(:,3:4))], [std(startRT_M(:,1:2))/sqrt(nSubj);std(startRT_M(:,3:4))/sqrt(nSubj)], data_startRT_M);
linegraph2(data_startRT_M,F.grey,xpos);
legend('Low volatility','High volatility','location','best');legend boxoff;

picName = 'fig3B';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f3,output_graph,'Resolution',300);

% ---- save data -------
T_Fig3A = array2table(betaStartGlm(2:end,:)', 'VariableNames', regressors_start);
T_Fig3A = addvars(T_Fig3A, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
writetable(T_Fig3A, savepath, 'Sheet', 'Fig3A');

T_Fig3B = array2table(startRT_M, 'VariableNames', condName);
T_Fig3B = addvars(T_Fig3B, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
writetable(T_Fig3B, savepath, 'Sheet', 'Fig3B');


%% manuscipt Fig.4A - plot learning accuracy (prior - last effort) in high/low noise
data_LN = prior_accuracy_LN(:);
data_HN = prior_accuracy_HN(:);

% compare distribution of learning accuracy in low vs.high noise
f4 = figure('color','w');
fsize = [0.1 0.1 0.2 0.26];
set(gcf,'units','normalized','position',fsize);
histogram(prior_accuracy_LN, 'BinWidth', 0.05, 'FaceAlpha', 0.4,'FaceColor', myColors(1,:), 'EdgeColor', 'k');hold on;
histogram(prior_accuracy_HN, 'BinWidth', 0.05, 'FaceAlpha', 0.3,'FaceColor', myColors(2,:), 'EdgeColor', 'k');
xline(0, '--k', 'LineWidth', 1); 
xlabel('Prior(t) - required effort(t-1)');
ylabel('Frequency');
title('Learning accuracy');
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k');
xlim([-0.5, 0.5]);
xticks(-0.4:0.2:0.4);
legend('Low noise','High noise');legend boxoff;

picName = 'fig4A1';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f4,output_graph,'Resolution',300);

% mean and sd of prior accuracy
mean_priorACC = [mean_priorAcc_LN mean_priorAcc_HN sd_priorAcc_LN sd_priorAcc_HN];
data_mean_priorACC(:,1,:) = [mean_priorACC(:,[1,3])];% there are a 3D matrix for individual data plot: data(:,1,ib); data(:,2,ib)
data_mean_priorACC(:,2,:) = [mean_priorACC(:,[2,4])];

colsbar = {[0, 0.4470, 0.7410], [0.4940, 0.1840, 0.5560], [0, 0.4470, 0.7410]+0.15,[0.4940, 0.1840, 0.5560]+0.15};%Deep blue,Muted purple
fS2 = figure('color','w');hold on;
fsize = [0.1 0.1 0.1 0.15]; set(gcf,'units','normalized','position',fsize);
set(gca, 'XTick',1:2,'XTickLabels',{'Mean','Variance'},'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); 
xpos = bargraph2([mean(mean_priorACC(:,1:2));mean(mean_priorACC(:,3:4))], [std(mean_priorACC(:,1:2))/sqrt(nSubj);std(mean_priorACC(:,3:4))/sqrt(nSubj)], data_mean_priorACC, colsbar);
linegraph2(data_mean_priorACC, F.grey, xpos);

picName = 'fig4A2';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(fS2,output_graph,'Resolution',300);

%% manuscipt Fig.4B - Prior learning kernel (without change trials)
% effect of required effort on prior
requEffLabels = {'t-3','t-2','t-1'};
betaPlateau_noChags = nan(length(requEffLabels)+1,nSubj,S.nBlockGroup);

for s = 1:nSubj
    % prior ~ beta3 * required eff (t-3) + beta2 * required eff (t-2) + beta1 * required eff (t-1)
    b = 0;
     for i = 1:2 % noise
        for j = 1:2 % vol
            b = b+1;

            idxPriorCon = intersect(idxAllTrialsNofirst, idxAllTrials(noise(:,s) == i & vol(:,s) == j));
            trialsNoChangesCon = intersect(setdiff(idxAllTrials,trChanges(s,:)),idxPriorCon);
            trialsrmnan = intersect(find(~isnan(prior_force(s,:))), find(~isnan(requEff_t3(:,s))));            
            trPriorCon_noChags = intersect(trialsrmnan, trialsNoChangesCon);

            designPlateau_noChags = zscore([zscore(requEff_t3(trPriorCon_noChags,s)) zscore(requEff_t2(trPriorCon_noChags,s)) zscore(requEff_t1(trPriorCon_noChags,s))]);
            betaPlateau_noChags(:,s,b) = glmfit(designPlateau_noChags,prior_force(s,trPriorCon_noChags),'normal');

        end
     end
end

% plot learning kernels
f4 = figure('Color','w');
fsize = [0.1 0.1 0.24 0.26]; set(gcf,'units','normalized','position',fsize);
set(gca,'XTick',1:length(requEffLabels),'XTickLabels',requEffLabels,'XTickLabelRotation',0,'FontSize',14);hold on;
title('Effect of required effort on prior');
ylabel('Beta weight');xlabel('Required effort');

cols = {[173,216,230]/255, [222,184,135]/255, [30, 129, 176]/255, [176, 77, 30]/255};
lines = {'-','-','-','-'};
for c = 1:4
    plot(mean(betaPlateau_noChags(2:end,:,c),2),'Color',cols{c},'LineStyle',lines{c},'LineWidth',2.5);hold on;

end

for c = 1:4
    errorbar(1:length(requEffLabels), mean(betaPlateau_noChags(2:end,:,c),2), std(betaPlateau_noChags(2:end,:,c),[],2)./sqrt(nSubj), ...
             'Color', cols{c}, 'LineStyle', lines{c}, 'LineWidth', 2.5, ...
             'CapSize', 10);
    hold on;
end

xlim([0.8, length(requEffLabels) + 0.2]); 
legend('LowV-LowN','HighV-LowN','LowV-HighN','HighV-HighN','Location','best'); legend boxoff;
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); 
ylim([0, 0.08]);

picName = 'fig4B';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f4,output_graph,'Resolution',300);

% ANOVA 
t = table([1:nSubj]', betaPlateau_noChags(end-2,:,1)', betaPlateau_noChags(end-2,:,2)', betaPlateau_noChags(end-2,:,3)', betaPlateau_noChags(end-2,:,4)', ...
    betaPlateau_noChags(end-1,:,1)', betaPlateau_noChags(end-1,:,2)', betaPlateau_noChags(end-1,:,3)', betaPlateau_noChags(end-1,:,4)', ...
    betaPlateau_noChags(end,:,1)', betaPlateau_noChags(end,:,2)', betaPlateau_noChags(end,:,3)', betaPlateau_noChags(end,:,4)', ...
    'VariableNames', {'Subject', 'A1B1C1', 'A1B2C1', 'A2B1C1', 'A2B2C1', 'A1B1C2', 'A1B2C2', 'A2B1C2', 'A2B2C2', 'A1B1C3', 'A1B2C3', 'A2B1C3', 'A2B2C3'}); % A-noise, B-volatility, C-trial

within = table({'A1';'A1';'A2';'A2';'A1';'A1';'A2';'A2';'A1';'A1';'A2';'A2'}, {'B1';'B2';'B1';'B2';'B1';'B2';'B1';'B2';'B1';'B2';'B1';'B2'}, {'C1';'C1';'C1';'C1';'C2';'C2';'C2';'C2';'C3';'C3';'C3';'C3'}, ...
    'VariableNames', {'Noise', 'Vol','Trial'});
rm = fitrm(t, 'A1B1C1-A2B2C3 ~ 1', 'WithinDesign', within);
ranovatbl = ranova(rm, 'WithinModel', 'Noise*Vol*Trial');
disp(ranovatbl);

multcompare(rm, 'Noise');
multcompare(rm, 'Vol');
multcompare(rm, 'Trial');
multcompare(rm, 'Noise', 'By', 'Trial','ComparisonType', 'bonferroni')

% averaging beta weight across four conditions for previous three trials,t-test against to zero
betaPlateau_noChags_avg = squeeze(mean(betaPlateau_noChags,3));
clear pAll statsAll stats
for i=1:size(betaPlateau_noChags_avg,1)
    [h,p,ci,stats] = ttest(betaPlateau_noChags_avg(i,:));
    pAll(i) = p;
    statsAll(i) = stats.tstat;
end
disp('p & t vale of t-3,t-2,t-1:');
disp(pAll);
disp(statsAll);

%% manuscipt Fig.4C - learning rate coeffi (without change trials)
pe_prior = nan(S.nTotalTrials,nSubj);
update_prior = nan(S.nTotalTrials,nSubj);
alpha_prior = nan(nSubj,S.nBlockGroup);

for s = 1:nSubj
    pe_prior(:,s) = requEff(:,s)- prior_force(s,:)';
    update_prior(2:end,s) = prior_force(s,2:end)' - prior_force(s,1:end-1)';

    % update (t+1) ~ alpha * pe (t)
    b = 0;
    for i = 1:2 % noise
        for j = 1:2 % vol

            b = b+1;
            idxForce_con = intersect(idxAllTrialsNofirst, idxAllTrials(noise(:,s) == i & vol(:,s) == j));
            trForce = intersect(idxForce_con, find(~isnan(update_prior(:,s))));

            % exclude change trial t and following trial t+1 (otherwise the learning rate would be larger than 1)
            idxChanges_t_t1 = sort([trChanges(s,:) trChanges_t1(s,:)]);
            trForce_noChages = intersect(setdiff(idxAllTrials,idxChanges_t_t1),trForce);

            mdl_force_noChages = glmfit(pe_prior(trForce_noChages(1:end-1),s),update_prior(trForce_noChages(2:end),s),'normal');
            alpha_prior(s,b) = mdl_force_noChages(2);

        end
    end
end

%---------% Plot alpha across conditions (without change trials)-----------
f4 = figure('color','w');
fsize = [0.1 0.1 0.2 0.25];set(gcf,'units','normalized','position',fsize);axis square;
set(gca,'XTick',1:2,'XTickLabels',{'Low noise','High noise'},'FontSize',16, 'LineWidth', 1.2,'XColor','k','YColor','k');hold on;
ylabel('Learning rate coefficient');
data_alpha_prior(:,1,:) = [alpha_prior(:,[1,3])];
data_alpha_prior(:,2,:) = [alpha_prior(:,[2,4])];
xpos = bargraph2([mean(alpha_prior(:,1:2));mean(alpha_prior(:,3:4))], [std(alpha_prior(:,1:2))/sqrt(nSubj);std(alpha_prior(:,3:4))/sqrt(nSubj)], data_alpha_prior);
linegraph2(data_alpha_prior,F.grey,xpos);
legend('Low volatility','High volatility','Location','best');legend boxoff;

picName = 'fig4C';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f4,output_graph,'Resolution',300);

% ANOVA
t = table([1:nSubj]', alpha_prior(:,1), alpha_prior(:,2), alpha_prior(:,3), alpha_prior(:,4),...
    'VariableNames', {'Subject', 'A1B1', 'A1B2', 'A2B1', 'A2B2'}); % A-noise, B-volatility
within = table({'A1';'A1';'A2';'A2'}, {'B1';'B2';'B1';'B2'}, ...
    'VariableNames', {'FactorA', 'FactorB'});
rm = fitrm(t, 'A1B1-A2B2 ~ 1', 'WithinDesign', within);
ranovatbl = ranova(rm, 'WithinModel', 'FactorA*FactorB');
disp(ranovatbl);

multcompare(rm, 'FactorA');% noise
multcompare(rm, 'FactorB');% volatility

% ---- save data -------
T_Fig4C = array2table(alpha_prior, 'VariableNames', condName);
T_Fig4C = addvars(T_Fig4C, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
writetable(T_Fig4C, savepath, 'Sheet', 'Fig4C');

%% manuscript Fig.5 - adjustment RT
% fig.5A - regression
f5 = figure('color','w');
fsize = [0.1 0.1 0.3 0.32];set(gcf,'units','normalized','position',fsize);
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); hold on;
title('Effects on log(adjustment RT)');
ylabel('Beta weight');hold on;
ylim([-0.3 0.5]);
set(gca,'XTick',1:length(regressors_adjust),'XTickLabels',regressors_adjust,'FontSize',16);
bargraph(mean(betaAdjustGlm(2:end,:),2),std(betaAdjustGlm(2:end,:),[],2)/sqrt(size(betaAdjustGlm(2:end,:),2)),betaAdjustGlm(2:end,:)');

clear pAll statsAll;
for i=1:length(glmLables_adjust)
    [h,p,ci,stats] = ttest(betaAdjustGlm(i,:));
    pAll(i) = p;
    statsAll(i) = stats.tstat;
end
disp('p & t value of adjust RT:');
disp(pAll);
disp(statsAll);

picName = 'fig5A';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f5,output_graph,'Resolution',300);

% fig.5B - adjustment RT for all trials
fsize = [0.1 0.1 .20 .28];
f5 = figure('color','w');hold on;
set(gca,'XTick',1:2,'XTickLabels',{'Low noise','High noise'},'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k');
set(gcf,'units','normalized','position',fsize);
ylabel('Adjustment RT (ms)');
ylim([0 950]);
data_adjust_RT_M(:,1,:) = [adjustRT_M(:,[1,3])];
data_adjust_RT_M(:,2,:) = [adjustRT_M(:,[2,4])];
xpos = bargraph2([mean(adjustRT_M(:,1:2));mean(adjustRT_M(:,3:4))], [std(adjustRT_M(:,1:2))/sqrt(nSubj);std(adjustRT_M(:,3:4))/sqrt(nSubj)], data_adjust_RT_M);
linegraph2(data_adjust_RT_M,F.grey, xpos);
legend('Low volatility','High volatility','location','best');legend boxoff;

picName = 'fig5B';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(f5,output_graph,'Resolution',300);

% ANOVA
t = table((1:nSubj)', adjustRT_M(:,1), adjustRT_M(:,2), adjustRT_M(:,3), adjustRT_M(:,4), ...
    'VariableNames', {'Subject', 'A1B1', 'A1B2', 'A2B1', 'A2B2'}); % t-1, A-noise, B-volatility
within = table({'A1';'A1';'A2';'A2'}, {'B1';'B2';'B1';'B2'}, ...
    'VariableNames', {'FactorA', 'FactorB'});
rm = fitrm(t, 'A1B1-A2B2 ~ 1', 'WithinDesign', within);
ranovatbl = ranova(rm, 'WithinModel', 'FactorA*FactorB');
disp(ranovatbl);
multcompare(rm, 'FactorB', 'By', 'FactorA') % compare volatility - B
multcompare(rm, 'FactorA', 'By', 'FactorB') % compare noise - A

% ---- save data -------
T_Fig5A = array2table(betaAdjustGlm(2:end,:)', 'VariableNames', regressors_adjust);
T_Fig5A = addvars(T_Fig5A, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
writetable(T_Fig5A, savepath, 'Sheet', 'Fig5A');

T_Fig5B = array2table(adjustRT_M, 'VariableNames', condName);
T_Fig5B = addvars(T_Fig5B, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
writetable(T_Fig5B, savepath, 'Sheet', 'Fig5B');

%% ------------------------------------------------------- Supplementary figures -------------------------------------------------
%%---------------------------------------------------------------------------------------------------------------------------------
%%---------------------------------------------------------------------------------------------------------------------------------

%% supplementary fig.S1 - raw force traces in low or high noise with best/middle/worst participant
myColorMap = colormap(sky(96));
repSubjs = {'s15','s27','s34'};
repSubjIdx = find(ismember(subjs,repSubjs));% best,worst,average

for s = repSubjIdx
    fS1 = figure('color','w');hold on;
    fsize = [0.1 0.1 0.20 0.58]; set(gcf,'units','normalized','position',fsize);

    subplot(2,1,1);
    set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k');hold on;
    xlabel('Time from trial onset (ms)');ylabel('Force (% MVC)');title('Low noise');
    xticks(0:500:1500);xticklabels(xticks * 2);% convert to ms
    ylim([-0.05 0.8]);yticks(0:0.2:0.8);
    i = 0;
    for tr = 1:96
        i = i +1;
        force_data = forceTrace(tr,:,s);
        if ~isnan(prior_idx(s,tr))
            plot(force_data,'linewidth',2,'color',myColorMap(i,:));hold on;
        end

    end

    subplot(2,1,2);hold on;
    xlabel('Time from trial onset (ms)');ylabel('Force (% MVC)');title('High noise');
    set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k');
    xticks(0:500:1500);xticklabels(xticks * 2);% convert to ms
    ylim([-0.05 0.8]);yticks(0:0.2:1);
    j = 0;
    for tr = 1+192:96+192
        j = j+1;
        force_data = forceTrace(tr,:,s);
        if ~isnan(prior_idx(s,tr))
            plot(force_data,'linewidth',2,'color',myColorMap(j,:));hold on;
        end

    end

    sgtitle(subjs{s},'fontsize',16, 'FontWeight', 'bold');

    picName = ['figS1_',subjs{s}];
    output_graph = fullfile(figurePath, [picName, '.tif']);
    exportgraphics(fS1,output_graph,'Resolution',300);

end


%% supplementary Fig.S2A - average force traces for change down/up trial t and t+1
forceTraceScale = nan(S.nTotalTrials,S.nSamples,nSubj);
nanCutoffConMin = min(nanCutoffAll,[],'all');
xIndex = 1:nanCutoffConMin-1;

allChanges_up = {};allChanges_up_t1 = {};allChanges_down = {};allChanges_down_t1 = {};
priorIdx_changeDown = nan(nSubj,1); priorIdx_changeDown_t1 = nan(nSubj,1); priorIdx_changeUp = nan(nSubj,1); priorIdx_changeUp_t1 = nan(nSubj,1); 

faceArea = 0.5;
bandAlpha = 0.25;

for s = 1:nSubj
    for tr = 1:S.nTotalTrials
        forceTraceScale(tr,:,s) = forceTrace(tr,:,s)./requEff(tr,s);
    end

    % find change trials for down and up
    idxChanges_up = find(diff(requEffNoN(:,s))>0) + 1;
    idxChanges_down = find(diff(requEffNoN(:,s))<0) + 1;

    allChanges_up{s,1} = intersect(idxAllTrialsNofirst,idxChanges_up);
    allChanges_down{s,1} = intersect(idxAllTrialsNofirst,idxChanges_down);
    allChanges_up_t1{s,1} = allChanges_up{s,1} + 1;
    allChanges_down_t1{s,1} = allChanges_down{s,1} + 1;

end

% average force trace for each subject
for s = 1:nSubj
    forceChangeDown(:,s) = mean(forceTraceScale(allChanges_down{s,1},xIndex,s),1,'omitnan')';
    forceChangeDown_t1(:,s) = mean(forceTraceScale(allChanges_down_t1{s,1},xIndex,s),1,'omitnan')';

    forceChangeUp(:,s) = mean(forceTraceScale(allChanges_up{s,1},xIndex,s),1,'omitnan')';
    forceChangeUp_t1(:,s) = mean(forceTraceScale(allChanges_up_t1{s,1},xIndex,s),1,'omitnan')';

    priorIdx_changeDown(s,1) = mean(prior_idx(s,allChanges_down{s,1}),'omitnan');
    priorIdx_changeDown_t1(s,1) = mean(prior_idx(s,allChanges_down_t1{s,1}),'omitnan');

    priorIdx_changeUp(s,1) = mean(prior_idx(s,allChanges_up{s,1}),'omitnan');
    priorIdx_changeUp_t1(s,1) = mean(prior_idx(s,allChanges_up_t1{s,1}),'omitnan');

end

% change down
sf2a = figure('color','w');hold on;
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); 
fsize = [0.1 0.1 0.2 0.26]; set(gcf,'units','normalized','position',fsize);

errorPatch(forceChangeDown,xIndex, myColors(1,:),faceArea);
errorPatch(forceChangeDown_t1,xIndex, myColors(2,:),faceArea);

priorIdxSE(priorIdx_changeDown,forceChangeDown,myColors(1,:),bandAlpha);
priorIdxSE(priorIdx_changeDown_t1,forceChangeDown_t1,myColors(2,:),bandAlpha);

plot(xIndex(:), mean(forceChangeDown(xIndex,:),2,'omitnan'),'linewidth',2,'color',myColors(1,:));hold on;
plot(xIndex(:), mean(forceChangeDown_t1(xIndex,:),2,'omitnan'),'linewidth',2,'color',myColors(2,:));hold on;
yline(1.0,'k--','linewidth',1.5);

xticks(0:500:1500);xticklabels(xticks * 2);% convert to ms
ylim([0 1.6]);yticks(0:0.5:1.5);
ylabel('Force./requEff (% MVC)');
xlabel('Time from trial onset (ms)');%ylabel('Force./requEff (MVC)');
title('Average change down trial');

legend('','','','','Prior at trial t','Prior at trial t+1','Trial t','Trial t+1','location','best');legend boxoff;

picName = 'figS2A_down';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(sf2a,output_graph,'Resolution',300);


% change up 
sf2b = figure('color','w');hold on;
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); 
fsize = [0.1 0.1 0.2 0.26]; set(gcf,'units','normalized','position',fsize);

errorPatch(forceChangeUp,xIndex, myColors(1,:),faceArea);
errorPatch(forceChangeUp_t1,xIndex, myColors(2,:),faceArea);

priorIdxSE(priorIdx_changeUp,forceChangeUp,myColors(1,:),bandAlpha);
priorIdxSE(priorIdx_changeUp_t1,forceChangeUp_t1,myColors(2,:),bandAlpha);

plot(xIndex(:), mean(forceChangeUp(xIndex,:),2,'omitnan'),'linewidth',2,'color',myColors(1,:));hold on;
plot(xIndex(:), mean(forceChangeUp_t1(xIndex,:),2,'omitnan'),'linewidth',2,'color',myColors(2,:));hold on;
yline(1.0,'k--','linewidth',1.5);

xticks(0:500:1500);xticklabels(xticks * 2);% convert to ms
ylim([0 1.6]);yticks(0:0.5:1.5);
ylabel('Force./requEff (% MVC)');
xlabel('Time from trial onset (ms)');%ylabel('Force./requEff (MVC)');
title('Average change up trial');

legend('','','','','Prior at trial t','Prior at trial t+1','Trial t','Trial t+1','location','best');legend boxoff;

picName = 'figS2A_up';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(sf2b,output_graph,'Resolution',300);

%% supplementary Fig.S2B - more example change down and following trials
worstSubjs = {'s07','s14','s19','s24','s27','s30','s32','s33'};
worstSubjIdx = find(ismember(subjs,worstSubjs));
pickTrials = nan(10,2);% subj & trial id
fS2 = figure('color','w');hold on;
fsize = [0.1 0.1 .5 .27];set(gcf,'units','normalized','position',fsize);

while true
    randSubjs = sort(randperm(nSubj,10));
    if sum(ismember(randSubjs,worstSubjIdx)) <=1
        break;
    end
end

for i = 1:length(randSubjs)
    subplot(2,5,i);hold on;
    s = randSubjs(i); %s34 % change down
    tr_changes = allChanges_down{s,1};

    tr = 1;
    while tr == 1 || isnan(prior_idx(s,tr))|| isnan(prior_idx(s,tr+1)) || isnan(stableOnset_idx(s,tr)) || isnan(stableOnset_idx(s,tr+1)) || forceTrace(tr,prior_idx(s,tr),s)-forceTrace(tr+1,prior_idx(s,tr+1),s)<0.06 || forceTrace(tr+1,prior_idx(s,tr+1),s)<0.15
        tr = tr_changes(randi(numel(tr_changes))); % random pick one change trial
    end

    pickTrials(s,1) = s;
    pickTrials(s,2) = tr;

    set(gca, 'LineWidth', 1.2,'fontsize',12,'XColor','k','YColor','k');

    ylim([0 0.7]);yticks(0:0.2:0.8);
    xticks(0:500:1500);xticklabels(xticks * 2);% convert to ms
    title([subjs{s},', trial ',num2str(tr)])

    force_data = forceTrace(tr,:,s);
    plot(force_data,'linewidth',2,'color',myColors(1,:));hold on;
    plot(prior_idx(s,tr),force_data(prior_idx(s,tr)),'*','MarkerSize',6,'LineWidth', 2,'color',myColors(1,:));
    yline(requEff(tr,s),'--k','color',myColors(1,:),'linewidth',1.5);
    yline(requEff(tr-1,s),'--k','linewidth',1.5);

    force_data_t1 = forceTrace(tr+1,:,s);
    plot(force_data_t1,'linewidth',2,'color',myColors(2,:));hold on;
    plot(prior_idx(s,tr+1),force_data_t1(prior_idx(s,tr+1)),'*','MarkerSize',6,'color',myColors(2,:),'LineWidth', 2);
    yline(requEff(tr+1,s),'--k','color',myColors(2,:),'linewidth',1.5);

end

sgtitle('Example change down trial','fontsize',14,'FontWeight', 'bold');

picName = 'figS2B_down';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(fS2,output_graph,'Resolution',300);

%% supplementary Fig.S2B - more example change up and following trials
worstSubjs = {'s14','s24','s27','s30','s33'};
worstSubjIdx = find(ismember(subjs,worstSubjs));
    sf2c = figure('color','w');hold on;
    fsize = [0.1 0.1 .5 .28];
    set(gcf,'units','normalized','position',fsize);
   
    while true
        randSubjs = sort(randperm(nSubj,10));
        if sum(ismember(randSubjs,worstSubjIdx)) <=1
            break;
        end
    end

for i = 1:length(randSubjs)
    subplot(2,5,i);hold on;
    s = randSubjs(i);
    tr_changes = allChanges_up{s,1};

    tr = 1;
    while tr == 1 || isnan(prior_idx(s,tr))|| isnan(prior_idx(s,tr+1)) || isnan(stableOnset_idx(s,tr)) || isnan(stableOnset_idx(s,tr+1)) || forceTrace(tr+1,prior_idx(s,tr+1),s)-forceTrace(tr,prior_idx(s,tr),s)<0.06 || forceTrace(tr,prior_idx(s,tr),s)<0.15 %|| plateau_idx(s,tr)>750
        tr = tr_changes(randi(numel(tr_changes))); % random pick one change trial
    end

    set(gca, 'LineWidth', 1.2,'fontsize',12,'XColor','k','YColor','k');

    ylim([0 0.7]);yticks(0:0.2:0.8);
    xticks(0:500:1500);xticklabels(xticks * 2);% convert to ms
    title([subjs{s},', trial ',num2str(tr)])

    force_data = forceTrace(tr,:,s);
    plot(force_data,'linewidth',2,'color',myColors(1,:));hold on;
    plot(prior_idx(s,tr),force_data(prior_idx(s,tr)),'*','MarkerSize',6,'LineWidth', 2,'color',myColors(1,:));
    yline(requEff(tr,s),'--k','color',myColors(1,:),'linewidth',1.5);
    yline(requEff(tr-1,s),'--k','linewidth',1.5);

    force_data_t1 = forceTrace(tr+1,:,s);
    plot(force_data_t1,'linewidth',2,'color',myColors(2,:));hold on;
    plot(prior_idx(s,tr+1),force_data_t1(prior_idx(s,tr+1)),'*','MarkerSize',6,'color',myColors(2,:),'LineWidth', 2);
    yline(requEff(tr+1,s),'--k','color',myColors(2,:),'linewidth',1.5);

end

 sgtitle('Example change up trial','fontsize',14,'FontWeight', 'bold');

 picName = 'figS2B_up';
 output_graph = fullfile(figurePath, [picName, '.tif']);
 exportgraphics(sf2c,output_graph,'Resolution',300);

%% manuscript fig.S3 - influence of sensorimotor feedback integration on Prior
% (prior - requEff) as a function of jump size bins
splitValue = [0.15 0.25]; % use 0.15/0.25 to make sure we have enough/relatively equal trial for each jump size bin
prior_bias_jumpBin_temp = nan(nSubj,3,2);
jumpSize_bin_temp = nan(nSubj,3,2);

for s = 1:nSubj

    for i = 1:2 % 1=change up, 2=change down
        idx_jumpBin1 = abs(jumpSize_changes{s,i})<=splitValue(1);
        prior_bias_jumpBin_temp(s,1,i) = mean(prior_changes{s,i}(idx_jumpBin1),'omitnan');
        jumpSize_bin_temp(s,1,i) = mean(jumpSize_changes{s,i}(idx_jumpBin1));

        idx_jumpBin2 = abs(jumpSize_changes{s,i})>splitValue(1) & abs(jumpSize_changes{s,i})<splitValue(2);
        prior_bias_jumpBin_temp(s,2,i) = mean(prior_changes{s,i}(idx_jumpBin2),'omitnan');
        jumpSize_bin_temp(s,2,i) = mean(jumpSize_changes{s,i}(idx_jumpBin2));

        idx_jumpBin3 = abs(jumpSize_changes{s,i})>=splitValue(2);
        prior_bias_jumpBin_temp(s,3,i) = mean(prior_changes{s,i}(idx_jumpBin3),'omitnan');
        jumpSize_bin_temp(s,3,i) = mean(jumpSize_changes{s,i}(idx_jumpBin3));

    end

end

% in a order of jump size: -0.3, -0.2, -0.1, 0.1, 0.2,0.3
prior_bias_jumpBin = [prior_bias_jumpBin_temp(:,3,2) prior_bias_jumpBin_temp(:,2,2) prior_bias_jumpBin_temp(:,1,2) prior_bias_jumpBin_temp(:,:,1)];
jumpSize_bin = [jumpSize_bin_temp(:,3,2) jumpSize_bin_temp(:,2,2) jumpSize_bin_temp(:,1,2) jumpSize_bin_temp(:,:,1)];
mean_jumpSize_level = mean(jumpSize_bin,1);

binLevels = cell(1,1);
for k = 1:length(mean_jumpSize_level)
    binLevels{k} = num2str(round(mean_jumpSize_level(k),2));
end

fS3 = figure('Color','w');
fsize = [0.1 0.1 0.26 0.30]; set(gcf,'units','normalized','position',fsize);
bargraph(mean(prior_bias_jumpBin,1),std(prior_bias_jumpBin,[],1)./sqrt(nSubj),prior_bias_jumpBin);hold on;
set(gca, 'XTick',1:length(binLevels),'XTickLabels',binLevels,'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k');
xlabel('Jump size'); 
ylabel('Prior(t) - required effort(t)');
ylim([-0.4 0.4]);
axis square;hold on;title('Priors on change trials');
plot([6 1], [mean_jumpSize_level(1) mean_jumpSize_level(end)],'k-','LineWidth', 1.5);

% plot real slope
lm = fitlm(1:length(binLevels), mean(prior_bias_jumpBin,1));
slope = lm.Coefficients.Estimate(2);
plot(1:length(binLevels), lm.Fitted, '-', 'Color',F.red,'LineWidth', 3);

% mean of slope
slope_prior_changes = zeros(nSubj,1);
for s = 1:nSubj
    lm = fitlm(mean_jumpSize_level, prior_bias_jumpBin(s,:));
    slope_prior_changes(s,1) = lm.Coefficients.Estimate(2); % Slope is the second coefficien
end
slope_mean = mean(slope_prior_changes);

% ideal slope
lm_ideal = fitlm(mean_jumpSize_level, [0.3 0.2 0.1 -0.1 -0.2 -0.3]);
slope_ideal = lm_ideal.Coefficients.Estimate(2);

picName = 'figS3';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(fS3,output_graph,'Resolution',300);

% ---- save data -------
T_FigS3 = array2table(prior_bias_jumpBin, 'VariableNames', binLevels);
T_FigS3 = addvars(T_FigS3, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
writetable(T_FigS3, savepath, 'Sheet', 'FigS3');


%% supplementary Fig.S4 - check fatigue effect in start RT regression
% Fig.S4A, start RT regression, separating noise and fatigue
fS4 = figure('color','w');
fsize = [0.1 0.1 0.28 0.32];
set(gcf,'units','normalized','position',fsize);
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); hold on;
title('Effects on log(start RT)');
ylabel('Beta weight');hold on;
ylim([-0.12 0.2]);
set(gca,'XTick',1:length(regressors_start_fatigue),'XTickLabels',regressors_start_fatigue,'FontSize',16);
bargraph(mean(betaStartGlm_fatigue(2:end,:),2),std(betaStartGlm_fatigue(2:end,:),[],2)/sqrt(size(betaStartGlm_fatigue(2:end,:),2)),betaStartGlm_fatigue(2:end,:)');
clear pAll statsAll;
for i=1:length(glmLables_start_fatigue)
    [h,p,ci,stats] = ttest(betaStartGlm_fatigue(i,:));
    pAll(i) = p;
    statsAll(i) = stats.tstat;
end
disp('p & t vale of start RT:');
disp(pAll);
disp(statsAll);

picName = 'figS4A';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(fS4,output_graph,'Resolution',300);

% Fig.S4B, fit linear regression over trials within each block group to check fatigue (same hand in each block group) 
firstTrialsSession = 1:96:384;
idx = 0;
blockGroupOrder = zeros(S.nTotalTrials,1);
for bg = 1:S.nBlockGroup
    idx = idx + 1;
    blockGroupOrder(firstTrialsSession(bg):firstTrialsSession(bg) + (S.nTrials-1)) = idx;
end

slope_startRT = nan(nSubj, S.nBlockGroup);
for bg = 1:S.nBlockGroup   
    trBlockGroup = intersect(idxAllTrialsNofirst,find(blockGroupOrder == bg));
    [slope] = fitLinearLine(startRT_log(:,trBlockGroup));
    slope_startRT(:,bg) = slope; % save slope for each subject in each blockGroup

end

% plot slope of start RT over trials in each block group
fsize = [0.1 0.1 .24 .28];
fS4 = figure('Name','slope of start RT','color','w');
hold on;
set(gca,'XTick',1:4,'XTickLabels',{'Block 1+2','Block 3+4','Block 5+6','Block 7+8'},'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k');
set(gcf,'units','normalized','position',fsize);

title('Slope of log(start RT) over trials');
slope_startRT = slope_startRT.*1000; % for visulisation

ylabel('Slope per block group (a.u.)');
ylim([-4 6]); 
bargraph(mean(slope_startRT), std(slope_startRT)/sqrt(nSubj), slope_startRT);
linegraph(slope_startRT,F.grey);

picName = 'figS4B';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(fS4,output_graph,'Resolution',300);

% ANOVA for four block groups
t = table((1:nSubj)', slope_startRT(:,1), slope_startRT(:,2), slope_startRT(:,3), slope_startRT(:,4),...
    'VariableNames', {'Subject', 'A1', 'A2', 'A3', 'A4'}); % A-block group
within = table({'A1';'A2';'A3';'A4'}, ...
    'VariableNames', {'FactorA'});
rm = fitrm(t, 'A1-A4 ~ 1', 'WithinDesign', within);
ranovatbl = ranova(rm, 'WithinModel', 'FactorA');
disp(ranovatbl);
multcompare(rm, 'FactorA');

% ---- save data -------
T_FigS4A = array2table(betaStartGlm_fatigue(2:end,:)', 'VariableNames', regressors_start_fatigue);
T_FigS4A = addvars(T_FigS4A, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
writetable(T_FigS4A, savepath, 'Sheet', 'FigS4A');

T_FigS4B = array2table(slope_startRT, 'VariableNames', {'Block 1+2','Block 3+4','Block 5+6','Block 7+8'});
T_FigS4B = addvars(T_FigS4B, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
writetable(T_FigS4B, savepath, 'Sheet', 'FigS4B');

%% supplementary Fig.S6AB - adjustment RT without changes
% Fig.S6A - regression
fS6 = figure('color','w');
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k');hold on;
title('Effects on log(adjustment RT)');
ylabel('Beta weight');hold on;
ylim([-0.2 0.4]);yticks([-40 -20 0 20 40]);
set(gca,'XTick',1:length(regressors_adjust),'XTickLabels',regressors_adjust,'FontSize',16);
fsize = [0.1 0.1 0.25 0.27]; set(gcf,'units','normalized','position',fsize);
bargraph(mean(betaAdjustGlm_noChags(2:end,:),2),std(betaAdjustGlm_noChags(2:end,:),[],2)/sqrt(size(betaAdjustGlm_noChags(2:end,:),2)),betaAdjustGlm_noChags(2:end,:)');
clear pAll statsAll;
for i=1:length(glmLables_adjust)
    [h,p,ci,stats] = ttest(betaAdjustGlm_noChags(i,:));
    pAll(i) = p;
    statsAll(i) = stats.tstat;
end
disp('p & t vale of adjust RT - noChages:');
disp(pAll);
disp(statsAll);

picName = 'figS6A';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(fS6,output_graph,'Resolution',300);

% Fig.S6B - adjustment RT without trials
fS6 = figure('color','w');hold on;
set(gca,'XTick',1:2,'XTickLabels',{'Low noise','High noise'},'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k');
fsize = [0.1 0.1 0.20 0.26]; set(gcf,'units','normalized','position',fsize);
ylabel('Adjustment RT (ms)');
ylim([0 900]);
data_adjust_RT_NoChags_M(:,1,:) = [adjustRT_NoChags_M(:,[1,3])];% there are a 3D matrix for individual data plot: data(ib,s,1:4)
data_adjust_RT_NoChags_M(:,2,:) = [adjustRT_NoChags_M(:,[2,4])];
xpos = bargraph2([mean(adjustRT_NoChags_M(:,1:2));mean(adjustRT_NoChags_M(:,3:4))], [std(adjustRT_NoChags_M(:,1:2))/sqrt(nSubj);std(adjustRT_NoChags_M(:,3:4))/sqrt(nSubj)], data_adjust_RT_NoChags_M);
linegraph2(data_adjust_RT_NoChags_M,F.grey, xpos);
legend('Low volatility','High volatility','location','best');legend boxoff;

picName = 'figS6B';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(fS6,output_graph,'Resolution',300);

% ANOVA
t = table((1:nSubj)', adjustRT_NoChags_M(:,1), adjustRT_NoChags_M(:,2), adjustRT_NoChags_M(:,3), adjustRT_NoChags_M(:,4), ...
    'VariableNames', {'Subject', 'A1B1', 'A1B2', 'A2B1', 'A2B2'}); % A-noise, B-volatility
within = table({'A1';'A1';'A2';'A2'}, {'B1';'B2';'B1';'B2'}, ...
    'VariableNames', {'FactorA', 'FactorB'});
rm = fitrm(t, 'A1B1-A2B2 ~ 1', 'WithinDesign', within);
ranovatbl = ranova(rm, 'WithinModel', 'FactorA*FactorB');
disp(ranovatbl);
multcompare(rm, 'FactorB', 'By', 'FactorA')
multcompare(rm, 'FactorA', 'By', 'FactorB')

% ---- save data -------
T_FigS6A = array2table(betaAdjustGlm_noChags(2:end,:)', 'VariableNames', regressors_adjust);
T_FigS6A = addvars(T_FigS6A, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
writetable(T_FigS6A, savepath, 'Sheet', 'FigS6A');

T_FigS6B = array2table(adjustRT_NoChags_M, 'VariableNames', condName);
T_FigS6B = addvars(T_FigS6B, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
writetable(T_FigS6B, savepath, 'Sheet', 'FigS6B');

%% supplementary Fig.S6CD - adjustment RT in early/late stage
adjustRT_afterChag_M = zeros(nSubj,6);
doChagsVersion = {'include','exclude'};
titleText = {'(including change trials)','(excluding change trials)'};
picName_all = {'figS6C','figS6D'};
sheetName = {'FigS6C','FigS6D'};

for m = 1:length(doChagsVersion)
    for s = 1:nSubj
        idxAfterChags = [];idxAfterChags_late = [];

        for k = 1:size(trChanges(s,:),2)

            if strcmp(doChagsVersion{m},'include')
                idx_temp = trChanges(s,k):trChanges(s,k)+4;% with changes
            elseif strcmp(doChagsVersion{m},'exclude')            
                idx_temp = trChanges(s,k)+1:trChanges(s,k)+4;% no changes
            end

            idxAfterChags = [idxAfterChags idx_temp];
            if vol(trChanges(s,k),s) == 1 % only for stable late stage
                idx_temp_late = trChanges(s,k)+5:trChanges(s,k)+9;
                idxAfterChags_late = [idxAfterChags_late idx_temp_late];
            end
        end

        c = 0;
        for i = 1:2 % noise
            for j = 1:3 % vol
                c = c+1;

                if j == 3 % for stable-late
                    trCon = intersect(idxAfterChags_late,idxAllTrials(noise(:,s) == i));
                else
                    trCon = intersect(idxAfterChags,idxAllTrials(noise(:,s) == i & vol(:,s) == j));
                end

                % early stable vs.vol vs. late stable in low/high noise
                adjustRT_afterChag_M(s,c) = mean(adjustRT(s,trCon),'omitnan');

            end
        end
    end

    fsize = [0.1 0.1 .24 .28];
    fS6 = figure('color','w');hold on;
    set(gca,'XTick',1:2,'XTickLabels',{'Low noise','High noise'},'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k');
    set(gcf,'units','normalized','position',fsize);
    ylabel('Adjustment RT (ms)');
    data_adjustRT_afterChag_M(:,:,1) = [adjustRT_afterChag_M(:,[1,4])]';% there are a 3D matrix for individual data plot: data(ib,s,1:3)
    data_adjustRT_afterChag_M(:,:,2) = [adjustRT_afterChag_M(:,[2,5])]';
    data_adjustRT_afterChag_M(:,:,3) = [adjustRT_afterChag_M(:,[3,6])]';
    bargraph3([mean(adjustRT_afterChag_M(:,1:3));mean(adjustRT_afterChag_M(:,4:6))], [std(adjustRT_afterChag_M(:,1:3))/sqrt(nSubj);std(adjustRT_afterChag_M(:,4:6))/sqrt(nSubj)], data_adjustRT_afterChag_M);
    legend('Low volatility (early)','High volatility','Low volatility (late)','location','best');legend boxoff;

    title(titleText{m});
    picName = picName_all{m};

    output_graph = fullfile(figurePath, [picName, '.tif']);
    exportgraphics(fS6,output_graph,'Resolution',300);

    % ANOVA
    t = table((1:nSubj)', adjustRT_afterChag_M(:,1), adjustRT_afterChag_M(:,2), adjustRT_afterChag_M(:,3),adjustRT_afterChag_M(:,4), adjustRT_afterChag_M(:,5), adjustRT_afterChag_M(:,6),...
        'VariableNames', {'Subject', 'A1B1', 'A1B2', 'A1B3','A2B1','A2B2','A2B3'});
    within = table({'A1';'A1';'A1';'A2';'A2';'A2'}, {'B1';'B2';'B3';'B1';'B2';'B3'}, ...
        'VariableNames', {'FactorA','FactorB'});
    rm = fitrm(t, 'A1B1-A2B3 ~ 1', 'WithinDesign', within);
    ranovatbl = ranova(rm, 'WithinModel', 'FactorA*FactorB');
    disp(ranovatbl);
    multcompare(rm, 'FactorB')
    multcompare(rm, 'FactorB', 'By', 'FactorA','ComparisonType', 'bonferroni')

    %----save data-----
    T_FigS6CD = array2table(adjustRT_afterChag_M, 'VariableNames', {'lowNoise_lowVolEarly','lowNoise_highVol','lowNoise_lowVolLate','highNoise_lowVolEarly','highNoise_highVol','highNoise_lowVolLate'});
    T_FigS6CD = addvars(T_FigS6CD, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
    writetable(T_FigS6CD, savepath, 'Sheet', sheetName{m});

end

close all;

%% supplementary Fig.S5 - motor noise

% ========================= define the offset of stable period ===========================
% find ending point of stable perior for trials going back to zero
zeroStableThres = -0.0001;
smoothKernel = 40;
threshStableForce = 0.1;

stableOffset_idx = findStableOffset(forceTrace,stableOnset_idx,smoothKernel,zeroStableThres,threshStableForce);

% check each force trace and fix by hand
[stableOffset_idx] = fixStableOffsetByHand(stableOffset_idx);

%% supplementary Fig.S5A - example force trace of offset in stable period
s = 26;
fsize = [0.1 0.1 0.2 0.26];
fS5 = figure('color','w');hold on;
set(gca, 'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k'); 
set(gcf,'units','normalized','position',fsize);
ylim([0 0.6]);yticks(0:0.2:0.6);
xticks(0:500:1500);xticklabels(xticks * 2);% convert to ms
xlabel('Time from trial onset (ms)');ylabel('Force (% MVC)');
title('Example trial');

tr = 86; % change down 
force_data = forceTrace(tr,:,s); 
plot(stableOnset_idx(s,tr),force_data(stableOnset_idx(s,tr)),'o','color',[0 0.8 0],'MarkerSize',14,'LineWidth', 2);hold on;
plot(stableOffset_idx(s,tr),force_data(stableOffset_idx(s,tr)),'rs','MarkerSize',14,'LineWidth', 2);
plot(force_data,'linewidth',2,'color',myColors(1,:));hold on;
legend('Onset','Offset');legend boxoff;

picName = ('figS5A');
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(fS5,output_graph,'Resolution',300);


%% supplementary Fig.S5B - motor noise by effort levels

% motor noise is standard deviation in the stable period (from stable onset to stable offset)
motorNoise = nan(nSubj,S.nTotalTrials);

for s = 1:nSubj

    for tr = 1:S.nTotalTrials
        force_data = forceTrace(tr,:,s);

        if ~isnan(stableOnset_idx(s,tr))

            stable_window = stableOnset_idx(s,tr):stableOffset_idx(s,tr);
            stable_force = force_data(stable_window);
            motorNoise(s,tr) = std(stable_force);

        end
    end
end

% plot motor noise by effort levels
effLevels = [0.2 0.3 0.4 0.5]; % use requEff without noise for each effort bin
binLevels = {'0.2','0.3','0.4','0.5'};
motorNoise_effBin = nan(nSubj,numel(effLevels));

for s = 1:nSubj

    for i = 1:numel(effLevels)
        idx_effBin = find(requEffNoN(:,s)==effLevels(i));
        motorNoise_effBin(s,i) = mean(motorNoise(s,idx_effBin),'omitnan');
    end

end

fS5 = figure('Color','w');
fsize = [0.1 0.1 0.2 0.26]; set(gcf,'units','normalized','position',fsize);
bargraph(mean(motorNoise_effBin,1),std(motorNoise_effBin,[],1)./sqrt(nSubj),motorNoise_effBin);hold on;
set(gca, 'XTick',1:length(binLevels),'XTickLabels',binLevels,'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k');
xlabel('Required effort'); 
ylabel('Motor noise');ylim([0 0.035]);
axis square;hold on;

lm = fitlm(1:length(binLevels), mean(motorNoise_effBin,1));
plot(1:length(binLevels), lm.Fitted, '-', 'Color',F.red,'LineWidth', 3);

picName = 'figS5B';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(fS5,output_graph,'Resolution',300);

% stats for slope
slope_motorNoise = zeros(nSubj,1);
for s = 1:nSubj
    lm = fitlm(1:length(binLevels), motorNoise_effBin(s,:));
    slope_motorNoise(s,1) = lm.Coefficients.Estimate(2); % Slope is the second coefficien
end

[h,p,ci,stats] = ttest(slope_motorNoise);
disp('p and t value of motor noise-effort levels slope: ');
disp(p);
disp(stats.tstat);

% ---- save data -------
T_FigS5B = array2table(motorNoise_effBin, 'VariableNames', binLevels);
T_FigS5B = addvars(T_FigS5B, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
writetable(T_FigS5B, savepath, 'Sheet', 'FigS5B');

%% supplementary Fig.S5C - effects on absolute prior updating
% model: prior updating(t) ~ motor noise(t-1) + required effort(t-1)
formula = 'priorUpdate ~ motorNoise_t1 + requEff_t1';
coefName = {'motorNoise_t1','requEff_t1'};
coefLabel = {'Motor noise (t-1)','Required effort (t-1)'};
nCoef = numel(coefName);
coefs = nan(nSubj,nCoef);

for s = 1:nSubj
    motorNoise_sub = motorNoise(s,1:end-1)';
    requEff_sub = requEff(1:end-1,s);
    update_sub = abs(update_prior(2:end,s));

    tr_update = intersect(find(~isnan(requEff_sub)), find(~isnan(update_sub)));
    validTrials = intersect(idxAllTrialsNofirst,find(~isnan(motorNoise_sub)));
    trUsed = intersect(validTrials,tr_update);

    motorNoise_sub_z = zscore(motorNoise_sub(trUsed));
    requEff_sub_z = zscore(requEff_sub(trUsed));
    update_sub_n = update_sub(trUsed);

    data_clean = table(motorNoise_sub_z,requEff_sub_z, update_sub_n,...
        'VariableNames',{'motorNoise_t1','requEff_t1','priorUpdate'});

    mdl = fitglm(data_clean,formula,'Distribution','normal','Options',statset('MaxIter',1000));

    for c = 1:nCoef
        coefs(s,c) = mdl.Coefficients.Estimate(c+1);
    end

end

% plot the beta weights
fS5 = figure('Color','w');
fsize = [0.1 0.1 0.2 0.28]; set(gcf,'units','normalized','position',fsize);
set(gca,'XTick',1:numel(coefLabel),'XTickLabels',coefLabel,'LineWidth', 1.2,'fontsize',16,'XColor','k','YColor','k');hold on;
title('Effects on absolute prior updating');ylabel('Effect size (a.u.)');

bargraph(mean(coefs),std(coefs)./sqrt(nSubj),coefs);
ylim([-0.02 0.04]);

[h,p,ci,stats] = ttest(coefs);
disp('p and t value of prior updating regression: ');
disp(p);
disp(stats.tstat);

picName = 'figS5C';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(fS5,output_graph,'Resolution',300);

% ---- save data -------
T_FigS5C = array2table(coefs, 'VariableNames', coefLabel);
T_FigS5C = addvars(T_FigS5C, subjs(:), 'Before', 1, 'NewVariableNames', 'Subject');
writetable(T_FigS5C, savepath, 'Sheet', 'FigS5C');

%% supplementary Fig.S5D - correlation of motor noise and learning rate
% plot correlation
fS5 = figure('color','w');
fsize = [0.1 0.1 0.34 0.24];set(gcf,'units','normalized','position',fsize);axis square;
mean_motorNoise = mean(motorNoise,2,'omitnan');
mean_alpha_prior = mean(alpha_prior,2);

set(gca,'FontSize',16,'XColor','k','YColor','k','LineWidth',1.2);hold on;
xlabel('Motor noise');ylabel('Learning rate');
plotCorrelation(mean_motorNoise,mean_alpha_prior);hold on;
xlim([0.01 0.025]);ylim([0.2 0.8]);yticks(0.2:0.2:0.8);

[rho,p] = corr(mean_motorNoise,mean_alpha_prior);
disp('p and r value of motor noise-learning rate correlation: ');
disp(p);
disp(rho);

picName = 'figS5D';
output_graph = fullfile(figurePath, [picName, '.tif']);
exportgraphics(fS5,output_graph,'Resolution',300);

close all;


% END OF SCRIPT
