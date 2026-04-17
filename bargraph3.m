function bargraph3(y,stdev,data)
% plots a bar graph of y with grey color scheme and suitable x axes and
% errorbars corresponding to std
%
% Miriam C Klein-Flügge
% August 2010

% adapted by Rong Bi, Sep 2025

barspacing = 1;

bar1=bar([1:barspacing:1+barspacing*(size(y,1)-1)],y);
xlim([0.5,size(y,1)+0.5]);
hold on;

% set(bar1(1),'FaceColor',[0.5 0.05 0.05]);
% set(bar1(2),'FaceColor',[0.6 0.4 0.4]);
% set(bar1(3),'FaceColor',[0 0.4 0.7]);
% set(bar1(4),'FaceColor',[0 0.3 0.6]);

% cols = {[0.4 0.2 0],[0.6 0.4 0.4],[0.2 0.5 0.8],[0 0.3 0.6]};
% cols = {[135 43 43]/255,[153 102 102]/255,[88 154 198]/255,[9 102 178]/255};
% cols = {[131, 169, 203]/255,[203, 165, 131]/255,[135 43 43]/255,[153 102 102]/255};
% cols = {[109, 174, 204]/255,[204, 139, 109]/255,[30, 129, 176]/255,[176, 77, 30]/255};
cols = {[109, 174, 204]/255,[222, 184, 135]/255,[30, 129, 176]/255,[176, 77, 30]/255};
% colsBr = {[135+40 43+40 43+40]/255,[153+40 102+40 102+40]/255,[88+40 154+40 198+40]/255,[9+40 102+40 178+40]/255};
% colsBr = {[109+40, 174+40, 204+40]/255,[204+40, 139+40, 109+40]/255,[30+40, 129+40, 176+40]/255,[176+40, 77+40, 30+40]/255};
colsBr = {[109+40, 174+40, 204+40]/255,[245,222,179]/255,[30+40, 129+40, 176+40]/255,[176+40, 77+40, 30+40]/255};


for i=1:3, set(bar1(i),'FaceColor',cols{i},'EdgeColor','none'); end

set(bar1(1),'FaceColor',[173,216,230]/255);
set(bar1(2),'FaceColor',[222,184,135]/255);
set(bar1(3),'FaceColor',[30, 129, 176]/255);
% set(bar1(4),'FaceColor',[176, 77, 30]/255);
% set(bar1(1),'FaceColor',[109, 174, 204]/255);
% set(bar1(2),'FaceColor',[204, 139, 109]/255);
% set(bar1(3),'FaceColor',[30, 129, 176]/255);
% set(bar1(4),'FaceColor',[176, 77, 30]/255);


% Plot individual data points
if nargin>2
    nsub = size(data,2);
    for ib = 1:size(y,1)
        
        scatter(repmat(ib,1,nsub)+.04*(rand(1,nsub)-0.5)-0.215,squeeze(data(ib,:,1)),30,'filled', ...
            'MarkerFaceColor',colsBr{1},'MarkerEdgeColor',cols{1}, ...
            'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.8);
 
        scatter(repmat(ib,1,nsub)+.04*(rand(1,nsub)-0.5),squeeze(data(ib,:,2)),30,'filled', ...
            'MarkerFaceColor',colsBr{2},'MarkerEdgeColor',cols{2}, ...
            'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.8);
        
        scatter(repmat(ib,1,nsub)+.04*(rand(1,nsub)-0.5)+0.215,squeeze(data(ib,:,3)),30,'filled', ...
            'MarkerFaceColor',colsBr{3}, 'MarkerEdgeColor',cols{3},...
            'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.8);
        
        % scatter(repmat(ib,1,nsub)+.04*(rand(1,nsub)-0.5)+0.27,squeeze(data(ib,:,4)),30,'filled', ...
        %     'MarkerFaceColor',colsBr{4}, 'MarkerEdgeColor',cols{4},...
        %     'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.8);
    end
end

for i=1:size(y,1)   
    line([1+barspacing*(i-1)-0.215,1+barspacing*(i-1)-0.215],[y(i,1)-stdev(i,1),y(i,1)+stdev(i,1)],'Color','k','LineWidth',2);
    line([1+barspacing*(i-1),1+barspacing*(i-1)],[y(i,2)-stdev(i,2),y(i,2)+stdev(i,2)],'Color','k','LineWidth',2);
    line([1+barspacing*(i-1)+0.215,1+barspacing*(i-1)+0.215],[y(i,3)-stdev(i,3),y(i,3)+stdev(i,3)],'Color','k','LineWidth',2);
    
end