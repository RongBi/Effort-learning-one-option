function [xpos] = bargraph2(y,stdev,data,cols)
% plots a bar graph of y with grey color scheme and suitable x axes and
% errorbars corresponding to std
%
% Miriam C Klein-Flügge
% August 2010

% adapted by Rong Bi 2026

barspacing = 1;

bar1=bar([1:barspacing:1+barspacing*(length(y)-1)],y);
xlim([0.5,length(y)+0.5]);
hold on;

% Two main and two slightly brighter colours
% bule and grey
% cols = {[0 0.4 0.7],[0.6 0.6 0.6],[0.2 0.4+0.2 0.7+0.2],[0.6+0.25 0.6+0.25 0.6+0.25]};
% cols = {[0 0.4 0.7],[0.850 0.325 0.098],[0 0.4 0.7]+0.2,[0.850 0.325 0.098]+0.25};% bule and orange
% cols = {[109, 174, 204]/255,[204, 139, 109]/255,[30, 129, 176]/255,[176, 77, 30]/255};

if nargin < 4
    cols = {[30, 129, 176]/255,[176, 77, 30]/255,[109, 174, 204]/255,[204, 139, 109]/255};
end

set(bar1(1),'FaceColor',cols{1});
set(bar1(2),'FaceColor',cols{2});
set(bar1(1),'EdgeColor','none');
set(bar1(2),'EdgeColor','none');

% Plot individual data points
xpos.left  = []; % nSubj x nNoise (ib)
xpos.right = [];

if nargin>2
    nsub = size(data,1);
    xpos.left  = zeros(nsub, size(y,1));
    xpos.right = zeros(nsub, size(y,1));

    for ib = 1:size(y,1)
        % jitter for per subject
        jit = 0.075*(rand(nsub,1)-0.5);

        xl = (1+barspacing*(ib-1)) - 0.15 + jit; 
        xr = (1+barspacing*(ib-1)) + 0.15 + jit;

        xpos.left(:,ib)  = xl;
        xpos.right(:,ib) = xr;

        scatter(xl, data(:,1,ib), 30, 'filled', ...
            'MarkerFaceColor', cols{3}, 'MarkerEdgeColor', cols{1}, ...
            'MarkerFaceAlpha', .6, 'MarkerEdgeAlpha', .8);

        % scatter high
        scatter(xr, data(:,2,ib), 30, 'filled', ...
            'MarkerFaceColor', cols{4}, 'MarkerEdgeColor', cols{2}, ...
            'MarkerFaceAlpha', .6, 'MarkerEdgeAlpha', .8);

        % scatter(repmat(ib,1,nsub)+.075*(rand(1,nsub)-0.5)-0.15,data(:,1,ib),30,'filled', ...
        %     'MarkerFaceColor',cols{3},'MarkerEdgeColor',cols{1}, ...
        %     'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.8);
        % 
        % scatter(repmat(ib,1,nsub)+.075*(rand(1,nsub)-0.5)+0.15,data(:,2,ib),30,'filled', ...
        %     'MarkerFaceColor',cols{4},'MarkerEdgeColor',cols{2}, ...
        %     'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.8);
    end
end

for i=1:length(y)   
    line([1+barspacing*(i-1)-0.15,1+barspacing*(i-1)-0.15],[y(i,1)-stdev(i,1),y(i,1)+stdev(i,1)],'Color','k','LineWidth',2);
    line([1+barspacing*(i-1)+0.15,1+barspacing*(i-1)+0.15],[y(i,2)-stdev(i,2),y(i,2)+stdev(i,2)],'Color','k','LineWidth',2);
    %line([1+barspacing*(i-1),1+barspacing*(i-1)],[y(i),y(i)+stdev(i)],'Color','k','LineWidth',2);
end
