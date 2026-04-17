function bargraph(y,stdev,data,cols)
% plots a bar graph of y with grey color scheme and suitable x axes and
% errorbars corresponding to std
%
% Miriam C Klein-Flügge
% August 2010

% adapted by Rong Bi 2026

barspacing = 1;

% Two main and two slightly brighter colours
if nargin<4
    % cols = {[0 0.4 0.7],[0.6 0.6 0.6],[0.2 0.4+0.2 0.7+0.2],[0.6+0.25 0.6+0.25 0.6+0.25]};
     cols = {[109, 174, 204]/255,[204, 139, 109]/255,[30, 129, 176]/255,[176, 77, 30]/255};
end

% bar([1:barspacing:1+barspacing*(length(y)-1)],y,'FaceColor',cols{1},'Edgecolor','none','BarWidth',0.4);
bar([1:barspacing:1+barspacing*(length(y)-1)],y,'FaceColor','[0.75 0.75 0.75 0.9]','Edgecolor','k','BarWidth',0.4,'Linewidth',1.5);

xlim([0.5,length(y)+0.5]);
hold on;

% Plot individual data points
if nargin>2
    nsub = size(data,1);
    for ib = 1:length(y)
        % scatter(repmat(ib,1,nsub)+.075*(rand(1,nsub)-0.5),data(:,ib),30,'filled', ...
        %     'MarkerFaceColor',cols{3},'MarkerEdgeColor',cols{1}, ...
        %     'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.8);
        scatter(repmat(ib,1,nsub)+.075*(rand(1,nsub)-0.5),data(:,ib),30,'filled', ...
            'MarkerFaceColor','None','MarkerEdgeColor',[0.5 0.5 0.5], ...
            'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.8);
     end
end

for i=1:length(y)
    line([1+barspacing*(i-1),1+barspacing*(i-1)],[y(i)-stdev(i),y(i)+stdev(i)],'Color','k','LineWidth',2);
    %line([1+barspacing*(i-1),1+barspacing*(i-1)],[y(i),y(i)+stdev(i)],'Color','k','LineWidth',2);
end
