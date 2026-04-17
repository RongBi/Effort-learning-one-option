function linegraph2(data,color,xpos)
% add line in bargraph2

% Rong Bi 2026

if nargin < 3
    color = [[108,117,125]/255, 0.6];
end
hold on;

nSubj = size(data,1);
% Assuming two noise groups (ib = 1 and 2) like your original code
for s = 1:nSubj
    % noise level 1 (ib=1): connect low vs high volatility
    plot([xpos.left(s,1), xpos.right(s,1)], [data(s,1), data(s,2)], '-', ...
        'Color', color, 'LineWidth', 1.5);

    % noise level 2 (ib=2)
    plot([xpos.left(s,2), xpos.right(s,2)], [data(s,3), data(s,4)], '-', ...
        'Color', color, 'LineWidth', 1.5);
end


%----------------------------
% barspacing = 1;
% nSubj = size(data,1); % data = nSubj x number of conditions
% x1BarCenters = zeros(1,2);
% x2BarCenters = zeros(1,2);
% 
% for ib=1:2
%     x1BarCenters(ib) = [1+barspacing*(ib-1)-0.15]; % keep consitent with bargraph2
%     x2BarCenters(ib) = [1+barspacing*(ib-1)+0.15];
% end
% 
% xBarCenters = [x1BarCenters x2BarCenters];
% 
% if nargin<2
%     color = [[108, 117, 125]/255,0.6]; % grey
% end
% 
% for s = 1:nSubj
%     plot(xBarCenters([1,3]),data(s,1:2),'-','Color',color,'LineWidth',2);
%     plot(xBarCenters([2,4]),data(s,3:4),'-','Color',color,'LineWidth',2);
% end
