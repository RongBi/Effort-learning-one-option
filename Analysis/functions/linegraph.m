function linegraph(data,color)
% add line in bargraph

% Rong Bi 2026

barspacing = 1;
nSubj = size(data,1); % data = nSubj x number of conditions
nBar = size(data,2);
xBarCenters = zeros(1,nBar);

for ib=1:nBar
    xBarCenters(ib) = [1+barspacing*(ib-1)]; % keep consitent with bargraph
end

for s = 1:nSubj
    plot(xBarCenters,data(s,:),'-','Color',color,'LineWidth',2);
end
