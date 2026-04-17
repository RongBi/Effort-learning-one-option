function [r,p] = plotCorrelation(xData,yData)

% Rong Bi 2026

plot(xData,yData,'.k','MarkerSize',12);hold on;
fitlin = robustfit(xData,yData);
plot(xData,fitlin(1)+fitlin(2)*xData,'-','Color','r','LineWidth',2);hold on;
[r,p] = corr(xData,yData); % Pearson (default)
