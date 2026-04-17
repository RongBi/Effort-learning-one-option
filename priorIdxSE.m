function priorIdxSE(data_prior,data_force, color,bandAlpha)

% plot the standard error of prior idx as a band
% Rong Bi. 05/09/2025

mu  = mean(data_prior,'omitnan');   
sem = std(data_prior,[],'omitnan') ./sqrt(size(data_prior,1));

ub = mu + sem;  lb = mu - sem;

y0 = 0; ymax = mean(data_force(round(mu),:),'omitnan');

xpoly = [lb ub ub lb];
ypoly = [y0 y0 ymax ymax];

hBand = patch(xpoly, ypoly, color, 'EdgeColor','none', 'FaceAlpha', bandAlpha);
uistack(hBand,'bottom'); 

plot(mu,ymax,'*','MarkerSize',10,'LineWidth', 2,'color',color);
hold on;