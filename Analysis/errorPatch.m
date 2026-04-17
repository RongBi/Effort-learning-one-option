function errorPatch(data,xIndex,color,face)
% plot the standard error as a shaded region for force trace 

% Rong Bi 
% June 2025

% patch is more efficient than fill
mean_data = mean(data(xIndex,:),2,'omitnan')';
sem_data = std(data(xIndex,:), [], 2, 'omitnan')./sqrt(size(data,2));
upper_bound = mean_data + sem_data';
lower_bound = mean_data - sem_data';

x = xIndex;
patch([x, fliplr(x)], [upper_bound, fliplr(lower_bound)], color, 'EdgeColor', 'none', 'FaceAlpha', face);
hold on;
