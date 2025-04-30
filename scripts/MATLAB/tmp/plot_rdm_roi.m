function plot_rdm_roi(data,c_code)
% data: R by S (R: number of regions, S: number of subjects)
mean_data = mean(data,2);
sem_data = std(data')'/sqrt(12);
% 
% % Sample data
% data = mDissSameFing_wCue;  % Replace this with your actual data
% 
% % Calculate mean and S.E.M.
% mean_data = mean(data);
% sem_data = std(data) / sqrt(size(data, 1));
% 
% Labels for the bars
% labels = {'L-SMA','L-PMv','L-PMd','L-M1','L-S1','L-SPLa','L-SPLp'};
% labels = {'L-SMA','L-PMv','L-PMd','L-M1','L-S1','L-SPLa','L-SPLp','L-DSVC','L-MT+','L-VSVC','L-EAC'};
labels = {'SMA','PMv','PMd','M1','S1','SPLa','SPLp','DSVC','MT+','VSVC','EAC'};

% Create bar plot
% figure;
bar(mean_data,c_code);
hold on;

% Add error bars
errorbar(mean_data, sem_data, 'k', 'linestyle', 'none');

% Customize plot
set(gca, 'XTick', 1:11, 'XTickLabel', labels);
% xlabel('Regions of Interest','Fontsize',14);
ylabel('Mean dissimilarity (a.u.)','Fontsize',14);
% title('Motor sequences within the same cue', 'Fontsize',18);
% title('First finger within the same cue', 'Fontsize',18);
% title('Cue with same motor sequence', 'Fontsize',18);

hold off;
box off;