%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code is the main program for analyzing the markers data of dorsal
% skin. The process consists of:
% 1. Initialize and load data
% 2. Compute the normalized neighborhood areas of boundary points
% 3. Compute Spearman correlation matrix
% 4. Community detection
% 5. Visualize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Initialize
clc
close all
clear

%% Load data
load GRF_data.mat
load TRC_data.mat

% Static posture and different speeds
Static_data = TRC_data(800:900,:); % Static
Walk07_data = TRC_data(3000:4500,:); % speed = 0.7m/s
Walk09_data = TRC_data(5900:7400,:); % speed = 0.9m/s
Walk11_data = TRC_data(8900:10000,:); % speed = 1.1m/s
Walk13_data = TRC_data(11500:13000,:); % speed = 1.3m/s
Walk15_data = TRC_data(14700:16000,:); % speed = 1.5m/s
Walk17_data = TRC_data(17200:18500,:); % speed = 1.7m/s

% Choose speed
Static_start = 30 ;
Static_end = 80 ;
Data_ready = [ Static_data(30:80,:) ; Walk13_data ] ;
Static_len = Static_end - Static_start + 1 ;

%% Transform the original data into (t,point_id,3d_position) rensor
num_points = size(Data_ready, 2)/3; % number of points (markers)
num_frames = size(Data_ready, 1);    % total frames
data_3d = zeros(num_frames, num_points, 3);

for t = 1:num_frames
    frame_data = Data_ready(t, :);
    for point = 1:num_points
        col_start = (point-1)*3 + 1;
        col_end = point*3;
        data_3d(t, point, :) = frame_data(col_start:col_end);
    end
    if mod(t, 10) == 0
        fprintf('%d/%d frames have been processed\n', t, num_frames);
    end
end

%% Filtering
cutoff_freq = 8; % cutoff frequency
sampling_rate = 100; % sampling frequency
before_filter_data_3d = data_3d ;

data_3d = lowpass_filter_markers(before_filter_data_3d, cutoff_freq, sampling_rate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2. Compute neighborhood areas

% boundary points
boundary = [1:9:91,5:9:95,2:4,92:94] ;
% inner points
inner = setdiff(1:num_points, boundary) ;

% Compute neighborhood areas of boundary points
Area = zeros(size(data_3d,1), num_points) ;
for i=1:size(data_3d,1)
    data_current = squeeze(data_3d(i, :, :));
    for j=1:num_points
        if ismember(j,boundary)
            Area(i,j) = 1e-5 ;
        else
            S1 = S_threepoint(data_current(j,:),data_current(j-4,:),data_current(j-5,:)) ;
            S2 = S_threepoint(data_current(j,:),data_current(j-5,:),data_current(j+4,:)) ;
            S3 = S_threepoint(data_current(j,:),data_current(j+4,:),data_current(j+5,:)) ;
            S4 = S_threepoint(data_current(j,:),data_current(j+5,:),data_current(j-4,:)) ;
            Area(i,j) = S1 + S2 + S3 + S4 ;
        end
    end
end

% Normalized neighborhood areas
Area_static = mean(Area(1:Static_len,:)) ;
Area_normalized = Area ./ Area_static ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3. Compute Spearman correlation matrix

Area_inner_walk = Area_normalized ;
Area_inner_walk(:,boundary) = [] ;
Area_inner_walk(1:Static_len,:) = [] ;

[G_Spearman, pval] =  corr(Area_inner_walk, Area_inner_walk, 'Type', 'Spearman') ;

% Test the symmetry of Spearman matrix and pval matrix
Asym_GSp = G_Spearman - G_Spearman' ;
Asym_pval = pval - pval' ;
if sum(abs(Asym_GSp(:))) > 1e-6 || sum(abs(Asym_pval(:))) > 1e-6
    warning('Asymmetric Spearman Correlation Matrix or pval Matrix!')
end

adj = ( G_Spearman + G_Spearman' ) / 2 ;
for i=1:length(inner)
    adj(i,i) = 0 ;
end

Hot_GSpearman(adj) ; % Heatmap of the Spearman correlation matrix

%% FDR-test and threshold
alpha = 0.05;   % FDR level
thre = 0.8;    % minimal correlation
mask = FDR_Threshold(alpha, thre, adj, pval, length(inner)) ;
adj(~mask) = 0 ;
adj_pm = adj ;
adj = abs(adj) ;

%% 4. Louvain Community Detection: BCT (Brain Connectivity Toolbox)
[M, Q] = community_louvain(adj) ;
G = graph(adj) ;


%% 5. Plot community detection results

% Sort
tbl = tabulate(M);
[~, sort_idx] = sort(tbl(:, 2), 'descend');  % descending order
sorted_communities = tbl(sort_idx, 1);

% Color map
n_community = sum(tbl(:, 2) > 1);
community_color = zeros(size(tbl, 1), 3);

if n_community <= 8
    colors = lines(n_community) ;
else
    colors = hsv(n_community) ;
end

% red and blue for two largest communities
if n_community >= 1
    colors(1, :) = [0, 0, 1];  % blue
end

if n_community >= 2
    colors(2, :) = [1, 0, 0];  % red
end

% Color for each xommunity
for i = 1:n_community
    community_idx = sorted_communities(i);
    community_color(tbl(:, 1) == community_idx, :) = colors(i, :) ;
end

node_colors = zeros(length(inner),3) ;
for i=1:length(inner)
    if tbl(M(i),2) > 1
        node_colors(i,:) = community_color(M(i),:) ;
    end
end

% Plot the network
figure ;
p = plot(G, 'NodeLabel', inner, 'Layout', 'force', 'UseGravity', true) ;
for i=1:length(inner)
    highlight(p, i, 'NodeColor', node_colors(i,:));
end

for i = 1:length(inner)
    for j = i+1:length(inner)
        if G_Spearman(i,j) > thre
            highlight(p, i, j, 'EdgeColor', 'k', 'LineWidth', 2.0); % black edge for positive weight
        elseif G_Spearman(i,j) < -thre
            highlight(p, i, j, 'EdgeColor', 'r', 'LineWidth', 2.0); % red edge for negative weight
        end
    end
end

%% Real position for communities

% Static posture
frame1_data = squeeze(data_3d(1, :, :)) ;

figure('Position', [100, 100, 400, 800], 'color', 'w');

x = frame1_data(:,1);
y = frame1_data(:,2);
z = frame1_data(:,3);

% Fit the back surface
[Y_grid, Z_grid] = meshgrid(linspace(min(y), max(y), 80), ...
                   linspace(min(z), max(z), 80));
X_fit = griddata(y, z, x, Y_grid, Z_grid, 'v4');
surf(1e-3*X_fit, 1e-3*Y_grid, 1e-3*Z_grid, 'FaceColor', [ 0.8 0.8 0.8 ], 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
hold on ;

scatter3(1e-3*(frame1_data(inner,1)+2), 1e-3*frame1_data(inner,2), 1e-3*frame1_data(inner,3), 60, node_colors, 'filled') ;
hold on ;
scatter3(1e-3*(frame1_data(boundary,1)+2), 1e-3*frame1_data(boundary,2), 1e-3*frame1_data(boundary,3), 30, 'k','filled') ;
scatter3(1e-3*(frame1_data(boundary,1)+2), 1e-3*frame1_data(boundary,2), 1e-3*frame1_data(boundary,3), 5, 'w', 'filled') ;

% View
view(55,15);
xlim(1e-3*[850 940]) ;
xticks(1e-3*[860 900 940]) ;
ylim(1e-3*[440 680]) ;
yticks(1e-3*[500 600]) ;
zlim(1e-3*[1100 1600]) ;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
% axis equal;
grid off ;
hold off ;