% Create the figure
figure('Position', [100, 100, 1200, 800], 'color', 'w');

% Upper: temporal curve
subplot(2,1,1);
curve_points = [15, 18, 56, 58];
line_styles = {'-', '--', ':', '-.'};

left_points = curve_points(1:2:end);
for i=1:length(left_points)
    point_idx = left_points(i);
    line_style = line_styles{i} ;
    color = [0, 0, 0.5 + 0.5*(i/length(curve_points))];
    plot(200:500,Area_normalized(200:500,point_idx),'Color', color, 'LineStyle', line_style, 'LineWidth', 1.5);
    hold on ;
end

right_points = curve_points(2:2:end);
for i=1:length(right_points)
    point_idx = right_points(i);
    line_style = line_styles{i} ;
    color = [0.5 + 0.5*(i/length(curve_points)), 0, 0];
    plot(200:500,Area_normalized(200:500,point_idx),'Color', color, 'LineStyle', line_style, 'LineWidth', 1.5);
    hold on ;
end

legend('Marker15', 'Marker56', ...
       'Marker18', 'Marker58', 'Location', 'best', 'NumColumns', 2);

% Specific frames
vertical_times = [308, 336, 363];
for i=1:length(vertical_times)
    line_handle = plot([vertical_times(i) vertical_times(i)], [0.9 1.12], ...
    'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    text(vertical_times(i)-2, max(Area_normalized(:)) * 0.97 , sprintf('t=%.2f', vertical_times(i)/100), ...
    'HorizontalAlignment', 'right', 'BackgroundColor', 'white');
end

xlim([200,500]) ;
ylim([0.9,1.12]) ;
grid on ;
hold off;


% Lower: heatmap
subplot(2,1,2);
t_plot = [308, 336, 363];

all_area_data = Area_normalized(:);
global_min = min(all_area_data);
global_max = max(all_area_data);

my_colormap = jet(256) ;

for i=1:length(t_plot)
    t = t_plot(i) ;
    frame_now_data = squeeze(data_3d(t,:,:));

    x = frame_now_data(:,1);
    y = frame_now_data(:,2);
    z = frame_now_data(:,3);
    
    % create subplot for each specific frame
    subplot(2,3,3+i);

    % Fit the back surface
    [Y_grid, Z_grid] = meshgrid(linspace(min(y), max(y), 100), ...
                       linspace(min(z), max(z), 100));
    X_fit = griddata(y, z, x, Y_grid, Z_grid, 'v4');

    color_data = Area_normalized(t,:);

    C_data = griddata(y, z, color_data, Y_grid, Z_grid, 'v4');

    % Plot
    if i==2
        surf(X_fit+5, Y_grid+55, Z_grid, 'CData', C_data, 'FaceAlpha', 0.5);
    elseif i==3
        surf(X_fit, Y_grid+70, Z_grid, 'CData', C_data, 'FaceAlpha', 0.5);
    else
        surf(X_fit, Y_grid, Z_grid, 'CData', C_data, 'FaceAlpha', 0.5);
    end
    shading interp;
    hold on;
    
    
    % Plot points
    if i==2
        scatter3(x+5, y+55, z, 35, color_data, 'filled');
    elseif i==3
        scatter3(x, y+70, z, 35, color_data, 'filled');
    else
        scatter3(x, y, z, 35, color_data, 'filled');
    end
    
    colormap(my_colormap);

    caxis_range = [global_min, global_max];

    % colorbar
    caxis(caxis_range);
    colorbar;
    
    % Specific markers
    if i==2
        scatter3(x(curve_points)+5, y(curve_points)+55, z(curve_points), 75, 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    elseif i==3
        scatter3(x(curve_points), y(curve_points)+70, z(curve_points), 75, 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    else
        scatter3(x(curve_points), y(curve_points), z(curve_points), 75, 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    end
    
    title(sprintf('t=%.2f s', t/100));
    
    grid off ;

    % View
    view(50,20);
    xlim([880 990]) ;
    xticks([900 940 980]) ;
    ylim([470 710]) ;
    zlim([1100 1560]) ;
    xlabel('x'); ylabel('y'); zlabel('z');
end