%% 绘制矩阵热图
function Hot_GSpearman(G_Spearman)

% 创建红-白-蓝映射
n = 256;
mid = floor(n/2);

% 按照左-中-右进行序号重排
new_order = zeros(1,size(G_Spearman,1)) ;
for i=1:length(new_order)
    if i<=29
        if mod(i,3) == 0
            new_order(i) = (i/3-1) * 7 + 5 ;
        else
            new_order(i) = floor(i/3) * 7 + mod(i,3) ;
        end
    elseif i<=38
        new_order(i) = 6 + (i - 30) * 7 ;
    else
        if mod(i-38,3) == 0
            new_order(i) = (i-38)/3 * 7 ;
        else
            new_order(i) = floor((i-38)/3) * 7 + mod(i-38,3) + 2 ;
        end
    end
end

Rearrange_G = G_Spearman(new_order, new_order) ;

Rearrange_G = G_Spearman ;

% 蓝到白 (负相关部分)
blue_to_white = [linspace(1,0,mid)', linspace(1,0,mid)', linspace(1,0.8,mid)'];
white_to_red = [linspace(1,0.8,mid)', linspace(1,0,mid)', linspace(1,0,mid)'];
redblue_map = [flipud(blue_to_white); white_to_red];

% 绘制热图
figure('color', 'w') ;
h = imagesc(Rearrange_G);
colormap(redblue_map);
colorbar;


% title('Correlation Matrix Heatmap');
xlabel('Marker Sequence');
ylabel('Marker Sequence');
axis square;

% 美化图形
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'XAxisLocation', 'top');
set(gca, 'xticklabel', [], 'xtick', []) ;
set(gca, 'yticklabel', [], 'ytick', []) ;
grid on;