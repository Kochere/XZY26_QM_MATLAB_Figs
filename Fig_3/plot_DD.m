clc; figure;

set(gcf, 'Position', [300, 300, 1200, 450]);

%% 设置子图位置，为中间的colorbar预留空间
% 计算子图的位置：左侧子图
ax1_pos = [0.1, 0.15, 0.35, 0.70];  % [left, bottom, width, height]
% 右侧子图
ax2_pos = [0.6, 0.15, 0.35, 0.70];  % 注意left值调整，为colorbar留出空间

%% 第一个子图
ax1 = subplot(1,2,1);
set(ax1, 'Position', ax1_pos);  % 设置第一个子图位置
colormap(balanced);
imagesc(p0span, f0span, dd1', 'AlphaData', ~isnan(dd1'));
set(gca, 'Color', 'white');

clim([-1e-2, 1e-2]);

ax = gca;
set(ax, 'FontSize', 25, 'LineWidth', 1, ...
    'XGrid', 'off', 'YGrid', 'off', 'ZGrid', 'off', ...
    'GridAlpha', 0.3, 'YDir', 'normal');
set(gca,'ticklength',[0.020 0.018]);
xlim([-0.5, 0.5]); ylim([0, 0.5]);
set(ax, 'XTick', -0.5:0.5:0.5, ...
        'YTick', 0:0.25:0.5);
xlabel('$p_0$', 'FontSize', 45, 'Interpreter', 'latex');
ylabel('$|f_0|$', 'FontSize', 45, 'Interpreter', 'latex');
grid off;

% 绘制边缘线
theta = linspace(0, 2*pi, 1000);
x_edge = 0.5 * cos(theta);
y_edge = 0.5 * abs(sin(theta));
hold on;
plot(x_edge, y_edge, '-', 'Color', [32/255, 162/255, 60/255], 'LineWidth', 2); 
% 添加标记文字
text(0.01, 0.94, '$(a) t=10$', 'Units', 'normalized', 'FontSize', 30, 'Interpreter', 'latex');

%% 第二个子图
ax2 = subplot(1,2,2);
set(ax2, 'Position', ax2_pos);  % 设置第二个子图位置
colormap(balanced);
imagesc(p0span, f0span, dd2', 'AlphaData', ~isnan(dd2'));
set(gca, 'Color', 'white');

clim([-1e-2, 1e-2]);

ax = gca;
set(ax, 'FontSize', 25, 'LineWidth', 1, ...
    'XGrid', 'off', 'YGrid', 'off', 'ZGrid', 'off', ...
    'GridAlpha', 0.3, 'YDir', 'normal');
set(gca,'ticklength',[0.020 0.018]);
xlim([-0.5, 0.5]); ylim([0, 0.5]);
set(ax, 'XTick', -0.5:0.5:0.5, ...
        'YTick', 0:0.25:0.5);
xlabel('$p_0$', 'FontSize', 45, 'Interpreter', 'latex');
grid off;

% 绘制边缘线
theta = linspace(0, 2*pi, 1000);
x_edge = 0.5 * cos(theta);
y_edge = 0.5 * abs(sin(theta));
hold on;
plot(x_edge, y_edge, '-', 'Color', [32/255, 162/255, 60/255], 'LineWidth', 2); 
% 添加标记文字
text(0.01, 0.94, '$(b) t=10^3$', 'Units', 'normalized', 'FontSize', 30, 'Interpreter', 'latex');

%% 创建中间的colorbar
% 计算colorbar位置：在两个子图之间
cbar_pos = [0.515, 0.15, 0.02, 0.70];  % [left, bottom, width, height]
c = colorbar('Position', cbar_pos, 'FontSize', 25);
clim([-1e-2, 1e-2]);
set(c, 'Ticks', [-0.01 0 0.01]);