%% Plot_Fig_S4.m
% Load QFI_results.mat and plot Figure S4: (a) B_x(t) for four parameters,
% (b) max of B_x(t) with segmented colors.

% clear; clc; close all;
%% Load data
% load('QFI_results.mat', 't', 'Nt', ...
%     'quantity_lam', 'quantity_om', 'quantity_ga', 'quantity_t', ...
%     'ratio_var_dLam', ...
%     'quantity_lam_max', 'quantity_om_max', 'quantity_ga_max', 'quantity_t_max');

%% Plotting parameters (adjustable)
figWidth = 1600;
figHeight = 700;
figPos = [50, 80, figWidth, figHeight];

axisFontSize = 22;
labelFontSize = 32;
legendFontSize = 29;
letterFontSize = 30;

letterPosLeft = [0.02, 0.96];
letterPosRight = [0.02, 0.96];

lineWidthMain = 2;
lineWidthMax = 1.8;

% Colors for parameters
color_lam = [228,26,28]./255;   % red
color_om  = [55,126,184]./255;  % blue
color_ga  = [77,175,74]./255;   % green
color_t   = [152,78,163]./255;  % purple

% Left subplot limits
xlimLeft = [0 30]; xticksLeft = [0 10 20 30];
ylimLeft = [1e-12 1e12]; yticksLeft = [1e-12 1e-6 1e0 1e6 1e12];

% Right subplot limits
xlimRight = [0 30]; xticksRight = [0 10 20 30];
ylimRight = [1e-2 1e12]; yticksRight = [1e-2 1e5 1e12];

%% Create figure
figure('Position', figPos, 'Color', 'w', 'Name', 'QFI Sensitivity Analysis');

% --- Left subplot: four B_x(t) curves ---
subplot(1, 2, 1);
semilogy(t, quantity_t, '-', 'Color', color_t, 'LineWidth', lineWidthMain); hold on;
semilogy(t, quantity_ga, '-', 'Color', color_ga, 'LineWidth', lineWidthMain);
semilogy(t, quantity_om, '-', 'Color', color_om, 'LineWidth', lineWidthMain);
semilogy(t, quantity_lam, '-', 'Color', color_lam, 'LineWidth', lineWidthMain);

set(gca, 'FontSize', axisFontSize, 'LineWidth', 1.5, 'TickLength', [0.020 0.018]);
xlim(xlimLeft); set(gca, 'XTick', xticksLeft);
ylim(ylimLeft); set(gca, 'YTick', yticksLeft);

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel('$\mathcal{B}_x(t)$', 'Interpreter', 'latex', 'FontSize', labelFontSize-2);

h1 = legend({'$x=t$','$x=\gamma$','$x=\omega$','$x=\lambda$'}, ...
    'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', legendFontSize);
set(h1, 'Box', 'off');

text(letterPosLeft(1), letterPosLeft(2), '(a)', 'Units', 'normalized', ...
    'FontSize', letterFontSize, 'FontWeight', 'bold', 'Interpreter', 'latex');

% --- Right subplot: max(B_x) with segmented colors ---
subplot(1, 2, 2);
semilogy(t, quantity_t_max, '-', 'Color', color_t, 'LineWidth', lineWidthMax); hold on;
semilogy(t, quantity_ga_max, '-', 'Color', color_ga, 'LineWidth', lineWidthMax);
semilogy(t, quantity_om_max, '-', 'Color', color_om, 'LineWidth', lineWidthMax);
semilogy(t, quantity_lam_max, '-', 'Color', color_lam, 'LineWidth', lineWidthMax);

set(gca, 'FontSize', axisFontSize, 'LineWidth', 1.5, 'TickLength', [0.020 0.018]);
xlim(xlimRight); set(gca, 'XTick', xticksRight);
ylim(ylimRight); set(gca, 'YTick', yticksRight);

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel('$\max_x\left\{\mathcal{B}_x(t)\right\}$', ...
    'Interpreter', 'latex', 'FontSize', labelFontSize-2);

h2 = legend({'$x=t$','$x=\gamma$','$x=\omega$','$x=\lambda$'}, ...
    'Location', 'northeast', 'Box', 'off', 'Interpreter', 'latex', 'FontSize', legendFontSize);

text(letterPosRight(1), letterPosRight(2), '(b)', 'Units', 'normalized', ...
    'FontSize', letterFontSize, 'FontWeight', 'bold', 'Interpreter', 'latex');

hold off;

fprintf('Plotting complete. \n');