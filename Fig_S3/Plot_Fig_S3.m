%% Plot_Fig_S3.m
% Plot derivatives and scaling law using workspace variables from Fig_S3.m.
% Uses: t, dO_dt, dO_dOm, dO_dLam, dO_dGa, LHS, RHS.

%% Visualization parameters (adjustable)
axisFontSize = 24;
labelFontSize = 36;
letterFontSize = 32;
legendFontSize = 32;
lineWidthVal = 2;

letterPosLeft = [0.85, 0.925];
letterPosRight = [0.935, 0.965];

figWidth = 1600;
figHeight = 700;
figPos = [50, 80, figWidth, figHeight];

% Colors
color_a = [77/255,182/255,172/255];    % dO/dt
color_b = [255/255,111/255,97/255];    % dO/dOm
color_c = [79/255,195/255,247/255];    % dO/dLam
color_d = [255/255,183/255,77/255];    % dO/dGa
color_e_lhs = [249/255,153/255,42/255];   % LHS (line)
color_e_rhs = [42/255,157/255,142/255];   % RHS (scatter)

% Axis limits
xlim_all = [0 30];
xticks_all = [0 10 20 30];

ylim_a = [-0.9 0.9]; yticks_a = [-0.9 -0.3 0.3 0.9];
ylim_b = [-2 1];     yticks_b = [-2 -1 0 1];
ylim_c = [-0.5 1];   yticks_c = [-0.5 0 0.5 1];
ylim_d = [-0.5 0.5]; yticks_d = [-0.5 0 0.5];
ylim_e = [-2 2];     yticks_e = [-2 -1 0 1 2];

scatterStride = 1000;
scatterSize = 12;

%% Create figure
figure('Position', figPos, 'Color', 'w', 'Name', 'H^{(1)}');

% Main layout: 1 row, 2 columns (left 60%, right 40%)
mainLayout = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
mainLayout.InnerPosition = [0.08, 0.12, 0.88, 0.82];

% Left: 2x2 subplots (a,b,c,d)
leftLayout = tiledlayout(mainLayout, 2, 2);
leftLayout.Layout.Tile = 1;
leftLayout.TileSpacing = 'compact';
leftLayout.Padding = 'compact';

% (a) dO/dt
ax_a = nexttile(leftLayout);
plot(t, dO_dt, 'Color', color_a, 'LineWidth', lineWidthVal);
set(gca, 'FontSize', axisFontSize, 'LineWidth', 1, 'TickLength', [0.020, 0.018]);
xlim(xlim_all); set(gca, 'XTick', xticks_all);
ylim(ylim_a); set(gca, 'YTick', yticks_a);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel('$\Delta_t \langle H^{(1)}\rangle (t)$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
text(letterPosLeft(1), letterPosLeft(2), '(a)', 'Units', 'normalized', ...
    'FontSize', letterFontSize, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (b) dO/dOm
ax_b = nexttile(leftLayout);
plot(t, dO_dOm, 'Color', color_b, 'LineWidth', lineWidthVal);
set(gca, 'FontSize', axisFontSize, 'LineWidth', 1, 'TickLength', [0.020, 0.018]);
xlim(xlim_all); set(gca, 'XTick', xticks_all);
ylim(ylim_b); set(gca, 'YTick', yticks_b);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel('$\Delta_\omega \langle H^{(1)}\rangle (t)$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
text(letterPosLeft(1), letterPosLeft(2), '(b)', 'Units', 'normalized', ...
    'FontSize', letterFontSize, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (c) dO/dLam
ax_c = nexttile(leftLayout);
plot(t, dO_dLam, 'Color', color_c, 'LineWidth', lineWidthVal);
set(gca, 'FontSize', axisFontSize, 'LineWidth', 1, 'TickLength', [0.020, 0.018]);
xlim(xlim_all); set(gca, 'XTick', xticks_all);
ylim(ylim_c); set(gca, 'YTick', yticks_c);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel('$\Delta_\lambda \langle H^{(1)}\rangle (t)$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
text(letterPosLeft(1), letterPosLeft(2), '(c)', 'Units', 'normalized', ...
    'FontSize', letterFontSize, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (d) dO/dGa
ax_d = nexttile(leftLayout);
plot(t, dO_dGa, 'Color', color_d, 'LineWidth', lineWidthVal);
set(gca, 'FontSize', axisFontSize, 'LineWidth', 1, 'TickLength', [0.020, 0.018]);
xlim(xlim_all); set(gca, 'XTick', xticks_all);
ylim(ylim_d); set(gca, 'YTick', yticks_d);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel('$\Delta_\gamma \langle H^{(1)}\rangle (t)$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
text(letterPosLeft(1), letterPosLeft(2), '(d)', 'Units', 'normalized', ...
    'FontSize', letterFontSize, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Right: single plot (e) scaling law
ax_e = nexttile(mainLayout, 2);
plot(t, LHS, '-', 'Color', color_e_lhs, 'LineWidth', lineWidthVal); hold on;
idx = 1:scatterStride:Nt;
plot(t(idx), RHS(idx), 'o', 'Color', color_e_rhs, 'LineWidth', 2, 'MarkerSize', scatterSize);
set(gca, 'FontSize', axisFontSize, 'LineWidth', 1, 'TickLength', [0.020, 0.018]);
xlim(xlim_all); set(gca, 'XTick', xticks_all);
ylim(ylim_e); set(gca, 'YTick', yticks_e);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
h = legend('$t\Delta_t \langle H^{(1)}\rangle$', ...
           '$\sum_\theta \theta\Delta_\theta \langle H^{(1)}\rangle$', ...
           'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', legendFontSize);
set(h, 'Box', 'off');
text(letterPosRight(1), letterPosRight(2), '(e)', 'Units', 'normalized', ...
    'FontSize', letterFontSize, 'FontWeight', 'bold', 'Interpreter', 'latex');