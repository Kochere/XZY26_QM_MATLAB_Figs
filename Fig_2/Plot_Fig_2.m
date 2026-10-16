%% Plot_Fig_2.m
% Plot Figure 2: variance, our bound, and QCRB for two initial states.
% Uses workspace variables: t, (Var_A, Bound_A, QCRB_A), (Var_B, Bound_B, QCRB_B)

%% Figure settings
figWidth = 297*1.8;
figHeight = 210*1.75;
figPos = [50, 80, figWidth, figHeight];
colorVar = [0 0 0];                    % black solid
colorBound = [255,127,15]/255;         % orange dotted
colorQCRB = [62,161,68]/255;           % green dashed
lineW = 1.5;
xxx = 0.02;
yyy = 0.89;

axisFontSize = 19*0.8;
labelFontSize = 25*0.8;
legendFontSize = 25*0.8;

xlims = [0 20];  xticks = [0 5 10 15 20];
ylims = [1e-4 1e12];  yticks = [1e-4 1e0 1e4 1e8 1e12];

%% Create figure
figure('Position', figPos, 'Color', 'w');

% --- Subplot (a): initial state |g> ---
ax1 = subplot(2,1,1);
semilogy(t, Var_A, '-', 'Color', colorVar, 'LineWidth', lineW); hold on;
semilogy(t, Bound_A, '-', 'Color', colorBound, 'LineWidth', lineW);
semilogy(t, QCRB_A, '--', 'Color', colorQCRB, 'LineWidth', lineW);
xlim(xlims); ylim(ylims);
set(gca, 'XTick', xticks, 'XTickLabel', {});
set(gca, 'YTick', yticks, 'FontSize', axisFontSize);
ylabel('Precision', 'Interpreter', 'latex', 'FontSize', labelFontSize);
text(xxx, yyy, '(a)', 'Units', 'normalized', ...
    'Interpreter', 'latex', 'FontSize', legendFontSize, 'FontWeight', 'bold');
hLeg = legend({'$\mathrm{Var}(\eta)$', 'Our bound', '$\mathrm{QCRB}$'}, ...
    'Interpreter', 'latex', 'FontSize', legendFontSize, 'Location', 'northeast');
set(hLeg, 'Box', 'off');
grid off; box on;

% --- Subplot (b): mixed initial state B ---
ax2 = subplot(2,1,2);
semilogy(t, Var_B, '-', 'Color', colorVar, 'LineWidth', lineW); hold on;
semilogy(t, Bound_B, '-', 'Color', colorBound, 'LineWidth', lineW);
semilogy(t, QCRB_B, '--', 'Color', colorQCRB, 'LineWidth', lineW);
xlim(xlims); ylim(ylims);
set(gca, 'XTick', xticks, 'FontSize', axisFontSize);
set(gca, 'YTick', yticks, 'FontSize', axisFontSize);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel('Precision', 'Interpreter', 'latex', 'FontSize', labelFontSize);
text(xxx, yyy, '(b)', 'Units', 'normalized', ...
    'Interpreter', 'latex', 'FontSize', legendFontSize, 'FontWeight', 'bold');
grid off; box on;

% --- Adjust subplot spacing ---
ax1.Units = 'normalized'; ax2.Units = 'normalized';
ax1.Position(2) = ax1.Position(2) - 0.04;
ax2.Position(2) = ax2.Position(2) + 0.04;
ax2.Position(4) = ax1.Position(4);