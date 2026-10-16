%% Plot_Fig_4.m
% Plot comparison of two observables: Var(lambda), Our bound, QCRB.
% Uses workspace variables from Fig_4.m: t, Var1, Bound1, QCRB1, Var2, Bound2, QCRB2.

%% Figure settings (identical to original)
figWidth = 297*1.8;
figHeight = 210*1.75;
figPos = [50, 80, figWidth, figHeight];
colorVar = [0 0 0];                    % black solid
colorBound = [255,127,15]/255;         % orange dotted
colorQCRB = [62,161,68]/255;           % green dash-dot
lineW = 1.5;
xxx = 0.02;
yyy = 0.89;

axisFontSize = 19*0.9;
labelFontSize = 25*0.9;
legendFontSize = 25*0.8;

xlims = [0 30];  xticks = [0 10 20 30];
ylims = [1e-4 1e12];  yticks = [1e-4 1e0 1e4 1e8 1e12];

%% Create figure
figure('Position', figPos, 'Color', 'w');

% --- Subplot (a): observable 1 (omega*sz⊗I) ---
ax1 = subplot(2,1,1);
semilogy(t, Var1, '-', 'Color', colorVar, 'LineWidth', lineW); hold on;
semilogy(t, Bound1, '-', 'Color', colorBound, 'LineWidth', lineW);
semilogy(t, QCRB1, '--', 'Color', colorQCRB, 'LineWidth', lineW);
xlim(xlims); ylim(ylims);
set(gca, 'XTick', xticks, 'XTickLabel', {});
set(gca, 'YTick', yticks, 'FontSize', axisFontSize);
ylabel('Precision', 'Interpreter', 'latex', 'FontSize', labelFontSize);
text(xxx, yyy, '(a)', 'Units', 'normalized', ...
    'Interpreter', 'latex', 'FontSize', legendFontSize, 'FontWeight', 'bold');
hLeg = legend({'$\mathrm{Var}(\lambda)$', 'Our bound', '$\mathrm{QCRB}$'}, ...
    'Interpreter', 'latex', 'FontSize', legendFontSize, 'Location', 'northeast');
set(hLeg, 'Box', 'off');
grid off; box on;

% --- Subplot (b): observable 2 (sx⊗sx) ---
ax2 = subplot(2,1,2);
semilogy(t, Var2, '-', 'Color', colorVar, 'LineWidth', lineW); hold on;
semilogy(t, Bound2, '-', 'Color', colorBound, 'LineWidth', lineW);
semilogy(t, QCRB2, '--', 'Color', colorQCRB, 'LineWidth', lineW);
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