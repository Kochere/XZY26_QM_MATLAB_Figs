%% Plot_Fig_S1.m
% Plot theoretical, numerical QFI and steady-state value.
% Uses workspace variables from Fig_S1.m: t, QFI_true, QFI_num, F_ss_theory.

%% Figure settings
figure('Position', [100 100 450 350]);
set(gca, 'YDir', 'normal');

% Plot theoretical (solid orange)
plot(t, QFI_true, '-', 'Color', [249/255,153/255,42/255], 'LineWidth', 2); hold on;

% Plot numerical (circles, teal)
plot(t, QFI_num, 'o', 'Color', [42/255,157/255,142/255], 'LineWidth', 2, 'MarkerSize', 8);

% Plot steady-state QFI (dashed black)
yline(F_ss_theory, 'k--', 'LineWidth', 1.5, 'Alpha', 0.5);

% Axes limits and ticks
xlim([0 8]);
set(gca, 'xtick', [0 4 8]);
ylim([0 2]);
set(gca, 'ytick', [0 0.5 1 1.5 2]);

% Font and line settings
set(gca, 'Fontsize', 14, 'linewidth', 1);
set(gca, 'ticklength', [0.020 0.018]);

% Labels
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\mathcal{F}_T(t)$', 'Interpreter', 'latex', 'FontSize', 18);

% Legend
h = legend('Theoretical', 'Numerical', '$\mathcal{F}_{T}^s$');
set(h, 'box', 'off');
set(h, 'location', 'southeast', 'fontsize', 16, 'interpreter', 'latex');