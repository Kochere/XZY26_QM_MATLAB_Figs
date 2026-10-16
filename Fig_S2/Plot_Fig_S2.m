%% Plot_Fig_S2.m
% Plot QFI for time, omega, lambda, gamma in a 2x2 subplot.
% Uses workspace variables from Fig_S2.m: t, QFI_t, QFI_om, QFI_lam, QFI_ga.

%% Figure settings
figure('Position', [60 80 850 700]);

% Color scheme (consistent with original)
color_t = [77/255,182/255,172/255];
color_om = [255/255,111/255,97/255];
color_lam = [79/255,195/255,247/255];
color_ga = [255/255,183/255,77/255];

% Common font sizes
font_num = 14;
font_axis = 20;
index_abcd = 19;
abcd_x = 0.01;
abcd_y = 0.925;

% (1,1): QFI for time
subplot(2,2,1);
plot(t, QFI_t, 'Color', color_t, 'LineWidth', 2);
xlim([0 30]); set(gca,'xtick',[0 10 20 30]);
ylim([0 6]); set(gca,'ytick',[0 2 4 6]);
set(gca,'Fontsize',font_num,'linewidth',1,'ticklength',[0.020 0.018]);
xlabel('$t$', 'Interpreter','latex','FontSize',font_axis);
ylabel('$\mathcal{F}_{t}(t)$', 'Interpreter','latex','FontSize',font_axis);
text(abcd_x, abcd_y, '(a)', 'Units','normalized','FontSize',index_abcd,...
    'FontWeight','bold','Interpreter','latex');

% (1,2): QFI for omega
subplot(2,2,2);
plot(t, QFI_om, 'Color', color_om, 'LineWidth', 2);
xlim([0 30]); set(gca,'xtick',[0 10 20 30]);
ylim([0 30]); set(gca,'ytick',[0 10 20 30]);
set(gca,'Fontsize',font_num,'linewidth',1,'ticklength',[0.020 0.018]);
xlabel('$t$', 'Interpreter','latex','FontSize',font_axis);
ylabel('$\mathcal{F}_{\omega}(t)$', 'Interpreter','latex','FontSize',font_axis);
text(abcd_x, abcd_y, '(b)', 'Units','normalized','FontSize',index_abcd,...
    'FontWeight','bold','Interpreter','latex');

% (2,1): QFI for lambda
subplot(2,2,3);
plot(t, QFI_lam, 'Color', color_lam, 'LineWidth', 2);
xlim([0 30]); set(gca,'xtick',[0 10 20 30]);
ylim([0 8]); set(gca,'ytick',[0 2 4 6 8]);
set(gca,'Fontsize',font_num,'linewidth',1,'ticklength',[0.020 0.018]);
xlabel('$t$', 'Interpreter','latex','FontSize',font_axis);
ylabel('$\mathcal{F}_{\lambda}(t)$', 'Interpreter','latex','FontSize',font_axis);
text(abcd_x, abcd_y, '(c)', 'Units','normalized','FontSize',index_abcd,...
    'FontWeight','bold','Interpreter','latex');

% (2,2): QFI for gamma
subplot(2,2,4);
plot(t, QFI_ga, 'Color', color_ga, 'LineWidth', 2);
xlim([0 30]); set(gca,'xtick',[0 10 20 30]);
ylim([0 5]); set(gca,'ytick',[0 1 2 3 4 5]);
set(gca,'Fontsize',font_num,'linewidth',1,'ticklength',[0.020 0.018]);
xlabel('$t$', 'Interpreter','latex','FontSize',font_axis);
ylabel('$\mathcal{F}_{\gamma}(t)$', 'Interpreter','latex','FontSize',font_axis);
text(abcd_x, abcd_y, '(d)', 'Units','normalized','FontSize',index_abcd,...
    'FontWeight','bold','Interpreter','latex');

% Optional: tighten layout
% set(gcf, 'Units', 'normalized');