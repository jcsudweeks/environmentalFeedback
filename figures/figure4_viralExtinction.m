% figure4_viralExtinction.m
% Viral extinction dynamics in the rg1_2CO parameter set at rho = 0.65.
% Starting from E_C^S, defectors are introduced; cooperators are driven
% into the E_U basin and virus goes extinct.
%
% Requires ../data/viralExtinctionData.mat -- generate it by running
% viralExtinctionData.m first if it isn't already in the repo.

load('../data/viralExtinctionData.mat')

csColor = [0,      0.4470, 0.6980];
euColor = [0,      0.6196, 0.4510];
VC_traj = pop_d(:,5);
VD_traj = pop_d(:,6);
fC_traj = VC_traj ./ max(VC_traj + VD_traj, 1e-10);

figure('Position', [100 100 700 450]); hold on
set(gcf, 'Renderer', 'painters')
colormap(flipud(parula))
xregion(0,           VC_boundary, 'FaceColor', euColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
xregion(VC_boundary, Inf,         'FaceColor', csColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
surface([VC_traj, VC_traj], [VD_traj, VD_traj], zeros(length(VC_traj), 2), ...
    [fC_traj, fC_traj], 'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', 1.5)
xline(VC_boundary, 'k:', 'LineWidth', 3)
threshold = mu1 / mu2;
VC_line   = linspace(0, threshold, 300);
plot(VC_line, threshold - VC_line, 'k--', 'LineWidth', 3)
plot([0, 3.5e5], [0, 3.5e5], 'k','Color',[.5 .5 .5], 'LineWidth', .5)
cb = colorbar;
cb.Label.String   = 'frequency of cooperators';
cb.Label.FontSize = 20;
clim([0.4 1])
axis([0, 3.5e5, 0, 1.7e5])
xlabel('V_C', 'FontSize', 20)
ylabel('V_D', 'FontSize', 20)
ax = gca;
ax.XAxis.Exponent = 5;
ax.YAxis.Exponent = 5;
box on

