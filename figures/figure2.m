% figure2.m
% Phase plane of the mixed limit cycle in (V_C, V_D) space for the
% rg1_2CO parameter set, with the single/double infection threshold
% V_C + V_D = mu1/mu2 marked.
%
% Requires ../data/mixedLCPhasePlaneData.mat -- generate it by running
% mixedLCPhasePlaneData.m first if it isn't already in the repo.

load('../data/mixedLCPhasePlaneData.mat', 'VCdata', 'VDdata', 'fCdata', 'rhoVec', 'threshold')


%cmap = flipud(gray);   % white (f_C=0, defectors) -> black (f_C=1, cooperators)
cmap = flipud(parula);
figure('Position', [100 100 700 450]); hold on
set(gcf, 'Renderer', 'painters')
colormap(cmap)

for k = 1:length(rhoVec)
    scatter(VCdata{k}, VDdata{k}, 3, fCdata{k}, 'filled', 'HandleVisibility', 'off')
end

% Threshold line: V_C + V_D = mu1/mu2
VC_line = linspace(0, threshold, 300);
VD_line = threshold - VC_line;
plot(VC_line, VD_line, 'k--', 'LineWidth', 1, ...
    'DisplayName', '$V_C + V_D = \mu_1/\mu_2$')

% Equal density line: V_D = V_C
plot([0, threshold], [0, threshold], 'Color',[.5 .5 .5], 'LineWidth', .5, ...
    'DisplayName', '$V_D = V_C$')
ylim([0 1.7*(10^5)])
ax = gca;
ax.YAxis.Exponent = 5;
xlabel('V_C',  'FontSize', 20)
ylabel('V_D',  'FontSize', 20)
%title('Mixed limit cycle: (V_C, V_D) phase plane', 'FontSize', 20)
cb = colorbar;
cb.Label.String = 'frequency of cooperators';
cb.Label.FontSize = 20
clim([0.5 1])

box on


