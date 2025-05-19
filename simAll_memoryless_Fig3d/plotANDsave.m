function plotANDsave(t, alpha, t_samp, kappa_samp, T_mem, bias, k_t0, scale, period, type)

f = figure ;
bounds = [min(t_samp),max(t_samp)] ;
% til = tiledlayout(3,1) ;
title(['T_{mem}=', num2str(T_mem), ' bias=', num2str(bias)])

nexttile([1 2])
plot(t(1:end-1),1-alpha)
ylabel('Memory (1-\alpha)')
xlim(bounds)
ylim([-0.1,max(1-alpha)+0.1])

nexttile([2 2])
plot(t_samp,kappa_samp), hold on
ylabel('Growth Rate')
xlabel('Time (h)')
xlim(bounds), hold off
% ylim([0,1.75])
% til.TileSpacing = 'tight' ;
if type == 2
    saveas(f,['fitFigures\', 'upshift', num2str(T_mem), 'bias', num2str(bias), 'kt0', num2str(k_t0), 'scale', num2str(scale), '.pdf'])
elseif type == 3
    saveas(f,['fitFigures\', 'downshift', num2str(T_mem), 'bias', num2str(bias), 'kt0', num2str(k_t0), 'scale', num2str(scale), '.pdf'])
elseif type == 4
    saveas(f,['fitFigures\', num2str(period), 'period_upstart', num2str(T_mem), 'bias', num2str(bias), 'kt0', num2str(k_t0), 'scale', num2str(scale), 'flucT.pdf'])
else
    saveas(f,['fitFigures\', num2str(period), 'period', num2str(T_mem), 'bias', num2str(bias), 'kt0', num2str(k_t0), 'scale', num2str(scale), 'flucT.pdf'])
end

end