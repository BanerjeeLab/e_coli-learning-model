% Goal: plotting experimental and simulation results nicely for paper for
% Figure 3a and 3c


close all
clear all


%% experimental data details
per = [0 2] ;
file = ["upshift_huijing.csv", "1h_pulse_huijing.csv"];
offset = [7+5.8 5] ;
start = [0 1.15] ;
exp_type = [2 1] ; % type of experiment (upshift, downshift, fluctuating)

bound_up = [0 7] ;

n_periods = [0; 6] ;


% T_mem = x(1) ; bias = x(2) ; k_t0 = x(3) ; scale = x(4) ;
x = [1.3 0.49913 1.85 0.3967] ;

h = 2^(-10) ; T = 15 ; t_shift_fluc = 12 ;

% choose i=1 for Fig 3a, i=2 for Fig 3c
i = 2 ;
period = per(i) ;
t_start = start(i) ;
type = exp_type(i) ;
b = bound_up(i) ;

if type == 1
    t0 = 0 ;
    X = readtable(['T_vs_tau_period', num2str(period), 'bias', num2str(x(2)), 'mem', num2str(x(1)), 'kt0', num2str(x(3)), '.csv']) ;
elseif type == 2
    t0 = -3 ;
    X = readtable(['T_vs_tau_upshift_bias', num2str(x(2)), 'mem', num2str(x(1)), 'kt0', num2str(x(3)), '.csv']) ;
end
t = X.Var1 ; z = X.Var2 ;

% importing experimental data
data = readtable(file(i)) ;
time_exp = data.Var1-offset(i) ;
kappa_exp = data.Var2 ;


% environment uncertainty estimate

N_mem = ceil(x(1) / h) ;
bias_N = x(2)*N_mem ;
n_points = ceil((T - (t0-x(1))) / h) ;
tc = linspace(t0-x(1),T,n_points) ;

conc = @(t,type) c(t,type,period,t_shift_fluc,t_start) ;

% calculating uncertainty
uncertainty = zeros(1,length(N_mem+1:n_points)) ;
j = 1 ;
for i = N_mem+1:n_points
    
    t_mem = tc(i-N_mem:i-1) ;
    k_n = sum(conc(t_mem,type)) ;
    pi_n = (k_n + bias_N) / (N_mem + bias_N) ;
    
    uncertainty(j) = (1 - abs(2*pi_n - 1)) + 0.01 ;
    j = j+1;
end
t_sub = tc(N_mem+1:n_points) ;
alpha = 1 - x(4) * uncertainty ; alpha = alpha(t_sub >= t0) ;



fig = figure ;
til = tiledlayout(3,1,'TileSpacing','none') ;

nexttile
if type == 2
    xbounds = [-3,6] ;
    ybounds = [0.1,1] ;
elseif type == 1
    xbounds = [0,12] ;
    ybounds = [0,0.9] ;
end
plot(t_sub(t_sub >= t0),1-alpha,'Color','black','LineWidth',2), hold on
ylabel('Memory, \it{k}', 'FontSize',18)
xlim(xbounds)
ylim([-0.1,max(1-alpha)+0.1])
xticks([])

nexttile([2 1])
ax = gca;
scatter(time_exp,kappa_exp,[],'white'), hold on
if type == 2
    rectangle('Position', [0, 0, 10, 1], 'FaceColor', [0.95, 0.95, 0.95], 'EdgeColor', 'none');    
    scatter(time_exp,kappa_exp,100,'filled',"MarkerFaceColor",[0.4940 0.1840 0.5560]) ;
    plot(t-6,z,'Color','black','LineWidth',2.5)
elseif type == 1
    for i = 0:5
        rectangle('Position', [2*i+1, ax.YLim(1), 1, diff(ax.YLim)], 'FaceColor', [0.95, 0.95, 0.95], 'EdgeColor', 'none');
    end
    rectangle('Position', [13, ax.YLim(1), 3, diff(ax.YLim)], 'FaceColor', [0.95, 0.95, 0.95], 'EdgeColor', 'none');
    scatter(time_exp,kappa_exp,100,'filled',"MarkerFaceColor",[0.4940 0.1840 0.5560]) ;
    plot(t,z,'Color','black','LineWidth',2.5)
end
ylabel('Growth Rate (h^{-1})', 'FontSize',18)
xlabel('Time (h)', 'FontSize',18)
ylim(ybounds)
xlim(xbounds)
hold off

%%

function k = c(t,type,period,t_shift_fluc,t_start)
    if type == 1
        % fluctuating low to high
        k = ((square(2*pi*t/period + pi)/2 + 0.5) - (square(2*pi*t/period + pi)/2 + 0.5) .* heaviside(t_start - t) + (square(2*pi*t/period)/2 + 0.5) .* heaviside(t - t_shift_fluc)) ; % fluctuating, remains low until after t_start and remains high after t_shift
    elseif type == 2
        % upshift
        k = heaviside(t - t_start) ;
    elseif type == 3
        % downshift
        k = 1 - heaviside(t - t_start) ;
    else
        % fluctuating high to low
        k = 1 - (square(2*pi*t/period + pi)/2 + 0.5) + (square(2*pi*t/period + pi)/2 + 0.5) .* heaviside(t_start - t) - (square(2*pi*t/period)/2 + 0.5) .* heaviside(t - t_shift_fluc) ; % fluctuating, remains high until after t_start and remains low after t_shift
    end
end