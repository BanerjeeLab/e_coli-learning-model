% Goal: simulate memoryless model and compare to dynamic memory model and
% experimental data
% Notes:
% - uses fde_pi12_pc func from Garrappa 2017 for simulation of memoryless
% model
% - uses csv of simulation results produced by dyn_mem_model_sim.m for
% comparing to dynamic memory model


close all
clear all

%% Parameters
global phiR_min phiR_max a_n a_t k_t0 kn0_high kn0_low mu_ns period t_shift t_shift_up
phiR_min = 0.049;
phiR_max = 0.55;
a_n = 1e-3;
a_t = 1e-4;

k_t0 = 4.05;
kn0_high = 1500;
kn0_low = 2.4;
mu_ns = 0;

period = 1;
t0 = 0; T = 8;
t_shift = 7 ;
t_shift_up = 1 ;


%% Defining functions
f = @(a) 1 ./ (1 + (a/a_n).^2);
f_prime = @(a) -2*a/a_n^2 ./ (1 + (a/a_n).^2).^2;
g = @(a) (a/a_t).^2 ./ (1 + (a/a_t).^2);
g_prime = @(a) 2*a/a_t^2 ./ (1 + (a/a_t).^2).^2;
f_R = @(a) (-f_prime(a).*g(a)*phiR_max + f(a).*g_prime(a)*phiR_min) ./ (-f_prime(a).*g(a) + f(a).*g_prime(a));

k_n0 = @(t) (kn0_high - kn0_low) * ((square(2*pi*t/period + pi)/2 + 0.5) + (square(2*pi*t/period)/2 + 0.5) .* heaviside(t-t_shift)) + kn0_low; % fluctuating, remains high after t_shift

%% Solving for initial conditions
phiR_ss = @(a) (k_n0(0)*f(a)*phiR_max + k_t0*g(a)*phiR_min + mu_ns) / (k_n0(0)*f(a) + k_t0*g(a));

func = @(x) phiR_ss(x) - f_R(x);
a0 = fzero(func,1e-4);
phi_R0 = f_R(a0);
y0 = [phi_R0; a0];

%% Defining first-order ODEs
% phi_R = y(1)  a = y(2)
dphiR_dt = @(t, y) k_t0 * g(y(2)) * (y(1) - phiR_min) * (f_R(y(2)) - y(1));
dAAdt = @(t, y) k_n0(t) * f(y(2)) * (phiR_max - y(1)) - k_t0 * g(y(2)) * (y(1) - phiR_min) + mu_ns;

f_fun = @(t, y) [dphiR_dt(t, y); dAAdt(t, y)];

kappa = @(y) k_t0 .* g(y(2,:)) .* (y(1,:) - phiR_min) - mu_ns;


%% environment uncertainty estimate - deterministic

global t_start t_end
h = 2^(-14) ;
T_mem = 1.9 ; % the length of memory used to estimate environment, in hrs
N_mem = ceil(T_mem / h) ;
T_delay = 0 ; % the length in delay between estimate and action, in hrs
N_delay = ceil(T_delay / h) ;
n_points = ceil((T - (t0-T_mem-T_delay)) / h) ;

t_start = 0 ; t_end = 7 ;
t = linspace(t0-T_mem-T_delay,T,n_points) ;

% c = @(t) heaviside(t - t_shift) ; % upshift
c = @(t) (square(2*pi*t/period + pi)/2 + 0.5) - (square(2*pi*t/period + pi)/2 + 0.5) .* heaviside(t_start - t) + (square(2*pi*t/period)/2 + 0.5) .* heaviside(t - t_end); % fluctuating, remains low until after t_shift


% calculating uncertainty
uncertainty = zeros(1,length(N_mem+N_delay+1:n_points)) ;

j = 1 ;
for i = N_mem+N_delay+1:n_points
    
    t_mem = t(i-N_mem-N_delay:i-N_delay-1) ;
    k_n = sum(c(t_mem)) ;
    pi_n = k_n / N_mem ;
    
    uncertainty(j) = (0.5 - abs(pi_n - 0.5)) / 0.5 ;
    j = j+1;
end

t_sub = t(N_mem+N_delay+1:n_points) ;

scale = 0.7 ;
alpha = 1 - scale * uncertainty ; alpha = alpha(t_sub >= t0) ;


%% Plot for paper (fluctuating)

% simulate memoryless model
h = 2^(-14) ;

alpha_nomem = [1; 1] ;
[t, y] = fde_pi12_pc(alpha_nomem, f_fun, t0, T, y0, h); % for simulating MOS


% import experimental data and var-order sim results
data = readtable('pulsing_nguyen.csv') ;
dat = readtable('pulsing_dyn_mem_model.csv') ; % importing results from previous simulation
t_var = dat.Var1 ;
kappa_var = dat.Var2 ;
phiR_var = dat.Var3 ;
a_var = dat.Var4 ;
%%

figure(1)
bounds = [0,8] ;
til = tiledlayout(14,1,'TileSpacing','none') ;

nexttile([2 1])
plot(t_sub,alpha*0, '-.','Color','black','LineWidth',2), hold on
plot(t_sub,1-alpha,'Color','black','LineWidth',2)
ylabel('Memory, \it{k}', 'FontSize',18)
xlim(bounds)
ylim([-0.1,max(1-alpha)+0.1])
xticks([])

nexttile([3 1])
plot(t,y(1,:),'-.','Color','black','LineWidth',2), hold on
plot(t_var,phiR_var,'Color','black','LineWidth',2)
ylabel('\phi_R', 'FontSize',18), hold off
xlim(bounds)
xticks([])

nexttile([3 1])
semilogy(t,y(2,:),'-.','Color','black','LineWidth',2), hold on
semilogy(t_var,a_var,'Color','black','LineWidth',2)
ylabel('\it{a}', 'FontSize',19), hold off
xlim(bounds)
xticks([])

nexttile([6 1])
ax = gca;
scatter(data.time-2.955,data.kappa,[],'white'), hold on
for i = 1:6
    rectangle('Position', [i-0.5, ax.YLim(1), 0.5, diff(ax.YLim)], 'FaceColor', [0.95, 0.95, 0.95], 'EdgeColor', 'none');
end
rectangle('Position', [6.5, ax.YLim(1), 1.5, diff(ax.YLim)], 'FaceColor', [0.95, 0.95, 0.95], 'EdgeColor', 'none');
scatter(data.time-2.955,data.kappa,100,'filled',"MarkerFaceColor",[0.4940 0.1840 0.5560]) ;
plot(t,kappa(y), '-.','Color','black','LineWidth',2)
plot(t_var,kappa_var,'Color','black','LineWidth',2)
ylabel('Growth Rate (h^{-1})', 'FontSize',18)
xlabel('Time (h)', 'FontSize',18)
xlim(bounds)
legend('','Nguyen et al., 2021','memoryless','variable-order memory', 'FontSize',14, 'Location','northwest')
hold off
