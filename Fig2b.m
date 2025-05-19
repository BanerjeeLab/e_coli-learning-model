% Goal: simulate memoryless model and compare to data to show it cannot
% explain observed growth control
% Note: uses fde_pi12_pc.m func from Garrappa 2017


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

k_n0 = @(t) (kn0_high - kn0_low) * heaviside(t - t_shift_up) + kn0_low; % upshift

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

%% Plot for paper (memoryless upshift)

% simulate
h = 2^(-12) ;

alpha = [1; 1] ;
[t, y] = fde_pi12_pc(alpha, f_fun, t0, T, y0, h); % for simulating MOS


% import experimental data from Nguyen 2021
data_up = readtable('upshift_nguyen.csv') ;
time_up = data_up.Var1/60+1 ;
kappa_up = data_up.Var2*log(2) ;

figure(1)
bounds = [0,5] ;

ax = gca;
scatter(time_up,kappa_up,[],'white'), hold on
rectangle('Position', [1, ax.YLim(1), 4, diff(ax.YLim)], 'FaceColor', [0.95, 0.95, 0.95], 'EdgeColor', 'none');
scatter(time_up,kappa_up,100,'filled',"MarkerFaceColor",[0.4940 0.1840 0.5560]) ;
plot(t,kappa(y), '--','Color','black','LineWidth',2)
ylabel('Growth Rate (h^{-1})', 'FontSize',22)
xlabel('Time (h)', 'FontSize',22)
xlim(bounds)
legend('','Nguyen et al., 2021','memoryless (\it{k}=0)', 'FontSize',18, 'Location','southeast')
hold off

% save figure
% saveas(gcf,'figures/memoryless_upshift_figure.pdf')
