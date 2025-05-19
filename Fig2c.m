% Goal: simulate constant order model to show growth control properties


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

period = 3;
t0 = 0; T = 4.5;
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


%% Plot for paper (fluctuating)

% simulate
h = 2^(-14) ;
al_range = [1, 0.7, 0.4] ;


figure(1)
set(gca, 'FontSize', 16);
bounds = [0,4.45] ;

% ax = gca;
rectangle('Position', [1.5, 0, 1.5, 2], 'FaceColor', [0.95, 0.95, 0.95], 'EdgeColor', 'none'); hold on

colors = [
    0, 0.5, 0.5;
    0.1, 0.6, 0.5;
    0.2, 0.7, 0.6;
];
linestyl = {':', '-.', '-'};

for i =1:length(al_range)
    alpha = [al_range(i); 1] ;
    [t, y] = fde_pi12_pc(alpha, f_fun, t0, T, y0, h); % for simulating MOS
    plot(t,kappa(y),'LineWidth',2.6,'DisplayName',num2str(1-al_range(i)),'Color',colors(i,:),'LineStyle',linestyl(i))
end

ylabel('Growth Rate (h^{-1})', 'FontSize',22)
xlabel('Time (h)', 'FontSize',22)
xlim(bounds)
lgd = legend('Location','southwest');
lgd.FontSize = 14 ;
lgd.Title.FontSize = 15 ;
lgd.Title.String = 'Memory strength, \it{k}'; hold off

% save figure
% saveas(gcf,'figures/constantorder_pulse_figure.pdf')

