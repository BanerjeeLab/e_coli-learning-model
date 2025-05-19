 % Goal: simulate dynamic memory model for different nutrient pulsing periods
 % and measure adaptation time tau. Simulates different periods in parallel


close all
clear all


%% experimental data details
per = [0.5 1 2 2.5 3 3.5 4 5 6] ;
exp_type = [1 1 1 1 1 1 1 1 1] ; % type of experiment (upshift, downshift, fluctuating)
start = [0 0 1.15 0 2 0 3 0 3.8] ;
bound_up = [6.75 6.5 7 6.25 7.5 8.75 6 7.5 9] ;


% setting parameters
t0 = 0 ; T = 12 ; h = 2^(-13) ;

% T_mem = x(1) ; bias = x(2) ; k_t0 = x(3) ; scale = x(4) ;
x = [1.3 0.49913 1.85 0.3967] ;


kap_cell = cell(1,length(per));
t_cell = cell(1,length(per));
tau_up = zeros(1,length(per)) ;
tau_down = zeros(1,length(per)) ;


if max(size(gcp)) == 0 % parallel pool needed
    parpool % create the parallel pool
end

parfor i=1:length(per)
    period = per(i) ;
    t_start = start(i) ; % 0
    type = exp_type(i) ;
    b = bound_up(i) ;
    

    % global phiR_min phiR_max a_n a_t mu_ns t_shift_up t_shift_fluc t_start t_end kappa_high kappa_low
    phiR_min = 0.049;
    phiR_max = 0.55;
    a_n = 1e-3;
    a_t = 1e-4;

    kappa_high = 0.81 ; % experimentally measured values
    kappa_low = 0.29 ;

    t_shift_fluc = 100 ; % fluctuating, change when want to stay up
    

    % environment estimate parameters
    t_end = 100 ; % change when want to stay up
    
    % Defining functions
    f = @(a) 1 ./ (1 + (a/a_n).^2);
    f_prime = @(a) -2*a/a_n^2 ./ (1 + (a/a_n).^2).^2;
    g = @(a) (a/a_t).^2 ./ (1 + (a/a_t).^2);
    g_prime = @(a) 2*a/a_t^2 ./ (1 + (a/a_t).^2).^2;
    f_R = @(a) (-f_prime(a).*g(a)*phiR_max + f(a).*g_prime(a)*phiR_min) ./ (-f_prime(a).*g(a) + f(a).*g_prime(a));

    
    % concentration functions
    k_n0 = @(t,type,kn0_high,kn0_low) kn0(t,type,kn0_high,kn0_low,period,t_shift_fluc,t_start) ;

    conc = @(t,type) c(t,type,period,t_shift_fluc,t_start) ;

    % Defining first-order ODEs
    % phi_R = y(1)  a = y(2)
    dphiR_dt = @(t, y, k_t0) k_t0 * g(y(2)) * (y(1) - phiR_min) * (f_R(y(2)) - y(1));
    dAAdt = @(t, y, type, k_t0, kn0_high, kn0_low) k_n0(t,type,kn0_high,kn0_low) * f(y(2)) * (phiR_max - y(1)) - k_t0 * g(y(2)) * (y(1) - phiR_min);

    ffun = @(t, y, type, k_t0, kn0_high, kn0_low) [dphiR_dt(t, y, k_t0); dAAdt(t, y, type, k_t0, kn0_high,kn0_low)];

    kappa = @(y, k_t0) k_t0 .* g(y(2,:)) .* (y(1,:) - phiR_min);

    zero_ss = @(a, kn0, k_t0) (kn0*f(a)*phiR_max + k_t0*g(a)*phiR_min) / (kn0*f(a) + k_t0*g(a)) - f_R(a) ; % function for solving for initial conditions


    [z, t] = fdeint_fit(x, kappa, ffun, t0,T,h,type,conc,g,f_R,kappa_high,kappa_low,zero_ss,period) ;
    
    if type == 2
        writematrix([t, z], ['T_vs_tau_upshift_bias', num2str(x(2)), 'mem', num2str(x(1)), 'kt0', num2str(x(3)), '.csv'])
    elseif type == 3
        writematrix([t, z], ['T_vs_tau_downshift_bias', num2str(x(2)), 'mem', num2str(x(1)), 'kt0', num2str(x(3)), '.csv'])
    elseif type == 4
        writematrix([t, z], ['T_vs_tau_period_upstart', num2str(period), 'bias', num2str(x(2)), 'mem', num2str(x(1)), 'kt0', num2str(x(3)), '.csv'])
    else
        writematrix([t, z], ['T_vs_tau_period', num2str(period), 'bias', num2str(x(2)), 'mem', num2str(x(1)), 'kt0', num2str(x(3)), '.csv'])
    end

    % determing relaxation timescale for fluctuating environments
    if (type == 1) || (type == 4)
        kap_cell{i} = z;
        t_cell{i} = t;

        fig = figure ;
        plot(t,z), hold on

        % finding relaxation time
        % x(1)=A x(2)=tau x(3)=b
        exp_model = @(x,t) x(1)*exp(-t/x(2)) + x(3) ;

        % UPSHIFT
        % identifying corresponding simulation trajectory
        if type == 4
            index = (t_cell{i} > b) & (t_cell{i} < b+period/2.5) ;
        else
            index = (t_cell{i} > b) & (t_cell{i} < b+period/2) ;
        end
        t_real = t_cell{i}(index) ;
        t_fit = t_real - min(t_real) ;
        g_fit = kap_cell{i}(index) ; 

        plot(t_real,g_fit)

        % exponential fit
        x0 = [-1,40,2] ;
        fit_params = lsqcurvefit(exp_model,x0,t_fit,g_fit) ;
        tau_up(i) = fit_params(2) ;

        plot(t_real,exp_model(fit_params,t_fit),'linewidth',2)

    end
end

% importing experimental data
data = readtable('TvTau_upHuijing.csv') ;
T_exp = data.Var1 ;
tau_exp = data.Var2 ;

f = figure ;
plot(per, tau_up, 'DisplayName','Upshift','linewidth',2,'color','black'), hold on
scatter(T_exp, tau_exp,100)
xlabel('period, T (h)','FontSize',18)
ylabel('adaptation time, \tau (h)','FontSize',18)
ylim([min(tau_exp) max(tau_exp)+1]), hold off
% saveas(f,['T_vs_tau_up_mu', num2str(x(2)), 'mem', num2str(x(1)), 'kt0', num2str(x(3)), '.pdf'])
