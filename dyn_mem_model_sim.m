% Goal: Simulate dynamic memory model in fluctuating environments.
% Simulation results are used in Fig2d.m
% Note: simulation uses the predictor-corrector approach described in
% Garrappa 2017. The function "fdeint" is a heavily modified version of
% fde_pi12_pc.m (from Garrappa 2017) which allows for simulation of
% variable-order FDEs (fde_pi12_pc.m only handles constant-order FDEs).
% Warning: due to the nature of the convolution across time, the algorithm
% used in this code is very slow (O(n^2)).

close all
clear all

%% Parameters
global phiR_min phiR_max a_n a_t k_t0 kn0_high kn0_low mu_ns period t_shift
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


%% pred_corr parameters
mu_tol = 1.0e-6 ;
mu = 1 ;
param = [] ;

%% Dynamic memory model simulation. This code 

% Check compatibility size of the problem with number of fractional orders
problem_size = size(y0,1) ;
alpha_start = [1, 1] ;

% Storage of initial conditions
ic.t0 = t0 ;
ic.y0 = y0 ;
ic.m_alpha = ceil(alpha_start) ;
for i = 1 : problem_size
    for j = 0 : ic.m_alpha(i)-1
        ic.m_alpha_factorial(i,j+1) = factorial(j) ;
    end
end

% Storage of information on the problem
Probl.ic = ic ;
Probl.f_fun = f_fun ;
Probl.problem_size = problem_size ;
Probl.param = param ;
Probl.alpha = alpha ;
Probl.y0 = y0 ;

%%
N = ceil((T-t0)/h) ; % number of timesteps being simulated
f_temp = f_vectorfield(t0,y0(:,1),Probl) ;

% Preallocation of some variables
y = zeros(Probl.problem_size,N+1) ;
fy = zeros(Probl.problem_size,N+1) ;

% Defining coefficient function
bn = @(n_t,n) (n+1).^alpha(n_t) - n.^alpha(n_t) ;
an = @(n_t,n) (n-1).^(alpha(n_t)+1) - 2*n.^(alpha(n_t)+1) + (n+1).^(alpha(n_t)+1) ;
a0 = @(n_t,n) (n-1).^(alpha(n_t)+1) - (n-alpha(n_t)-1).*n.^alpha(n_t) ;

METH.mu = mu ; METH.mu_tol = mu_tol ;
METH.bn = bn ; METH.h = h ;
METH.an = an ; METH.a0 = a0 ;

% Initializing solution and process of computation
t = t0 + (0 : N)*h ;
y(:,1) = y0(:,1) ;
fy(:,1) = f_temp ;
[y, fy] = fdeint(1, N, t, y, fy, METH, Probl ) ;

% Evaluation solution in T when T is not in the mesh
if T < t(N+1)
    c = (T - t(N))/h ;
    t(N+1) = T ;
    y(:,N+1) = (1-c)*y(:,N) + c*y(:,N+1) ;
end

t = t(1:N+1) ; y = y(:,1:N+1) ;

%% Export data

% export simulation results
n = 150 ; % sample spacing

phiR_samp = y(1,1:n:end) ;
a_samp = y(2,1:n:end) ;
kappa_samp = kappa([phiR_samp; a_samp]) ;
alpha_samp = alpha(1:n:end) ;
t_samp = t(1:n:end) ;

dat = [t_samp; kappa_samp; phiR_samp; a_samp]' ;
writematrix(dat, 'pulsing_dyn_mem_model.csv')

%%

% =========================================================================
% =========================================================================
function [y, fy] = fdeint(nxi, nxf, t, y, fy, METH, Probl)

for n = nxi : nxf
    
    n_t = n ;
    
    % Evaluation of the predictor
    y_pred = zeros(Probl.problem_size,1) ;
    Phi = 0 ;
    j_beg = 0 ;
    for j = j_beg : n-1
        Phi = Phi + METH.bn(n_t,n-j).*fy(1,j+1) ;
    end
    St = StartingTerm(t(n+1),Probl.ic) ;
    y_pred(1) = St(1) + Phi./gamma(Probl.alpha(n_t)+1) .* METH.h.^Probl.alpha(n_t) ;
    y_pred(2) = y(2,n) + METH.h.*fy(2,n) ;
    f_pred = f_vectorfield(t(n+1),y_pred,Probl) ;
    
    % Evaluation of the corrector
    if METH.mu == 0
        y(:,n+1) = y_pred ;
        fy(:,n+1) = f_pred ;
    else
        j_beg = nxi ;
        Phi = zeros(Probl.problem_size,1) ;
        for j = j_beg : n-1
            if n-j == 0
                Phi = Phi + fy(:,j+1) ;
            else
                Phi(1) = Phi(1) + METH.an(n_t,n-j+1).*fy(1,j+1) ;
                Phi(2) = Phi(2) + 2*fy(2,j+1) ;
            end
        end
        Phi_n = zeros(Probl.problem_size,1) ;
        Phi_n(1) = St(1) + ...
            (METH.a0(n_t,n+1).*fy(1,1) + Phi(1)).*METH.h.^Probl.alpha(n_t)./gamma(Probl.alpha(n_t)+2) ;
        Phi_n(2) = St(2) + ...
            (METH.a0(1,n+1).*fy(2,1) + Phi(2)).*METH.h./gamma(1+2) ;
        yn0 = y_pred ; fn0 = f_pred ;        
        stop = 0 ; mu_it = 0 ;
        while ~stop
            yn1 = zeros(Probl.problem_size,1) ;
            yn1(1) = Phi_n(1) + METH.h.^Probl.alpha(n_t)/gamma(Probl.alpha(n_t)+2).*fn0(1) ;
            yn1(2) = Phi_n(2) + METH.h./gamma(1+2).*fn0(2) ;
            mu_it = mu_it + 1 ;
            if METH.mu == Inf
                stop = norm(yn1-yn0,inf) < METH.mu_tol ;
                if mu_it > 100 && ~stop
                    warning('MATLAB:fde12:NonConvegence',...
                        strcat('It has been requested to run corrector iterations until convergence but ', ...
                        'the process does not converge to the tolerance %e in 100 iteration'),METH.mu_tol) ;
                    stop = 1 ;
                end
            else
                stop = mu_it == METH.mu ;
            end
            fn1 = f_vectorfield(t(n+1),yn1,Probl) ;
            yn0 = yn1 ; fn0 = fn1 ;
        end
        y(:,n+1) = yn1 ;
        fy(:,n+1) = fn1 ;
    end
end

end


% =========================================================================
% =========================================================================
function f = f_vectorfield(t,y,Probl)

if isempty(Probl.param)
    f = feval(Probl.f_fun,t,y) ;
else
    f = feval(Probl.f_fun,t,y,Probl.param) ;
end

end


% =========================================================================
% =========================================================================
function ys = StartingTerm(t,ic)

ys = zeros(size(ic.y0,1),1) ;
for k = 1 : max(ic.m_alpha)
    if length(ic.m_alpha) == 1
        ys = ys + (t-ic.t0)^(k-1)/ic.m_alpha_factorial(k)*ic.y0(:,k) ;
    else
        i_alpha = find(k<=ic.m_alpha) ;
        ys(i_alpha,1) = ys(i_alpha,1) + (t-ic.t0)^(k-1)*ic.y0(i_alpha,k)./ic.m_alpha_factorial(i_alpha,k) ;
    end
end

end