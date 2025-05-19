function [kappa_samp, t_samp] = fdeint_fit(x, kappa, ffun, t0,T,h,type,c,g,f_R,kappa_high,kappa_low,zero_ss,period)

T_mem = x(1) ; bias = x(2) ; k_t0 = x(3) ; scale = x(4) ; mem_min = 0.01 ;

% solving for initial conditions
fun = @(x) ss_func(x,k_t0, zero_ss,g,f_R,kappa_high,kappa_low) ;
options = optimoptions('fsolve','MaxIterations',1000) ;
x_solve = fsolve(fun,[1e-4, 2, 1e-3, 10],options);
if type < 3
    a0 = x_solve(1) ;
else
    a0 = x_solve(3) ;
end
phi_R0 = f_R(a0) ; kn0_low = x_solve(2) ; kn0_high = x_solve(4) ;
y0 = [phi_R0; a0];

f_fun = @(t, y) ffun(t, y, type, k_t0, kn0_high, kn0_low) ;


% environment uncertainty estimate

N_mem = ceil(T_mem / h) ;
bias_N = bias*N_mem ;
n_points = ceil((T - (t0-T_mem)) / h) ;
t = linspace(t0-T_mem,T,n_points) ;

% calculating uncertainty
uncertainty = zeros(1,length(N_mem+1:n_points)) ;
j = 1 ;
for i = N_mem+1:n_points
    
    t_mem = t(i-N_mem:i-1) ;
    k_n = sum(c(t_mem,type)) ;
    pi_n = (k_n + bias_N) / (N_mem + bias_N) ;
    
    uncertainty(j) = (1 - abs(2*pi_n - 1)) + mem_min ;
    j = j+1;
end

t_sub = t(N_mem+1:n_points) ;
alpha = 1 - scale * uncertainty ; alpha = alpha(t_sub >= t0) ;



% pred_corr parameters
mu_tol = 1.0e-6 ;
mu = 0 ;
param = [] ;

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

N = ceil((T-t0)/h) ;  % number of timesteps being simulated
if N > length(alpha)
    alpha(end+1) = alpha(end) ; % making sure alpha is correct length
end

% Storage of information on the problem
Probl.ic = ic ;
Probl.f_fun = f_fun ;
Probl.problem_size = problem_size ;
Probl.param = param ;
Probl.alpha = alpha ;
Probl.y0 = y0 ;

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
[y, ~] = fdeint(1, N, t, y, fy, METH, Probl ) ;

% Evaluation solution in T when T is not in the mesh
if T < t(N+1)
    c = (T - t(N))/h ;
    t(N+1) = T ;
    y(:,N+1) = (1-c)*y(:,N) + c*y(:,N+1) ;
end

t = t(1:N+1) ; y = y(:,1:N+1) ;


% Returning only y values corresponding to input time values

% if type == 1
% downsampling to ease computation cost
n_samp = 100 ;
t_samp = t(1:n_samp:end)' ;
y_samp = y(:,1:n_samp:end) ;
kappa_samp = kappa(y_samp, k_t0)' ;
% else
%     t_samp = t' ;
%     kappa_samp = kappa(y, k_t0)' ;
% end

% Will use later if fitting
% % finding corresponding y value for each input t value
% A = repmat(t_samp, [1 length(time)]) ;
% [~,closestIndex] = min(abs(A-time')) ;
% y = kappa_samp(closestIndex) ;


% Saving plot of result
plotANDsave(t, alpha, t_samp, kappa_samp, T_mem, bias, k_t0, scale, period, type)

end