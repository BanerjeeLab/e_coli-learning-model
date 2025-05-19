function F = ss_func(x, k_t0, zero_ss,g,f_R,kappa_high,kappa_low)
    % x(1) = a_low, x(2) = kn0_low, x(3) = a_high, x(4) = kn0_high
    F(1) = zero_ss(x(1),x(2),k_t0) ;
    F(2) = zero_ss(x(3),x(4),k_t0) ;
    F(3) = k_t0*g(x(3))*(f_R(x(3))-0.049) - kappa_high ;
    F(4) = k_t0*g(x(1))*(f_R(x(1))-0.049) - kappa_low ;
end