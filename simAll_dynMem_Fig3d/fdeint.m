function [y, fy] = fdeint(nxi, nxf, t, y, fy, METH, Probl)

for n = nxi : nxf
    
    n_t = n ;
    
    % Evaluation of the predictor
    y_pred = zeros(Probl.problem_size,1) ;
    Phi = 0 ;
    j_beg = 0 ;
    for j = j_beg : n-1
        Phi = Phi + METH.bn(j+1,n-j).*fy(1,j+1)./gamma(Probl.alpha(j+1)+1) .* METH.h.^Probl.alpha(j+1) ;
    end
    St = StartingTerm(t(n+1),Probl.ic) ;
    y_pred(1) = St(1) + Phi ;
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
                Phi(1) = Phi(1) + METH.an(j+1,n-j+1).*fy(1,j+1).*METH.h.^Probl.alpha(j+1)./gamma(Probl.alpha(j+1)+2);
                Phi(2) = Phi(2) + 2*fy(2,j+1) ;
            end
        end
        Phi_n = zeros(Probl.problem_size,1) ;
        Phi_n(1) = St(1) + ...
            (METH.a0(n_t,n+1).*fy(1,1) + Phi(1)) ; % not sure if this n_t should be changed to j+1
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