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