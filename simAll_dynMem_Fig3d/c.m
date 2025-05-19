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