function k = kn0(t,type,kn0_high,kn0_low,period,t_shift_fluc,t_start)
if type == 1
    % fluctuating low to high
    k = (kn0_high - kn0_low) * ((square(2*pi*t/period + pi)/2 + 0.5) - (square(2*pi*t/period + pi)/2 + 0.5) .* heaviside(t_start - t) + (square(2*pi*t/period)/2 + 0.5) .* heaviside(t - t_shift_fluc)) + kn0_low; % fluctuating, remains low until after t_start and remains high after t_shift
elseif type == 2
    % upshift
    k = (kn0_high - kn0_low) * heaviside(t - t_start) + kn0_low ;
elseif type == 3
    % downshift
    k = kn0_high - (kn0_high - kn0_low) * heaviside(t - t_start) ;
else
    % fluctuating high to low
    k = (kn0_high - kn0_low) * ((square(2*pi*t/period)/2 + 0.5) + (square(2*pi*t/period + pi)/2 + 0.5) .* heaviside(t-t_shift_fluc)) + kn0_low;
end
end