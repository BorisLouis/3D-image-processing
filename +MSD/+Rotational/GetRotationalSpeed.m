function [vR] = GetRotationalSpeed(Gr, tau)

    [peaks,locs] = findpeaks(Gr, tau);
    if length(locs) >= 2
        period = mean(diff(locs));
        vR = 2*pi/period; % rad/s
    else
        vR = NaN; % no clear oscillation
    end

end

