function [vR] = GetRotationalSpeed(Gr, tau)

    [peaks,locs, w, p] = findpeaks(medfilt1(Gr, 100), tau,'MinPeakProminence',0.1);
    if length(locs) >= 2
        period = median(diff(locs));
        vR = 360/(period*4); % grad/s
    else
        vR = NaN; % no clear oscillation
    end

end

