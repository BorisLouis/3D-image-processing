function [R,lags] = GetAutoCorrelation(Diff,fs, ExpTime)
   
    t = [ExpTime:ExpTime:size(Diff, 1)*ExpTime]';
    ToRemove = isnan(Diff);
    Diff(ToRemove) = [];
    t(ToRemove) = [];
    x = Diff(:);
    N = length(x);

    lags = 0:ExpTime:max(t);
    R = zeros(size(lags));

    for k = 1:length(lags)
        tau = lags(k);
        num = 0;
        den = 0;

        % For each time point, find matches at t+tau
        for i = 1:N
            % Find index where t_j ≈ t_i + tau
            [~, j] = min(abs(t - (t(i) + tau)));
            if abs(t(j) - (t(i) + tau)) < ExpTime/2
                num = num + x(i)*x(j);
                den = den + 1;
            end
        end

        if den > 0
            R(k) = num / den;
        else
            R(k) = NaN; % no pairs found
        end
    end

    % Normalize
    if ~isempty(R)
        R = R / R(1);
    else
        R = NaN;
        lags = NaN;
    end


    % tau = [ExpTime:ExpTime:size(Diff,1)*ExpTime]';
    % T_win = 1; 
    % N_win = round(fs*T_win);  
    % 
    % trend = movmean(Diff, N_win, 'Endpoints','shrink');
    % r = Diff - trend;
    % r = r - mean(r);  
    % 
    % %%% autocorrelation
    % [acf, lags] = xcorr(tau, r, 'unbiased');
    % acf = acf(lags>=0);
    % tau = [ExpTime:ExpTime:size(acf,1)*ExpTime];
    % Gr = acf./acf(1);
end

