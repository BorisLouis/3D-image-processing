function [Gr,tau] = GetAutoCorrelation(Diff,fs)
    T_win = 1; 
    N_win = round(fs*T_win);  

    trend = movmean(Diff, N_win, 'Endpoints','shrink');
    r = Diff - trend;
    r = r - mean(r);  

    %%% autocorrelation
    [acf, lags] = xcorr(r, 'unbiased');
    acf = acf(lags>=0);
    tau = lags(lags>=0)/fs;
    Gr = acf./acf(1);
end

