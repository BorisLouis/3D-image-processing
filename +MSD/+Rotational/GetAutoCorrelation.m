function [Gr,tau] = GetAutoCorrelation(Diff,fs, ExpTime)
    T_win = 1; 
    N_win = round(fs*T_win);  

    trend = movmean(Diff, N_win, 'Endpoints','shrink');
    r = Diff - trend;
    r = r - mean(r);  

    %%% autocorrelation
    [acf, lags] = xcorr(r, 'unbiased');
    acf = acf(lags>=0);
    tau = [ExpTime:ExpTime:size(acf,1)*ExpTime];
    Gr = acf./acf(1);
end

