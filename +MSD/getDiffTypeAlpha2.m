% This function aim at determining which type of diffusion we have by
% determining the anomaly diffusion coefficient alpha. According to the
% following source : https://tinevez.github.io/msdanalyzer/tutorial/MSDTuto_confined.html#4
% This can be done by fitting msd curve in a log-log scale and obtaining
% the slope as alpha.
%
% if alpha >1 : super-diffusive regime (e.g. active transport, drift)
% if alpha = 1: brownian motion
% if alpha < 1: sub-diffusive to confine
%
% INPUT:
%   - msd: msd curve for the data
%   - expTime: exposure time, time between two frame
% OUPUT:
%   - alpha: abnormal diffusion coefficient

function [alpha] = detDiffTypeAlpha2(msd,expTime,stepsize)

    if 0.014 < stepsize/4
        if size(msd, 1) > 20
            FitRange = 10;
        else
            FitRange = round(size(msd, 1)./2);
        end
    elseif 0.014 < stepsize
        if size(msd, 1) < 100
            FitRange = round(0.16*size(msd, 1));
        else
            FitRange = round(0.07*size(msd, 1));
        end
    else
        if size(msd, 1) < 100
            FitRange = round(0.45*size(msd, 1));
        else
            FitRange = round(0.20*size(msd, 1));
        end
    end
    if FitRange >= 4
        t = (1:FitRange)*expTime;
        % toFit = log(msd(1:FitRange));
        % fitPar = fit(log(t(:)),toFit(:),'a*x+b');
        % alpha = fitPar.a;

        %% faster to solve analytically then to do linear fit
        x = log(t(:));
        y = log(msd(1:FitRange));
        xm = mean(x);
        ym = mean(y);
        alpha = sum((x - xm) .* (y - ym)) / sum((x - xm).^2);
    else
        alpha = NaN;
    end
end