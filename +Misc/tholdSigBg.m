function [tHold] = tholdSigBg(bg,sig)
%THOLDSIGBG finds optimal treshhold that separates two distributions one
%comming from background and another from signal
%   Detailed explanation goes here

sigV = sig(:);
bgV  = bg(:);
assert(mean(bgV)<mean(sigV),'strange your background looks larger than your signal')

% find CDF of signal and CCDF of background
[sigCDF] = Plotting.getCDF(sigV);
[~, bgCCDF]  = Plotting.getCDF(bgV);

% find optimal tHold that separates two distributions
allX = [sigCDF.x; bgCCDF.x];
minX = min(allX);
maxX = max(allX);
globXq = linspace(double(minX),double(maxX),length(allX)*3);
sigVq = interp1(double(sigCDF.x),double(sigCDF.y),globXq);
bgVq  = interp1(double(bgCCDF.x),double(bgCCDF.y),globXq);
diffVq = abs(sigVq-bgVq);
%not sure what was this for:
if mean(sigCDF.x(sigCDF.y>0.95))>2*max(bg)
    tHold = max(bg);
else
    [~,idx] = nanmin(diffVq);
    val = max(bgV);
    tHold = val;
end
%try to take the threshold that remove 99.9% of background

Minprominence = 5*10^(-7);
[data, p] = islocalmin(smoothdata(diff(sigCDF.y(find(sigCDF.x > 120, 1, "first"):end))), 'MinProminence', Minprominence);
while all(data == 0)
    Minprominence = Minprominence ./10;
    [data, p] = islocalmin(smoothdata(diff(sigCDF.y(find(sigCDF.x > 120, 1, "first"):end))), 'MinProminence', Minprominence);
end

Idx = find(data == 1, 1, "first");
tHold = sigCDF.x(Idx);
%previous
%tHold = mean(bg)+3*std(bg);
if tHold < 110
    test = abs(diff(medfilt1(sig, 50)));
    test(test < 10) = [];
    tHold = 100 + mean(test);
end
end

