%Split chip of each camera into an upper and lower image with each two
%images

function [chan1, chan2, idx] = splitMultiModalCamera(im)
    horizmean = mean(im,2);
    horizmeancut = horizmean(300:end-300);
    minimum = min(horizmeancut);
    idx = find(horizmean == minimum);

    chan1 = im(1:idx, :);
    chan2 = im((idx+1: end), :);
end