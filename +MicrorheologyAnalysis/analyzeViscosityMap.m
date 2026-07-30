function S = analyzeViscosityMap(ViscMap, varargin)
%ANALYZEVISCOSITYMAP  Quantify the spatial structure of a masked viscosity map.
%
%   S = ANALYZEVISCOSITYMAP(ViscMap) computes a set of scalar and vector
%   descriptors of the spatial heterogeneity of a 2-D viscosity map, where
%   pixels outside the hand-drawn ROI (and any excluded regions such as the
%   nucleus) are NaN. Only non-NaN ("in-ROI") pixels are analysed.
%
%   S = ANALYZEVISCOSITYMAP(ViscMap, 'Name', Value, ...) allows options:
%
%     'PixelSize'         Physical size of one pixel (e.g. um/px). Default 1
%                         (results reported in pixel units).
%     'ClipNegative'      Floor small negative noise values (common near the
%                         detection limit) to 0 before analysis. Default true.
%     'PatchPercentile'   Percentile of in-ROI values used to define the
%                         headline "high-viscosity patch" mask. Default 75.
%     'ThresholdSweep'    Vector of percentiles used to build the patch-size
%                         spectrum (num. patches / area vs. threshold).
%                         Default 50:5:95.
%     'NumGrayLevels'     Number of quantization levels used for the GLCM
%                         texture features. Default 16.
%     'NVariogramPairs'   Number of random pixel pairs used to estimate the
%                         empirical semivariogram (spatial correlation
%                         length). Default 200000.
%     'VariogramNBins'    Number of distance bins for the semivariogram.
%                         Default 40.
%     'VariogramMaxLagFrac'  Fraction of the max observed pixel-pair distance
%                         used to set the bin range (needs to be large enough
%                         to see any "hole effect" turnaround). Default 0.9.
%     'VariogramSmoothWindow'  Moving-average window (in bins) used to
%                         smooth the semivariogram before detecting the
%                         first peak/turning point. Default 3.
%     'VariogramPeakLookahead'  Number of subsequent bins that must show a
%                         sustained decline to confirm a genuine peak (vs.
%                         sampling noise). Default 3.
%     'PlotResults'       If true, produce a diagnostic figure. Default false.
%     'Title'             Title used on the diagnostic figure / stored in S.Meta.
%
%   OUTPUT  S is a struct with the following groups of fields:
%
%   S.Meta          Options used, ROI pixel count, image size, pixel size.
%
%   S.Basic         Distributional statistics of the in-ROI values:
%                   N, Mean, Median, SD, CV (SD/Mean), MAD, IQR, P5, P95,
%                   RobustRange (P95-P5), Skewness, Kurtosis.
%                   -> "How homogeneous is the map, overall?"
%
%   S.Heterogeneity Gini            Gini coefficient of the in-ROI values
%                                   (0 = perfectly even, ->1 = concentrated
%                                   in a few pixels).
%                   Entropy         Shannon entropy of the value histogram,
%                                   normalised to [0,1] by log2(nbins)
%                                   (1 = maximally spread distribution).
%                   -> "How unevenly is viscosity distributed, ignoring
%                      spatial position?"
%
%   S.Spatial       MoranI          Global Moran's I spatial autocorrelation
%                                   (4-connected neighbours). ~1 = smoothly
%                                   varying/clustered, ~0 = spatially random
%                                   ("salt and pepper"), <0 = checkerboard.
%                   LocalVarMean    Mean local variance in a 5x5 window
%                                   (fine-scale roughness/graininess).
%                   CorrLength_px   Spatial correlation length (range of the
%                                   empirical semivariogram), in pixels and
%                                   in physical units if PixelSize is set.
%                                   For ROIs with a cut-out nucleus (or any
%                                   hole), the semivariogram/correlation curve
%                                   is automatically truncated at its first
%                                   peak/turning point before the range is
%                                   estimated, so the "hole effect" bump
%                                   caused by the missing region does not
%                                   bias the result (see Variogram below).
%                   Variogram       Struct with the (truncated) curve used
%                                   for CorrLength_px: .dist, .gamma,
%                                   .gammaSmooth, .corr (= 1-gamma/sill,
%                                   directly plottable as "correlation vs.
%                                   distance"), .sill, .holeEffectDetected,
%                                   and .full (untruncated curve, for
%                                   reference/debugging).
%                   GLCM            Haralick texture features (Contrast,
%                                   Correlation, Energy, Homogeneity, Entropy)
%                                   from a mask-aware gray-level co-occurrence
%                                   matrix.
%                   -> "Are similar values clustered together, and at what
%                      length scale?"
%
%   S.Patches       At the headline threshold (PatchPercentile):
%                   Threshold, NumPatches, PatchAreas_px, MeanArea_px,
%                   MedianArea_px, CVArea, MaxArea_px, AreaFraction,
%                   Prominence (per patch, peak-value minus local
%                   surrounding ring value), MeanProminence, MaxProminence,
%                   PatchDensity_per1000px, ClarkEvansR (nearest-neighbour
%                   index of patch centroids: <1 clustered, ~1 random,
%                   >1 regularly spread).
%                   Sweep           Same summary computed across
%                                   ThresholdSweep percentiles, for
%                                   plotting a "patch-size spectrum".
%                   -> "How many patches, how big, how prominent, how are
%                      they arranged?"
%
%   S.Spread        RadialTrend     Correlation between distance-from-ROI-
%                                   centroid and value (>0 = higher viscosity
%                                   toward the periphery, <0 = toward the
%                                   centre).
%                   CentroidOffsetNorm  Distance between the value-weighted
%                                   centroid and the geometric ROI centroid,
%                                   normalised by the ROI's equivalent radius
%                                   (0 = hot spots centred on the ROI, larger
%                                   = hot spots skewed to one side).
%                   -> "Where in the cell is the signal concentrated?"
%
%   Example:
%       load ViscosityMap.mat            % contains ViscMap
%       S = analyzeViscosityMap(ViscMap, 'PlotResults', true, 'Title', 'Exp1');
%
%   Requires Image Processing Toolbox (bwconncomp, regionprops, imdilate).
%
%   See also COMPAREVISCOSITYMAPS.

% ---- parse inputs -------------------------------------------------------
p = inputParser;
addParameter(p, 'PixelSize', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'ClipNegative', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PatchPercentile', 75, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'ThresholdSweep', 50:5:95, @isnumeric);
addParameter(p, 'NumGrayLevels', 16, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'NVariogramPairs', 200000, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'VariogramNBins', 40, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'VariogramMaxLagFrac', 0.9, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'VariogramSmoothWindow', 3, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'VariogramPeakLookahead', 3, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'PlotResults', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Title', '', @ischar);
parse(p, varargin{:});
opt = p.Results;

% ---- prepare image & mask -----------------------------------------------
ViscMap = double(ViscMap);
mask = ~isnan(ViscMap);
if ~any(mask(:))
    error('analyzeViscosityMap:noROI', 'ViscMap contains no non-NaN (in-ROI) pixels.');
end

img = ViscMap;
if opt.ClipNegative
    img(mask & img < 0) = 0;
end

vals = img(mask);
N = numel(vals);

% ==========================================================================
% S.Basic
% ==========================================================================
Basic = struct();
Basic.N        = N;
Basic.Mean     = mean(vals);
Basic.Median   = median(vals);
Basic.SD       = std(vals);
Basic.CV       = Basic.SD / Basic.Mean;
Basic.MAD      = median(abs(vals - Basic.Median));
Basic.P5       = localPercentile(vals, 5);
Basic.P95      = localPercentile(vals, 95);
Basic.IQR      = localPercentile(vals, 75) - localPercentile(vals, 25);
Basic.RobustRange = Basic.P95 - Basic.P5;
[Basic.Skewness, Basic.Kurtosis] = localSkewKurt(vals);

% ==========================================================================
% S.Heterogeneity  (value distribution only, no spatial information)
% ==========================================================================
Heterogeneity = struct();
Heterogeneity.Gini    = localGini(vals);
Heterogeneity.Entropy = localShannonEntropy(vals, 32);

% ==========================================================================
% S.Spatial  (texture / spatial autocorrelation)
% ==========================================================================
Spatial = struct();
Spatial.MoranI       = localMoranI(img, mask);
Spatial.LocalVarMean = localLocalVariance(img, mask, 5);
[Spatial.CorrLength_px, vgram] = localVariogramRange(img, mask, opt);
Spatial.CorrLength_phys = Spatial.CorrLength_px * opt.PixelSize;
Spatial.Variogram    = vgram; % .dist, .gamma, .sill  (for plotting)
Spatial.GLCM         = localGLCM(img, mask, opt.NumGrayLevels);

% ==========================================================================
% S.Patches  (connected-component / blob analysis)
% ==========================================================================
Patches = localPatchAnalysis(img, mask, opt.PatchPercentile);
Patches.Sweep = localPatchSweep(img, mask, opt.ThresholdSweep);

% ==========================================================================
% S.Spread  (where in the ROI the signal sits)
% ==========================================================================
Spread = localRadialProfile(img, mask);

% ==========================================================================
% S.Meta
% ==========================================================================
Meta = struct();
Meta.ImageSize   = size(ViscMap);
Meta.ROIPixels   = N;
Meta.PixelSize   = opt.PixelSize;
Meta.Options     = opt;
Meta.Title       = opt.Title;
Meta.AnalysisDate = datestr(now);

S = struct('Meta', Meta, 'Basic', Basic, 'Heterogeneity', Heterogeneity, ...
           'Spatial', Spatial, 'Patches', Patches, 'Spread', Spread);

% ---- optional diagnostic plot -------------------------------------------
if opt.PlotResults
    localPlotDiagnostics(img, mask, S, opt);
end

end % analyzeViscosityMap


% ==========================================================================
% LOCAL FUNCTIONS
% ==========================================================================

function q = localPercentile(x, p)
% Linear-interpolation percentile, no Statistics Toolbox required.
% p may be a scalar or vector of percentiles (0-100).
x = sort(x(:));
n = numel(x);
if n == 1, q = repmat(x(1), size(p)); return; end
r = (p/100) * (n - 1) + 1;
lo = floor(r); hi = ceil(r);
lo(lo < 1) = 1;
hi(hi > n) = n;
w = r - lo;
q = (1 - w) .* x(lo) + w .* x(hi);
end


function [sk, ku] = localSkewKurt(x)
% Population skewness / (excess) kurtosis, no Statistics Toolbox required.
mu = mean(x);
sd = std(x, 1); % population SD (N divisor)
sk = mean((x - mu).^3) / sd^3;
ku = mean((x - mu).^4) / sd^4 - 3; % excess kurtosis
end


function g = localGini(x)
% Gini coefficient of a non-negative vector.
x = sort(x(:));
x(x < 0) = 0; % Gini undefined for negative values; should be pre-clipped
n = numel(x);
cumx = cumsum(x);
if cumx(end) == 0
    g = 0;
    return;
end
g = (n + 1 - 2 * sum(cumx) / cumx(end)) / n;
end


function H = localShannonEntropy(x, nbins)
% Shannon entropy of the value histogram, normalised to [0,1].
counts = histcounts(x, nbins);
p = counts / sum(counts);
p = p(p > 0);
H = -sum(p .* log2(p)) / log2(nbins);
end


function I = localMoranI(img, mask)
% Global Moran's I using 4-connected (rook) neighbours, ROI-aware.
xbar = mean(img(mask));
dev = img - xbar;

num = 0; Wsum = 0;

% Right neighbour
a = mask(:, 1:end-1) & mask(:, 2:end);
da = dev(:, 1:end-1); db = dev(:, 2:end);
num  = num  + 2 * sum(da(a) .* db(a));
Wsum = Wsum + 2 * nnz(a);

% Down neighbour
a2 = mask(1:end-1, :) & mask(2:end, :);
da2 = dev(1:end-1, :); db2 = dev(2:end, :);
num  = num  + 2 * sum(da2(a2) .* db2(a2));
Wsum = Wsum + 2 * nnz(a2);

denom = sum(dev(mask).^2);
n = nnz(mask);
I = (n / Wsum) * (num / denom);
end


function v = localLocalVariance(img, mask, win)
% Mean local variance in a win-by-win window, ignoring out-of-ROI pixels.
x = img; x(~mask) = 0;
m = double(mask);
k = ones(win);

sumX  = conv2(x, k, 'same');
sumX2 = conv2(x.^2, k, 'same');
cnt   = conv2(m, k, 'same');

cnt(cnt == 0) = NaN;
meanLocal   = sumX ./ cnt;
meanSqLocal = sumX2 ./ cnt;
varLocal = meanSqLocal - meanLocal.^2;
varLocal(varLocal < 0) = 0; % guard against floating point noise

valid = mask & (cnt >= 0.9 * win * win); % require mostly-full windows
v = mean(varLocal(valid), 'omitnan');
end


function [rangeEst, vgram] = localVariogramRange(img, mask, opt)
% Empirical semivariogram via random pixel-pair sampling.
%
% For a simply-connected ROI, semivariance gamma(d) rises monotonically
% with distance d and plateaus at the sill (~sample variance); "range" is
% the distance at which it first reaches ~95% of the sill.
%
% For an ROI with a hole in it (e.g. a cut-out nucleus), gamma(d) typically
% rises, reaches a first peak, DROPS, then rises again ("hole effect")
% because pixel pairs at those larger distances are geometrically forced
% to sit on opposite sides of the hole rather than being freely sampled
% from the whole domain. That second hump is a shape artefact of the ROI
% geometry, not additional spatial structure, so both the plotted curve
% and the range estimate below are truncated at the first genuine peak.

npairs = opt.NVariogramPairs;
[ys, xs] = find(mask);
vals = img(mask);
n = numel(vals);

npairs = min(npairs, round(n * (n - 1) / 2)); % can't exceed unique pairs
i1 = randi(n, npairs, 1);
i2 = randi(n, npairs, 1);
same = (i1 == i2);
i2(same) = mod(i2(same), n) + 1;

dist = sqrt((ys(i1) - ys(i2)).^2 + (xs(i1) - xs(i2)).^2);
sqd  = (vals(i1) - vals(i2)).^2;

nbins = opt.VariogramNBins;
maxLag = opt.VariogramMaxLagFrac * max(dist);
edges = linspace(0, maxLag, nbins + 1);
binIdx = discretize(dist, edges);

gammaRaw = nan(nbins, 1);
for k = 1:nbins
    inBin = (binIdx == k);
    if any(inBin)
        gammaRaw(k) = 0.5 * mean(sqd(inBin));
    end
end
centers = (edges(1:end-1) + edges(2:end)) / 2;

% Fill sparsely-populated bins (few pixel pairs at large lags) and smooth
% before looking for a turning point, so sampling noise isn't mistaken for
% the hole-effect peak.
gammaFilled = fillmissing(gammaRaw, 'linear', 'EndValues', 'nearest');
w = max(1, opt.VariogramSmoothWindow);
gammaSmooth = movmean(gammaFilled, w);

% First genuine peak: a bin higher than the previous one AND higher than
% the next `lookahead` bins (i.e. followed by a sustained decline). If no
% such peak exists, the curve is a normal monotonic-to-plateau shape and
% nothing is truncated.
lookahead = opt.VariogramPeakLookahead;
kPeak = numel(gammaSmooth);
holeDetected = false;
for k = 2:(numel(gammaSmooth) - lookahead)
    following = gammaSmooth(k+1 : k+lookahead);
    if gammaSmooth(k) > gammaSmooth(k-1) && all(gammaSmooth(k) > following)
        kPeak = k;
        holeDetected = true;
        break;
    end
end

% Truncate the curve (for plotting) and the range/sill estimate (for the
% headline CorrLength number) to only the short-range, first-peak portion.
distT        = centers(1:kPeak);
gammaT       = gammaFilled(1:kPeak);
gammaSmoothT = gammaSmooth(1:kPeak);

sill = gammaSmoothT(end); % plateau/peak value reached within the trusted range
target = 0.95 * sill;
reached = find(gammaT >= target, 1, 'first');
if isempty(reached)
    rangeEst = distT(end); % doesn't clearly plateau within the trusted range
else
    rangeEst = distT(reached);
end

vgram = struct();
vgram.dist = distT;
vgram.gamma = gammaT;
vgram.gammaSmooth = gammaSmoothT;
vgram.corr = 1 - gammaT / sill; % correlation-style curve: ~1 at d=0, decays
vgram.sill = sill;
vgram.holeEffectDetected = holeDetected;
vgram.full = struct('dist', centers, 'gamma', gammaFilled); % untruncated, for reference
end


function G = localGLCM(img, mask, Ng)
% Mask-aware gray-level co-occurrence matrix (distance 1, averaged over
% 0/45/90/135 degrees) and standard Haralick features.
vals = img(mask);
edges = localPercentile(vals, linspace(0, 100, Ng + 1));
edges = unique(edges);
if numel(edges) < 3
    edges = linspace(min(vals), max(vals) + eps, Ng + 1);
end
Ng = numel(edges) - 1;

Q = nan(size(img));
Q(mask) = discretize(vals, edges);
Q(isnan(Q)) = -1; % sentinel for "outside ROI" (discretize/NaN both land here)

glcm = zeros(Ng, Ng);
offsets = [0 1; 1 1; 1 0; 1 -1];
[H, W] = size(img);
for k = 1:size(offsets, 1)
    dy = offsets(k, 1); dx = offsets(k, 2);
    y0a = max(1, 1 - dy); y1a = min(H, H - dy);
    x0a = max(1, 1 - dx); x1a = min(W, W - dx);
    a = Q(y0a:y1a, x0a:x1a);
    b = Q(y0a + dy:y1a + dy, x0a + dx:x1a + dx);
    valid = (a > 0) & (b > 0);
    ai = a(valid); bi = b(valid);
    idx = (ai - 1) * Ng + bi;
    counts = accumarray(idx(:), 1, [Ng * Ng, 1]);
    glcm = glcm + reshape(counts, Ng, Ng)';
end
glcm = glcm + glcm'; % symmetric
P = glcm / sum(glcm(:));

[ii, jj] = ndgrid(1:Ng, 1:Ng);
muI = sum(ii(:) .* P(:));
muJ = sum(jj(:) .* P(:));
sigI = sqrt(sum(((ii(:) - muI).^2) .* P(:)));
sigJ = sqrt(sum(((jj(:) - muJ).^2) .* P(:)));

G = struct();
G.Contrast     = sum(P(:) .* (ii(:) - jj(:)).^2);
G.Energy       = sum(P(:).^2);
G.Homogeneity  = sum(P(:) ./ (1 + (ii(:) - jj(:)).^2));
G.Correlation  = sum(P(:) .* (ii(:) - muI) .* (jj(:) - muJ)) / (sigI * sigJ);
nz = P(:) > 0;
G.Entropy      = -sum(P(nz) .* log2(P(nz)));
G.NumGrayLevels = Ng;
end


function Patches = localPatchAnalysis(img, mask, pctile)
% Connected-component analysis of the "high viscosity" mask defined at the
% given percentile of in-ROI values.
vals = img(mask);
thr = localPercentile(vals, pctile);
bw = (img > thr) & mask;

CC = bwconncomp(bw, 8);
n = CC.NumObjects;
roiArea = nnz(mask);

Patches = struct();
Patches.Threshold = thr;
Patches.NumPatches = n;

if n == 0
    Patches.PatchAreas_px = [];
    Patches.MeanArea_px = NaN; Patches.MedianArea_px = NaN;
    Patches.CVArea = NaN; Patches.MaxArea_px = NaN;
    Patches.AreaFraction = 0;
    Patches.Prominence = []; Patches.MeanProminence = NaN; Patches.MaxProminence = NaN;
    Patches.PatchDensity_per1000px = 0;
    Patches.ClarkEvansR = NaN;
    return;
end

rp = regionprops(CC, 'Area', 'Centroid');
areas = [rp.Area]';
Patches.PatchAreas_px = areas;
Patches.MeanArea_px = mean(areas);
Patches.MedianArea_px = median(areas);
Patches.CVArea = std(areas) / mean(areas);
Patches.MaxArea_px = max(areas);
Patches.AreaFraction = nnz(bw) / roiArea;
Patches.PatchDensity_per1000px = 1000 * n / roiArea;

% Prominence: peak value in patch minus median value of a ring of pixels
% immediately surrounding the patch (still inside the ROI). This measures
% how much a patch stands out from its own local surroundings, independent
% of the global threshold.
se = strel('disk', 3);
prom = nan(n, 1);
for k = 1:n
    patchMask = false(size(img));
    patchMask(CC.PixelIdxList{k}) = true;
    peakVal = max(img(patchMask));
    ring = imdilate(patchMask, se) & mask & ~patchMask;
    if any(ring(:))
        prom(k) = peakVal - median(img(ring));
    end
end
Patches.Prominence = prom;
Patches.MeanProminence = mean(prom, 'omitnan');
Patches.MaxProminence = max(prom);

% Clark-Evans nearest-neighbour index of patch centroids: compares the
% observed mean nearest-neighbour distance to that expected under complete
% spatial randomness over the ROI area. <1 clustered, ~1 random, >1 regular.
if n >= 2
    C = cat(1, rp.Centroid); % [x y]
    Patches.ClarkEvansR = localClarkEvans(C, roiArea);
else
    Patches.ClarkEvansR = NaN;
end

end


function R = localClarkEvans(centroids, roiArea)
n = size(centroids, 1);
Dx = centroids(:, 1) - centroids(:, 1)';
Dy = centroids(:, 2) - centroids(:, 2)';
D = sqrt(Dx.^2 + Dy.^2);
D(1:n+1:end) = Inf;
nnDist = min(D, [], 2);
obsMean = mean(nnDist);
density = n / roiArea;
expMean = 0.5 / sqrt(density);
R = obsMean / expMean;
end


function sweep = localPatchSweep(img, mask, pctiles)
% Patch count / area statistics across a range of thresholds, for
% comparing the overall "patch-size spectrum" between maps.
vals = img(mask);
roiArea = nnz(mask);
np = numel(pctiles);

numPatches   = nan(np, 1);
meanArea     = nan(np, 1);
areaFraction = nan(np, 1);

for k = 1:np
    thr = localPercentile(vals, pctiles(k));
    bw = (img > thr) & mask;
    CC = bwconncomp(bw, 8);
    numPatches(k) = CC.NumObjects;
    if CC.NumObjects > 0
        rp = regionprops(CC, 'Area');
        meanArea(k) = mean([rp.Area]);
    else
        meanArea(k) = NaN;
    end
    areaFraction(k) = nnz(bw) / roiArea;
end

sweep = struct('Percentile', pctiles(:), 'NumPatches', numPatches, ...
               'MeanArea_px', meanArea, 'AreaFraction', areaFraction);
end


function Spread = localRadialProfile(img, mask)
% Radial trend and centroid offset describing where in the ROI signal
% is concentrated.
[ys, xs] = find(mask);
vals = img(mask);

cy = mean(ys); cx = mean(xs); % geometric ROI centroid
d = sqrt((ys - cy).^2 + (xs - cx).^2);

cc = corrcoef(d, vals);
radialTrend = cc(1, 2);

% value-weighted centroid ("centre of mass" of viscosity signal)
wy = sum(ys .* vals) / sum(vals);
wx = sum(xs .* vals) / sum(vals);
offset = sqrt((wy - cy)^2 + (wx - cx)^2);
equivRadius = sqrt(nnz(mask) / pi);

Spread = struct();
Spread.RadialTrend = radialTrend;
Spread.CentroidOffsetNorm = offset / equivRadius;
Spread.ROICentroid = [cx, cy];
Spread.WeightedCentroid = [wx, wy];
Spread.EquivRadius_px = equivRadius;

% Binned radial profile (10 bins) for plotting
nbins = 10;
edges = linspace(0, max(d), nbins + 1);
binIdx = discretize(d, edges);
profile = nan(nbins, 1);
for k = 1:nbins
    inBin = (binIdx == k);
    if any(inBin), profile(k) = mean(vals(inBin)); end
end
Spread.RadialProfile = profile;
Spread.RadialProfileCenters = (edges(1:end-1) + edges(2:end)) / 2;
end


function localPlotDiagnostics(img, mask, S, opt)
figure('Name', ['Viscosity map diagnostics: ' opt.Title], 'Position', [100 100 1200 800]);

imgShow = img; imgShow(~mask) = NaN;

subplot(2,3,1);
imagesc(imgShow, 'AlphaData', mask); axis image; colorbar;
title(sprintf('%s\nmean %.3g \\pm %.3g', opt.Title, S.Basic.Mean, S.Basic.SD));

subplot(2,3,2);
histogram(img(mask), 40);
xlabel('Viscosity'); ylabel('Count');
title(sprintf('Gini=%.2f, Entropy=%.2f, CV=%.2f', S.Heterogeneity.Gini, S.Heterogeneity.Entropy, S.Basic.CV));

subplot(2,3,3);
thr = S.Patches.Threshold;
bw = (img > thr) & mask;
imagesc(bw); axis image; colormap(gca, gray);
title(sprintf('Patches @ P%g: n=%d, meanArea=%.1fpx', ...
    opt.PatchPercentile, S.Patches.NumPatches, S.Patches.MeanArea_px));

subplot(2,3,4);
plot(S.Spatial.Variogram.dist, S.Spatial.Variogram.corr, 'o-');
yline(0, '--');
xlabel('Distance (px)'); ylabel('Spatial correlation');
if S.Spatial.Variogram.holeEffectDetected
    holeStr = ' (cut before hole-effect rebound)';
else
    holeStr = '';
end
title(sprintf('Corr. length \\approx %.1f px%s', S.Spatial.CorrLength_px, holeStr));

subplot(2,3,5);
plot(S.Spread.RadialProfileCenters, S.Spread.RadialProfile, 'o-');
xlabel('Distance from ROI centroid (px)'); ylabel('Mean viscosity');
title(sprintf('Radial trend r=%.2f, centroid offset=%.2f', ...
    S.Spread.RadialTrend, S.Spread.CentroidOffsetNorm));

subplot(2,3,6);
plot(S.Patches.Sweep.Percentile, S.Patches.Sweep.NumPatches, 'o-');
xlabel('Threshold percentile'); ylabel('Number of patches');
title('Patch-size spectrum');

end
