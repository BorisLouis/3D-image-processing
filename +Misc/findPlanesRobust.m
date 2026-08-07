function [BW, boxes, info] = findPlanesRobust(imRaw, nPlanes, varargin)
%FINDPLANESROBUST  Locate the nPlanes rectangular planes on the camera chip.
%
%   Fallback for the case where the plane background is only marginally
%   different from the chip background (bad signal to noise) and/or where a
%   left-right illumination gradient runs over the chip.  In that situation
%   NO single global threshold can cut out the four planes.
%
%   The trick: do not rely on the plane background at all.  The emitters
%   (nanoparticles) only exist inside the planes, and a top-hat filter makes
%   them visible independently of any slowly varying background.  Their
%   column/row occupancy gives a coarse box per plane; the box edges are then
%   snapped onto the real intensity step of a flat-fielded copy of the image.
%
%   [BW, boxes] = Misc.findPlanesRobust(Raw, 4)
%
%   BW    : logical mask, true inside the nPlanes rectangles
%   boxes : nPlanes x 4  [xStart xEnd yStart yEnd]
%   info  : struct with the intermediate results (for debugging / plotting)
%
%   Name-value options (defaults tuned on 650 x 2048 chips, 4 planes):
%     'PartRadius'    10   ~ size of a diffraction limited spot [px]
%     'PartMinArea'    6   min. area of an accepted emitter [px]
%     'MinPlaneWidth'150   smallest acceptable plane width [px]
%     'MaxGapInPlane'101   largest gap in the emitter profile that is still
%                          considered "inside a plane".  MUST be smaller
%                          than the gap between two neighbouring planes.
%     'EdgeSearch'    80   +/- window used to snap an edge onto the real step
%
%   -> put this file in your +Misc package folder.

if nargin < 2 || isempty(nPlanes), nPlanes = 4; end

p = inputParser;
p.addParameter('PartRadius',    10);
p.addParameter('PartMinArea',    6);
p.addParameter('MinPlaneWidth',150);
p.addParameter('MaxGapInPlane',101);
p.addParameter('EdgeSearch',    80);
p.parse(varargin{:});
o = p.Results;

im = double(imRaw);
[nRow, nCol] = size(im);

%% 1. emitter map --------------------------------------------------------
% the top-hat removes everything that varies slower than the structuring
% element, so both the illumination gradient and the plane offset disappear
% and only the nanoparticles are left
tp   = imtophat(im, strel('disk', o.PartRadius));
medT = median(tp(:));
sigT = 1.4826 * median(abs(tp(:) - medT));      % robust sigma (MAD)
if sigT <= 0, sigT = std(tp(:)); end

P = false(nRow, nCol);
for k = [6 5 4 3]                               % relax if very few emitters
    P = bwareaopen(tp > medT + k*sigT, o.PartMinArea);
    if nnz(P) > 50*nPlanes, break; end
end

%% 2. coarse column bands from the emitter positions ---------------------
colProf  = movmean(sum(P, 1), 51);
[xs, xe] = segProfile(colProf, nPlanes, o.MaxGapInPlane, o.MinPlaneWidth);

%% 3. if a plane has too few emitters: fill it in using the geometry -----
% the planes all have the same width and a regular pitch, so a missing one
% can be reconstructed from the ones that were found
if numel(xs) < nPlanes && numel(xs) >= 2
    W = median(xe - xs);
    C = (xs + xe)/2;
    pitch = median(diff(C));
    while numel(C) < nPlanes
        if C(1) - pitch > W/2
            C = [C(1)-pitch, C];                %#ok<AGROW>
        else
            C = [C, C(end)+pitch];              %#ok<AGROW>
        end
    end
    C  = sort(C);
    xs = max(round(C - W/2), 1);
    xe = min(round(C + W/2), nCol);
end
if numel(xs) ~= nPlanes
    warning('findPlanesRobust:nPlanes', ...
        'found %d plane(s) instead of %d - check PartRadius / MaxGapInPlane.', ...
        numel(xs), nPlanes);
end

%% 4. flat field ---------------------------------------------------------
% an opening with a horizontal SE LONGER than one plane erases the planes and
% keeps only the slowly varying chip background; subtracting it removes the
% illumination gradient that makes a global threshold fail
L    = 2*round(0.75*median(xe - xs)) + 1;
L    = max(L, 301);
bgS  = imopen(im, strel('rectangle', [1 L]));
bgS  = imgaussfilt(bgS, 25);
flat = im - bgS;

%% 5. snap the column edges onto the real intensity step -----------------
rowProf  = movmean(sum(P, 2).', 51);
[ys, ye] = segProfile(rowProf, 1, o.MaxGapInPlane, o.MinPlaneWidth);
y0 = ys(1);  y1 = ye(1);

xProf = movmean(medfilt1(median(flat(y0:y1, :), 1), 25), 9);
gx    = gradient(xProf);
gThr  = 0.25 * robPrctile(abs(gx), 99.5);
for k = 1:numel(xs)
    xs(k) = snapEdge(gx, xs(k), +1, o.EdgeSearch, gThr);
    xe(k) = snapEdge(gx, xe(k), -1, o.EdgeSearch, gThr);
    if xe(k) > 2015
        xe(k) = 2015;
    end
end

%% 6. top / bottom edge: half maximum of the flat-fielded row profile ----
cols = false(1, nCol);
for k = 1:numel(xs), cols(xs(k):xe(k)) = true; end

for PlaneIdx = 1:nPlanes
    yProf = movmean(medfilt1(median(flat(:, xs(PlaneIdx):xe(PlaneIdx)), 2).', 25), 9);
    base  = max(median(yProf(1:20)), median(yProf(end-19:end)));
    ctr   = round((y0 + y1)/2);
    plate = median(yProf(max(ctr-50,1) : min(ctr+50,nRow)));
    half  = base + 0.5*(plate - base);
    yTop(PlaneIdx) = find(yProf(1:ctr) < half, 1, 'last');
    if isempty(yTop(PlaneIdx)), yTop(PlaneIdx) = 1; else, yTop(PlaneIdx) = yTop(PlaneIdx) + 1; end
    yBot(PlaneIdx) = find(yProf(ctr:end) > half, 1, 'last');
    if isempty(yBot(PlaneIdx)), yBot(PlaneIdx) = nRow; else, yBot(PlaneIdx) = ctr + yBot(PlaneIdx) - 2; end
end



%% 7. build the mask -----------------------------------------------------
BW    = false(nRow, nCol);
boxes = zeros(numel(xs), 4);
for k = 1:numel(xs)
    BW(yTop(k):yBot(k), xs(k):xe(k)) = true;
    boxes(k, :) = [xs(k) xe(k) yTop(k) yBot(k)];
end

info = struct('particleMask', P, 'flat', flat, ...
              'colProfile', colProf, 'xProfile', xProf, 'yProfile', yProf);
end
% ------------------------------------------------------------------------
function [s, e] = segProfile(prof, nWanted, maxGap, minLen)
% threshold a 1-D profile, close the holes inside a plane, keep the wide runs
s = []; e = [];
for frac = [0.12 0.08 0.05 0.20 0.30]
    b  = prof > frac*max(prof);
    b  = imclose(b, true(1, maxGap));
    b  = imopen (b, true(1, 51));
    d  = diff([false, b, false]);
    s_ = find(d ==  1);
    e_ = find(d == -1) - 1;
    keep = (e_ - s_) >= minLen;
    s_ = s_(keep);  e_ = e_(keep);
    if numel(s_) == nWanted, s = s_; e = e_; return; end
    if numel(s_) > numel(s),  s = s_; e = e_; end   % keep the best so far
end
end
% ------------------------------------------------------------------------
function idx = snapEdge(g, idx0, sgn, win, minSlope)
% move an approximate edge onto the steepest intensity step next to it
a = max(idx0 - win, 1);
b = min(idx0 + win, numel(g));
[v, k] = max(sgn * g(a:b));
if v > minSlope, idx = a + k - 1; else, idx = idx0; end
end
% ------------------------------------------------------------------------
function v = robPrctile(x, pct)
% percentile without the Statistics toolbox
x = sort(x(:));
v = x(min(max(round(pct/100*numel(x)), 1), numel(x)));
end