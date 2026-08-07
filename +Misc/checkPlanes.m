function [ok, BW, boxes] = checkPlanes(BW, nPlanes, minExtent, maxAreaRatio)
%CHECKPLANES  Sanity check on a binarised chip image.
%   Returns ok = true only when BW contains exactly nPlanes blobs that look
%   like the imaging planes (large, similar in size, nearly rectangular).
%   BW is returned with all the small junk removed, so only the plane blobs
%   are kept.
%
%   [ok, BW, boxes] = Misc.checkPlanes(BW, 4)
%
%   boxes : nPlanes x 4  [xStart xEnd yStart yEnd]   (empty when ok == false)
%
%   -> put this file in your +Misc package folder.

if nargin < 2 || isempty(nPlanes),      nPlanes      = 4;    end
if nargin < 3 || isempty(minExtent),    minExtent    = 0.75; end   % area / boundingbox area
if nargin < 4 || isempty(maxAreaRatio), maxAreaRatio = 1.6;  end   % largest / smallest plane

ok    = false;
boxes = [];

CC = bwconncomp(BW);
if CC.NumObjects == 0, BW = false(size(BW)); return; end

S = regionprops(CC, 'Area', 'BoundingBox', 'Extent');
A = [S.Area];

% a plane covers at least ~1% of the chip and is of the same order of
% magnitude as the biggest object that was found
isBig = A > 0.25*max(A) & A > 0.01*numel(BW);

BW = ismember(labelmatrix(CC), find(isBig));   % drop the small junk
if ~any(isBig), return; end

E  = [S(isBig).Extent];
Ab = A(isBig);
ok = (sum(isBig) == nPlanes) && all(E > minExtent) && (max(Ab)/min(Ab) < maxAreaRatio);

if ok
    bb    = reshape([S(isBig).BoundingBox], 4, []).';       % [x y w h]
    boxes = [ceil(bb(:,1)), ceil(bb(:,1))+bb(:,3)-1, ...
             ceil(bb(:,2)), ceil(bb(:,2))+bb(:,4)-1];
    boxes = sortrows(boxes, 1);
end
end