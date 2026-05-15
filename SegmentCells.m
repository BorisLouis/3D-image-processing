%% segment_cells.m
% Segments cells from membrane staining (CAM1 .tif image) using watershed,
% then extracts per-cell videos from CAM2 .his movie.
%
% Folder structure expected:
%   baseDir\
%     subfolderX\
%       CAM1\*.tif   (single-frame membrane staining)
%       CAM2\*.his   (multi-frame movie)
%
% Output structure:
%   outputDir\
%     subfolderX\
%       segmentation_QC.tiff       (colour overlay for visual check)
%       debug_steps.tiff           (intermediate steps, only if debugMode=true)
%       cell_01\cell_01.tiff
%       cell_02\cell_02.tiff
%       ...

clear; clc;

%% --- USER SETTINGS -----------------------------------------------------------
baseDir   = 'D:\Data Hannah\widefield_multicolor\20260318';
outputDir = 'D:\Data Hannah\widefield_multicolor\20260318_output';

% --- Segmentation parameters (tune these if results are poor) -----------------
gaussSigma    = 2;       % smoothing sigma before thresholding (px)
tophatRadius  = 40;      % ~half expected cell diameter in pixels
                         %   larger = better illumination correction for big cells
minCellArea   = 2000;    % minimum cell area (px^2) - raise to reject more small debris
maxCellArea   = 150000;  % maximum cell area (px^2) - lower to reject fused-cell blobs
dilateRadius  = 4;       % membrane gap-closing radius (px) - raise if membranes have gaps
hMinima       = 0.04;    % h-minima suppression (0-1)
                         %   raise (e.g. 0.08-0.15) if cells are over-split
borderClearPx = 5;       % ignore cells within this many px of image border
nFrames = 500;

% --- Debug mode: saves intermediate step images for the FIRST subfolder -------
debugMode     = 0;    % set false to skip after you are happy with results

%% --- DISCOVER SUBFOLDERS -----------------------------------------------------
subDirs = dir(baseDir);
subDirs = subDirs([subDirs.isdir] & ~startsWith({subDirs.name}, '.'));

if isempty(subDirs)
    error('No subfolders found in %s', baseDir);
end
fprintf('Found %d subfolders.\n', numel(subDirs));

%% --- PROCESS EACH SUBFOLDER --------------------------------------------------
for s = 1:numel(subDirs)

    subName = subDirs(s).name;
    cam1Dir = fullfile(baseDir, subName, 'CAM1');
    cam2Dir = fullfile(baseDir, subName, 'CAM2');

    fprintf('\n=== Processing: %s ===\n', subName);

    % --- Find input files -----------------------------------------------------
    cam1Files = dir(fullfile(cam1Dir, '*.tif'));
    if isempty(cam1Files)
        cam1Files = dir(fullfile(cam1Dir, '*.tiff'));
    end
    cam2Files = dir(fullfile(cam2Dir, '*.his'));

    if isempty(cam1Files)
        warning('No .tif/.tiff file in CAM1 for %s - skipping.', subName);  continue;
    end
    if isempty(cam2Files)
        warning('No .his file in CAM2 for %s - skipping.', subName);  continue;
    end

    cam1Path = fullfile(cam1Dir, cam1Files(1).name);
    cam2Path = fullfile(cam2Dir, cam2Files(1).name);

    % --- Read CAM1 & CAM2 -----------------------------------------------------
    fprintf('  Reading CAM1: %s\n', cam1Files(1).name);
    memImg = imread(cam1Path);

    fprintf('  Reading CAM2: %s\n', cam2Files(1).name);
    %v = Load.Movie.his.getFrame(cam2Path, 1:500);
    fprintf('  CAM2 frames: %d\n', nFrames);

    % --- Segment ---------------------------------------------------------------
    fprintf('  Segmenting...\n');
    doDebug = debugMode && (s == 1);   % only debug first folder
    qcDir   = fullfile(outputDir, subName);
    if ~exist(qcDir, 'dir'), mkdir(qcDir); end

    [labelMap, debugInfo] = segmentCells(memImg, gaussSigma, tophatRadius, ...
                                         minCellArea, maxCellArea, ...
                                         dilateRadius, hMinima, borderClearPx);

    % Save debug intermediates for first subfolder
    if doDebug
        saveDebugSteps(debugInfo, fullfile(qcDir, 'debug_steps.tiff'));
        fprintf('  Debug steps saved to: %s\n', qcDir);
    end

    nCells = max(labelMap(:));
    fprintf('  Detected %d cells.\n', nCells);

    if nCells == 0
        warning('No cells detected in %s - skipping video extraction.', subName);
        continue;
    end

    % QC overlay
    saveLabelQC(memImg, labelMap, fullfile(qcDir, 'segmentation_QC.tiff'));

    % --- Extract per-cell videos -----------------------------------------------
    props = regionprops(labelMap, 'BoundingBox');

    for c = 1:nCells
        cellName   = sprintf('cell_%02d', c);
        cellOutDir = fullfile(outputDir, subName, cellName);
        if ~exist(cellOutDir, 'dir'), mkdir(cellOutDir); end
        outTiff = fullfile(cellOutDir, [cellName '.tiff']);

        bb   = round(props(c).BoundingBox);
        rMin = max(1, bb(2));
        rMax = min(size(labelMap,1), bb(2)+bb(4)-1);
        cMin = max(1, bb(1));
        cMax = min(size(labelMap,2), bb(1)+bb(3)-1);

        cellMask = labelMap(rMin:rMax, cMin:cMax) == c;

        tiffObj = Tiff(outTiff, 'w');
        tagStruct.ImageLength     = rMax-rMin+1;
        tagStruct.ImageWidth      = cMax-cMin+1;
        tagStruct.BitsPerSample   = 16;
        tagStruct.SamplesPerPixel = 1;
        tagStruct.RowsPerStrip    = 16;
        tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagStruct.Photometric     = Tiff.Photometric.MinIsBlack;
        tagStruct.Compression     = Tiff.Compression.None;

        disp(append('    cell ', num2str(c), '/', num2str(nCells), ' Loading frames'))
        allFrames = Load.Movie.his.getFrame(cam2Path, 1:nFrames);
        disp(append('    cell ', num2str(c), '/', num2str(nCells), ' Cutting ROI'))
        allFrames = allFrames(rMin:rMax, cMin:cMax, :);
        allFrames(repmat(~cellMask, 1, 1, nFrames)) = 0;

        hh = waitbar(0, 'Writing frames...');
        for f = 1:nFrames
            waitbar(f/nFrames, hh, sprintf('Saving frame %d/%d', f, nFrames));
            if f == 1
                tiffObj.setTag(tagStruct);
                tiffObj.write(uint16(allFrames(:,:,f)));
            else
                tiffObj.writeDirectory();
                tiffObj.setTag(tagStruct);
                tiffObj.write(uint16(allFrames(:,:,f)));
            end
        end
        tiffObj.close();
        close(hh);
    end

    fprintf('  Saved %d cell videos to: %s\n', nCells, fullfile(outputDir, subName));
end

fprintf('\nAll done.\n');


%% =============================================================================
%  LOCAL FUNCTIONS
%% =============================================================================

function [labelMap, dbg] = segmentCells(img, gaussSigma, tophatRadius, ...
                                         minArea, maxArea, ...
                                         dilateRadius, hMin, borderPx)
    img2 = imSegmentation.adaptiveThresh(img, 4, 0.55, 10);
    img2 = bwareaopen(~img2, 150);
    se = strel('disk', 12);
    img3 = imclose(img2, se);
    bwThin  = bwskel(img3, 'MinBranchLength', 20);
    bwFinal = imdilate(bwThin, strel('disk', 2));
    interior = ~bwFinal;
    dist = bwdist(bwFinal);
    dist = -dist;
    dist = imhmin(dist, 2);
    dist(bwFinal) = Inf;
    ws = watershed(dist);
    props = regionprops(ws, 'Area');
    areas = [props.Area];
    
    % Build a size image (each pixel gets the area of its region)
    sizeImg = zeros(size(ws));
    for k = 1:numel(areas)
        sizeImg(ws == k) = areas(k);
    end
    sizeImg(or(sizeImg > 25000, sizeImg < 10000)) = 0;

    labelMap = bwlabeln(sizeImg);
    dbg = [];
end


function saveDebugSteps(dbg, outPath)
% Save each intermediate segmentation step as a page in a multi-page TIFF.
% Open with FIJI: File > Import > TIFF Virtual Stack
    fields = {'img','imgCorr','imgNoVes','memEnhanced','memMask','cellInterior','distMap'};
    titles = {'1_Normalised','2_IllumCorr','3_NoVesicles','4_MemEnhanced', ...
              '5_MemMask','6_CellInterior','7_DistMap'};

    for i = 1:numel(fields)
        data = mat2gray(double(dbg.(fields{i})));
        img8 = uint8(data * 255);
        if i == 1
            imwrite(img8, outPath, 'tiff', 'WriteMode','overwrite', ...
                    'Description',titles{i}, 'Compression','none');
        else
            imwrite(img8, outPath, 'tiff', 'WriteMode','append', ...
                    'Description',titles{i}, 'Compression','none');
        end
    end

    % Append colour label map as final page
    nCells = double(max(dbg.labelMap(:)));
    if nCells > 0
        rgb = label2rgb(dbg.labelMap, lines(nCells), 'k');
        imwrite(rgb, outPath, 'tiff', 'WriteMode','append', ...
                'Description','8_LabelMap', 'Compression','none');
    end
end


function saveLabelQC(img, labelMap, outPath)
% Save colour-outline overlay TIFF for visual quality control.
    img8   = uint8(255 * mat2gray(double(img)));
    rgb    = repmat(img8, 1, 1, 3);
    nCells = double(max(labelMap(:)));
    if nCells == 0, imwrite(rgb, outPath); return; end

    cmap = lines(nCells);
    for c = 1:nCells
        perim = bwperim(labelMap == c);
        for ch = 1:3
            layer = rgb(:,:,ch);
            layer(perim) = uint8(cmap(c,ch) * 255);
            rgb(:,:,ch)  = layer;
        end
    end
    imwrite(rgb, outPath);
end


function img = readHIS(filePath, frameIdx)
% Read one frame from a Hamamatsu .his file.
    fid = fopen(filePath, 'rb', 'l');
    if fid < 0, error('Cannot open: %s', filePath); end

    magic = fread(fid, 4, '*char')';
    if ~strcmp(strtrim(magic(1:2)), 'IM')
        fclose(fid);
        error('Not a HIS file: %s', filePath);
    end

    fseek(fid,  4, 'bof'); commentLen = fread(fid, 1, 'uint16');
    fseek(fid,  6, 'bof'); width      = fread(fid, 1, 'uint16');
    fseek(fid,  8, 'bof'); height     = fread(fid, 1, 'uint16');
    fseek(fid, 14, 'bof'); bitDepth   = fread(fid, 1, 'uint16');

    headerSize    = 64 + commentLen;
    bytesPerPixel = bitDepth / 8;
    frameSize     = width * height * bytesPerPixel;

    fseek(fid, headerSize + (frameIdx-1)*frameSize, 'bof');

    if bitDepth == 16
        raw = fread(fid, width*height, 'uint16=>uint16');
    else
        raw = fread(fid, width*height, 'uint8=>uint16');
    end
    fclose(fid);

    img = reshape(raw, [width, height])';
end


function n = getHISFrameCount(filePath)
% Return the total number of frames in a .his file.
    info = dir(filePath);
    if isempty(info), error('File not found: %s', filePath); end

    fid = fopen(filePath, 'rb', 'l');
    fseek(fid,  4, 'bof'); commentLen = fread(fid, 1, 'uint16');
    fseek(fid,  6, 'bof'); width      = fread(fid, 1, 'uint16');
    fseek(fid,  8, 'bof'); height     = fread(fid, 1, 'uint16');
    fseek(fid, 14, 'bof'); bitDepth   = fread(fid, 1, 'uint16');
    fclose(fid);

    headerSize    = 64 + commentLen;
    bytesPerPixel = bitDepth / 8;
    frameSize     = width * height * bytesPerPixel;
    n             = floor((info.bytes - headerSize) / frameSize);
end