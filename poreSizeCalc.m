clear;
clc;
close all;

%% User Input
pxSize = 100; %in nm
Threshold = 0.7; %number between 0 and 1 (% of max intensity)
bigPores = 25; %in px ("draw a line" between small pores and big pores);
nBins = 20; %For Log histogram
% used during testing 2, normal should 244
nFrame = 244; %n of frame to analyze
%% Loading Data
% conversion from pixel to area
pxArea = pxSize*pxSize*1e-6; %in �m^2
bigPores = bigPores*pxArea;
mainFolderName = uigetdir;
assert(ischar(mainFolderName),'User canceled the selection of file, excecution aborted');
%extract the name of the current folder
idx = strfind(mainFolderName,filesep) ;
currentFolderName = mainFolderName(idx(end)+1:end) ;
%Remove dots from the name
currentFolderName = regexprep(currentFolderName,'\.','_');

%Extract the part of the folder that is a tif file
Folder_Content = dir(mainFolderName);
index2Images   = contains({Folder_Content.name},'.tif');
images2Analyze = Folder_Content(index2Images);
%% Displaying First image
stack2Load  = 1;
path2Stacks = strcat(images2Analyze(stack2Load).folder,filesep);
p2file      = strcat(path2Stacks,images2Analyze(stack2Load).name);
warning('off');
fileInfo    = loadMovie.tif.getinfo(p2file);
IM     = loadMovie.tif.getframes(p2file, 10); %Loading on of the frame
warning('on');

% smoothing the image
IM = imgaussfilt(IM,2);
%Adaptive threshold (usually better for high conc)
T = adaptthresh(IM,Threshold,'ForegroundPolarity','dark','Statistic','mean');
% binarization
BWadapt = imbinarize(IM,T);
% removing small pores
BWadapt = bwareaopen(BWadapt,4);
% smoothing
se = strel('disk',2);
BWadapt = imclose(BWadapt,se);

%No Adaptive Threshold (usually better for low conc)
BW = imbinarize(IM);
holes = ~BW;
holes = bwareaopen(holes,9);
holesAdapt = ~BWadapt;

% plot and save the figure
H0 = figure;
hold(gca,'on')
subplot(1,3,1)

imagesc(IM)
title('Raw image')
axis image

subplot(1,3,2)
imagesc(holesAdapt)
title(['Adaptive, sen: ' num2str(Threshold,2)])
axis image

subplot(1,3,3)
imagesc(holes)
title('Automated Threshold')
axis image
hold off

fileName0 = sprintf('%s%sIm0-Pores_Threshold-%d.png',mainFolderName,filesep,Threshold*100);
saveas(H0,fileName0)


answer = questdlg('Which method should be use?', ...
	'Threshold method', ...
	'Adaptive','Automated','None - cancel','Adaptive');
% Handle response
switch answer
    case 'Adaptive'
        disp([answer ' coming right up.'])
        method2Use = 'adaptThreshold';
    case 'Automated'
        disp([answer ' coming right up.'])
        method2Use = 'normThreshold';
    case 'None - cancel'
        disp([answer ' we get out now'])
        return
end

%% Looping through the Data

h = waitbar(0);
nImStacks = size(images2Analyze,1);
bins = zeros(nBins,nImStacks); %Store Bins
occurrences = bins; % Store occurences

allData = [];
for j = 1:nImStacks
    hMessage = sprintf('Loading image stack number %d/%d',j,nImStacks);
    waitbar(0,h,hMessage);
    %Data loading
    path2Stacks = strcat(images2Analyze(j).folder,filesep);
    tmpName = images2Analyze(j).name;
    p2file      = strcat(path2Stacks,tmpName);
    warning('off','all')
    fileInfo    = loadMovie.tif.getinfo(p2file);
    
    % segmented images
    segIm = false(fileInfo.Width, fileInfo.Length, nFrame);
    
    warning('on','all')
    tNframes = fileInfo.Frame_n;
    assert(tNframes>=nFrame,'you dont have the expected number of frames')
    
     
    % init data that contains all infor for a single tif file
    tifStackData = [];
    
    hMessage = sprintf('Analysis of Stack Number %d/%d',j,nImStacks);
    %loop through the frames of the current stack
    nIM = nFrame;
    for i=1:nIM
        waitbar(i/nIM,h,hMessage);
        % Loading image number i
        warning('off','all')
        IM     = loadMovie.tif.getframes(p2file, i);
        warning('on','all')
        disp(['loaded frame: ' num2str(i)])
        I  = double(IM);
        % Normalize the image
        I  = I./max(max(I));
        % Gaussian filtering (smooth image)
        I  = imgaussfilt(I,3);
        
        %%%%%%%%%%%%%%%%%%%%NOTE TO RAFA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TO DO: I tried to make it work with the leas amount of code
        %repeated for the two method, hoping we could use the same processing
        %e.g. looping through all the pore in both case. However it seems
        %quite long if we take High concentration data that usually give
        %better results with adaptive threshold and then loop through the
        %pore. (e.g. couple of minute and it was not done for a single
        %frame).
        %%%%%%%%%%%%%%%%%%%%END NOTE TO RAFA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        switch answer
            case 'Adaptive'
                disp([answer ' coming right up.'])
                % Get an adaptive threshold (~depends on local intensity levels)
                T  = adaptthresh(I,Threshold,'ForegroundPolarity','dark','Statistic','mean');
                % Binarization of image and cleaning
                BWadapt = imbinarize(I,T);
                BWadapt = bwareaopen(BWadapt,4);
                se = strel('disk',2);
                BWadapt = imclose(BWadapt,se);
                BW_pores = ~BWadapt;

            case 'Automated'
                disp([answer ' coming right up.'])
                BW = imbinarize(I);
                BW_pores = ~BW;
                BW_pores = bwareaopen(BW_pores,9); 
            otherwise
                error('Unexpected!')
        end
        
        segIm(:,:,i) = BW_pores;
        %Remove pores that are in touch with the edges(uncertainty about
        %size)
        %Need to change this so it does not always clear pores
%         BW_pores = imclearborder(BW_pores);
        [L,n] = bwlabel(BW_pores);
              
        % regions data is: width, size, solidity, isBorder
        regData = zeros(n,4);
%         hi = waitbar(0,'iterating over regions');
        
        % the par for is not super clean at the moment as it will also init
        % the workes, I can fix that later on
        
        % Iterating over each region
        parfor k = 1:n %parfor can be placed here
            
            tmpBW = L==k;
            % this one is the time consuming step -  are there any faster
            % options?
            % calculating width
            [ fWidth ] = SDcalc.fastWidthBW( tmpBW );
            
            % calculating area
            tmpSize = sum(sum(tmpBW));
            
            % getting boundary
            [B] = bwboundaries(tmpBW,'noholes');
            % make sure we get the largest shape in case there are children
            Blength = cellfun(@length,B);
            [~, idx] = maxk(Blength,2);
            B = B(idx);
            Blength = Blength(idx);
            assert(length(B)==1,'problems');
            boundary = B{1};
            
            % calculating solidity
            % test for colinearity
            tmp = sum(abs(diff(boundary)));
            isColi = any(tmp==0);
            if and(Blength >3,~isColi)
                [ vals, names ] = SDcalc.solidity( boundary' );
                sol = vals{1};
            else
                % if it is small or colinear then sol is 1
                sol = 1;
            end
            
            % checking if contour is at the border
            test = imclearborder(tmpBW);
            isBorder = sum(test(:)) == 0;
            
            % clean way of indexing so parfor works
            tmpOut = [fWidth, tmpSize, sol, isBorder];
            regData(k,:) = tmpOut;
            
            disp(['Done for region ' num2str(k) '/' num2str(n)...
                   '; with b: ' num2str(length(boundary))])
%             waitbar(k / n)
        end
        disp(['Done for Image ' num2str(i) '/' num2str(nIM)])
        disp('---------------------NEXT FRAME ----------')
%         close(hi)
%         %Get properties of the pores on the image.
%         stats = regionprops(BW_pores,'Area','MajorAxisLength',...
%         'MinorAxisLength');
        regData(:,1:2) = regData(:,1:2).*pxArea;
        % add info about image index
        regData = cat(2,ones(n,1).*i,regData);
        tifStackData = cat(1,tifStackData, regData);
    end
    
    % add info about tif index
    tifStackData = cat(2,ones(size(tifStackData,1),1).*j,tifStackData);
    
    disp(['Done for Tif file ' num2str(j) '/' num2str(nImStacks) ', now saving'])
    
    % now we save segmented images
    tifName = sprintf('%s%sSegmented_%s',mainFolderName,filesep,tmpName);
    t = Tiff(tifName, 'w');
    setTag(t,'ImageLength',size(segIm,1))
    setTag(t,'ImageWidth',size(segIm,2))
    setTag(t,'Photometric',Tiff.Photometric.MinIsBlack)
    setTag(t,'BitsPerSample',1)
    setTag(t,'SamplesPerPixel',1)
    setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
    t.write(segIm(:,:,1))

    for i = 2:size(segIm,3)
%         disp(i)
        t.writeDirectory
        setTag(t,'ImageLength',size(segIm,1))
        setTag(t,'ImageWidth',size(segIm,2))
        setTag(t,'Photometric',Tiff.Photometric.MinIsBlack)
        setTag(t,'BitsPerSample',1)
        setTag(t,'SamplesPerPixel',1)
        setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
        t.write(segIm(:,:,i))
    end
    t.close 

    disp('---------------------NEXT TIF ----------')
        
    allData = cat(1,allData, tifStackData);

end
close(h);

%%
T = array2table(allData,...
    'VariableNames',{'TifIDX','ImageIDX','Width','Area', 'Solidity','IsAtBorder'});

save([ path2Stacks 'poreProps.mat'],'T')

h = msgbox('The Data were succesfully saved !', 'Success');
return
%% Plotting
totalV = T.Area;
totalV = totalV(:);
[CDF,CCDF] = Misc.getCDF(totalV);

figure()
plot(CCDF.x,CCDF.y)
a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF for AREA',' for all tif files in folder'})
xlabel('Pore size')
ylabel('CCDF - prob [0 1]')
a.FontSize = 14;
 
 
%%
% H1 = figure(10);
% hold (gca,'on')
% title(currentFolderName)
% set(gca,'YScale','log');
% set(gca,'XScale','log');
% ylabel('Number of Pores');
% xlabel('Pore area');
% leg = cell(1,nImStacks);
% for i = 1:nImStacks
%     plot(bins(:,i),occurrences(:,i));
%     leg{i} = sprintf('Stack %d',i);
% 
% end
% 
% legend(leg)
% hold (gca,'off')
% 
% 
% H2 = figure(11);
% hold(gca,'on')
% title(sprintf('%s - Overview',currentFolderName))
% scatter(median(bins,2),median(occurrences,2))
% errorbar(median(bins,2),median(occurrences,2),std(occurrences,1,2),'LineStyle','none')
% set(gca,'YScale','log');
% set(gca,'XScale','log');
% ylabel('Number of Pores');
% xlabel('Pore area');
% legend('median','Standard deviation')
% 
% hold(gca,'off')
% 
% histData(1).medBins = median(bins,2);
% histData(1).medOcc  = median(occurrences,2);
% histData(1).STD     = std(occurrences,1,2);
% %% Saving figures & Data
% 
% 
% fileName1 = sprintf('%s%s%s-AllCurves',mainFolderName,'\',currentFolderName);
% savefig(H1,fileName1)
% 
% fileName2 = sprintf('%s%s%s-AverageCurve',mainFolderName,'\',currentFolderName);
% savefig(H2,fileName2)
% 
% fileNameMat = sprintf('%s%s%s-histData',mainFolderName,'\',currentFolderName);
% save(fileNameMat,'histData');
% 

%%

