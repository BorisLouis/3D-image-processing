FolderName = 'S:\Indra\20250710\Cam1\All data\mSi';
Folder = dir(FolderName);


for i = 3:size(Folder, 1)
    File = append(FolderName, filesep, Folder(i).name);
    load(append(File, '\SegmentMovie\SegmentMovie5'));
    Frame = Load.Movie.tif.getframes(append(File, '\calibrated1\calibratedPlane1.tif'), 500);
    Mask = Mask{500,1};
    Fig = figure();
    
    subplot(1,3,1)
    imagesc(Frame);
    colormap("gray")
    axis image
    title('Raw data')
    subplot(1,3,2)
    imagesc(Mask);
    axis image
    title('Mask')
    subplot(1,3,3)
    imshowpair(Frame, Mask);
    axis image
    title("overlay")
    
    SaveName = append(File, filesep, 'SegmentTestframe500.png');
    saveas(Fig, SaveName)
end