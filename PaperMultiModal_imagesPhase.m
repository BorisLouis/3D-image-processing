clc 
close all

lim = [0 6 10 13];
for i = lim
    filename1 = append('E:\Multimodal tracking\20250724\alldata\50mgmL_sample2_min', num2str(i), '__1\50mgmL_sample2_min', num2str(i), '__1_MMStack_Pos0.ome');
    Frame1 = Load.Movie.ome.load(.getframes(filename1, 50);

    Fig1 = figure();
    imagesc(Frame1(25:end-25, 25:end-25));
    colormap('gray')
    title(append('minute', num2str(i)));
end

for i = 1:8
    filename2 = append('S:\Dual Color\20251015_localisation_error\10ms_expTime\sample_3\calibrated1\calibratedPlane4.tif');
    Frame2 = Load.Movie.tif.getframes(filename2, 50);

    Fig2 = figure();
    imagesc(fliplr(flipud(Frame2(267:317,82:132))));
    particle = Frame2(267:317,82:132);
    clim([median(particle, 'all') 250])
    colormap('hot')
    title(append('Channel 2 - Plane', num2str(i)))
    saveas(Fig2, append('C:\Users\steve\Downloads', filesep, 'Ch2_Plane', num2str(i), '.svg'))
end