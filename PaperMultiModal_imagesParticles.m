clc 
close all

lim = [150 160 250 500 250 200 170 160 150];
for i = 1:8
    filename1 = append('S:\Dual Color\20251015_localisation_error\10ms_expTime\sample_3\calibrated1\calibratedPlane', num2str(i), '.tif');
    Frame1 = Load.Movie.tif.getframes(filename1, 50);

    Fig1 = figure();
    imagesc(Frame1(412:462,33:83))
    particle = Frame1(412:462,33:83);
    clim([median(particle, 'all') 250])

    colormap('hot')
    title(append('Channel 1 - Plane', num2str(i)))
    saveas(Fig1, append('C:\Users\steve\Downloads', filesep, 'Ch1_Plane', num2str(i), '.svg'))
end

for i = 1:8
    filename2 = append('S:\Dual Color\20251015_localisation_error\10ms_expTime\sample_3\calibrated1\calibratedPlane', num2str(i), '.tif');
    Frame2 = Load.Movie.tif.getframes(filename2, 50);

    Fig2 = figure();
    imagesc(fliplr(flipud(Frame2(267:317,82:132))));
    particle = Frame2(267:317,82:132);
    clim([median(particle, 'all') 250])
    colormap('hot')
    title(append('Channel 2 - Plane', num2str(i)))
    saveas(Fig2, append('C:\Users\steve\Downloads', filesep, 'Ch2_Plane', num2str(i), '.svg'))
end