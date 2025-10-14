close all
for i = 1:800
    Ch1(:,:,i) = Load.Movie.tif.getframes('S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol\glycerol_95\glycerol_95_sample_1\calibrated1\calibratedPlane4.tif',i);
    Ch2(:,:,i) = Load.Movie.tif.getframes('S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol\glycerol_95\glycerol_95_sample_1\calibrated2\calibratedPlane4.tif',i);
end

for j = [172]
    fig = figure()
    subplot(1,2,1)
    imagesc(Ch1(180:210, 135:165, j))
    colormap('hot')
    clim([350 550])
    axis square
    hold on

    subplot(1,2,2)
    imagesc(Ch2(165:195, 135:165, j))
    colormap('hot')
    axis square
    clim([350 550])

    filename = append('C:\Users\Windows 11\Downloads', filesep, 'Frame_', num2str(j), '.svg');
    saveas(fig, filename)
end