clear
close all
clc

file.path  = 'C:\Users\steve\OneDrive\Documenten\TestData Indra\20250613_test_camera_cal';
file.ext   = '.his';
info.runMethod = 'run';
info.nChan = 1; %Number of images in 1 channel from 1 camera
info.method = 'Fluorescence'; %Darkfield Phase, Fluorescence,....
%% 
calib = Core.Channel2DCalibration(file,info);

calib.retrieveMovies;
calib.calcIndivCal; 
calib.calcCombinedCal;
close all