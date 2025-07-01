clear
close all
clc

file.path  = 'F:\Indra\260625_cam1\merged';
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