%%%%%%%%%%%%%%%%%%%%%% PLANE CALIBRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: Folder containing folder containing .ome.tif file that are 2D    %
% calibration.                                                            %
%                                                                         %
% Results: 2D Calibration will be saved (ROI + position of planes in z    %
%                                                                         %
% if the offTarget returned to the command window is >20 nm there might be%
% something wrong with either the analysis,the dataset or the alignment of%
% the camera                                                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

file.path  = 'E:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\2DCal';
file.ext   = '.ome.tif';
info.runMethod = 'run';
info.nChan = 4; %Number of images in 1 channel from 1 camera (mostly 4)
info.method = 'Darkfield Phase'; %Darkfield Phase, Fluorescence,....
%% 
calib = Core.MPPlaneCalibration(file,info);

calib.retrieveMovies;
calib.calcIndivCal; 
calib.calcCombinedCal;

calib.showCal(1)
calib.offTarget;
calib.save;
close all