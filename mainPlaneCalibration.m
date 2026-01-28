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

file.path  = 'C:\Users\steve\data_no_onedrive\Data Hannah\202512117\2D-cal';
file.ext   = '.ome.tif';
info.runMethod = 'run';
info.nChan = 4; %Number of images in 1 channel from 1 camera (mostly 4)
info.Method = 'Phase'; %Darkfield Phase, Fluorescence,....
info.Modality = 'Fluorescence'; %Darkfield phase
%% 
calib = Core.MPPlaneCalibration(file,info);

calib.retrieveMovies;
calib.calcIndivCal; 
calib.calcCombinedCal;
close all
calib.showCal(1)
calib.offTarget;
calib.save;
