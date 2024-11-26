clc

file.path  = 'C:\Users\Windows 11\OneDrive - KU Leuven\Documents\KU Leuven\PhD\data\Multicolor Project\Testdata_rotational\20241115_AuBPs_2DCal\2DCal_200nm_PS';

file.ext   = '.ome.tif';
info.runMethod = 'run';
info.nChan = 4;
%% 
calib = Core.MPPlaneCalibration(file,info);

calib.retrieveMovies;
calib.calcIndivCal; 
calib.calcCombinedCal;

calib.showCal(1)
calib.offTarget;
calib.save;