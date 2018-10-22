%%
clc 
close all
clear 

path2ZCal = 'E:\Data\Leuven Data\2018\06-June\27\ZCal - NormObjCorr';
path2SRCal = 'E:\Data\Leuven Data\2018\06-June\27\2DCal - normObjCorrPSFE';

path2File = 'E:\Data\Leuven Data\2018\06-June\27\ZCal - NormObjCorr\zStackFluoBeads200_PIC1_270618__3';
path2Cal = 'E:\Data\Leuven Data\2018\06-June\27\2DCal - normObjCorrPSFE\zStackFluoBeads200_S3_270618__2';

detectParam.delta = 6;
detectParam.chi2 = 80;
info.runMethod = 'load';
info.type = 'normal';
%%
calib = Core.MPCalibration(path2Cal);

%%

MPTrackMov = Core.MPTrackingMovie(path2File,calib.getCal,info,path2SRCal,path2ZCal);

%% Detection

MPTrackMov.giveInfo
%find candidate
MPTrackMov.findCandidatePos(detectParam);

%fit position
MPTrackMov.SRLocalizeCandidate;

%% Data correction
rot = true;
refPlane = 3;
MPTrackMov.applySRCal(rot,refPlane);
%% e-Z transformation
MPTrackMov.applyZCal;

%% Plane consolidation

MPTrackMov.consolidatePlanes
%% Super resolve
val2Use = 'bestFocus';
MPTrackMov.superResolve(val2Use);
%% plot
frames = 1:500;

MPTrackMov.showCorrLoc();

%% showFrame

%MPTrackMov.showFrame(80,5);
%MPTrackMov.showParticle;

%% tracking
trackParam.euDistXY = 400;
trackParam.euDistZ  = 400;
MPTrackMov.trackParticle(trackParam);
traces = MPTrackMov.getTraces;
%% plot
%MPTrackMov.showTraces;

%% eval accuracy
MPTrackMov.evalAccuracy

%% get RMSD

[RMSD,D] = MPTrackMov.getRMSD('3D');

%% Intensity and SNR
[int,SNR] = MPTrackMov.getAvgIntensity;

%% getTraces 3D
frame = 1:85;
idx2Trace = 2;
ROIradius = 12;
frameRate = 10;
scaleBar  = 1; %in nm 
MPTrackMov.getTracesMovie(frame,idx2Trace,ROIradius,frameRate,scaleBar);
%MPTrackMov.getTraces3DMovie(frame,idx2Trace,ROIradius,frameRate);
%MPTrackMov.getPartMovie(frame,idx2Trace,ROIradius,frameRate);
