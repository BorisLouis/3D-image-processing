%%
clc 
close all
clear 

path2ZCal = 'E:\Data\Leuven Data\2018\06-June\27\ZCal - NormObjCorr';
path2SRCal = 'E:\Data\Leuven Data\2018\06-June\27\2DCal - normObjCorrPSFE';

path2File = 'E:\Data\Leuven Data\2018\06-June\27\trackingCal - normObjCorr\TrackingX\FluoBeads200_TrackingX_320ms__1';
path2Cal = 'E:\Data\Leuven Data\2018\06-June\29\2DCal\zStackFLuoBeads_2D3DS3__1';

detectParam.delta = 6;
detectParam.chi2 = 80;
%%
calib = Core.MPCalibration(path2Cal);

%%

MPTrackMov = Core.MPTrackingMovie(path2File,calib.getCal,path2SRCal,path2ZCal);

%% Detection

MPTrackMov.giveInfo
%find candidate
MPTrackMov.findCandidatePos(detectParam);

%fit position
MPTrackMov.SRLocalizeCandidate;

%% Data correction
rot = true;
refPlane = 5;
MPTrackMov.applySRCal(rot,refPlane);
%% e-Z transformation
MPTrackMov.applyZCal;

%% Plane consolidation

MPTrackMov.consolidatePlanes
%% Super resolve
val2Use = 'bestFocus';
MPTrackMov.superResolve(val2Use);
%% plot
frames = 1:100;

MPTrackMov.showCorrLoc(frames);

%% showFrame

MPTrackMov.showFrame(500);

%% tracking
trackParam.euDistXY = 250;
trackParam.euDistZ  = 500;
MPTrackMov.trackParticle(trackParam);
%% plot
MPTrackMov.showTraces;