%% main
%This code aims to be used by the end-user for analysis of HIV particles in
%cells as well as segmentation of various cell structures. This code will
%mostly call the main functions needed while most of the work will occur in
%the background. This is to make this code as understandable and short as
%possible
clear
close all
clc

%% User Input
path = 'D:\Documents\Unif\PhD\2021-Data\08 - August\Aline HIV\Cell'; %give empty brackets [], to open file selection
ext = '.tif'; %expected extension of the movie(s);
info.runMethod = 'load'; %'load'or 'run', if load is chosen it will try to load previously calculated data(e.g localized particles)
info.fitMethod = 'phasor'; %'Gauss' or 'phasor'
info.zMethod   = 'Intensity';
%if it exist. run will always re-run the analysis and erase previous data.

%% Data Loading

HIVData = Core.HIVCellMovie(ext,info,path);

HIVData.getExtraInfo();

%% Particle detection
%rough detection of particle
detectParam.delta = 6; % ROI around the molecule for testing
detectParam.chi2  = 60;% Threshold for detections ([24-80],24 for single molecules, up to 80 for brighter objects
HIVData.findCandidatePos(detectParam);
if ~isempty(HIVData.candidatePos)
    HIVData.showCandidate(1,5);
end
%% Localization
%Fitting for accurate localization of the detected particle
HIVData.SRLocalizeCandidate();

%% consolidation of particle position
%check that particle are detected in more than one z-slice
HIVData.consolidatePlanes();

%% Super-resolve in 3D
HIVData.superResolve();


%% Segmentation of Lamina
membrane = 'lamina';
sensitivity = 0.4;
HIVData.segmentLamina(sensitivity,true);
HIVData.showMembrane(membrane);
%% Segmentation of NUP
sensitivity = 0.4;%put higher when signal is worst (e.g. expansion images)
HIVData.segmentNUP(sensitivity,true);
membrane = 'NUP';
HIVData.showMembrane(membrane);
%% Segmentation of Lipid
membrane = 'lipid';
%HIVData.segmentLipid();
%HIVData.showMembrane(membrane);
%% Segmentation of red lipid
HIVData.segmentRedLipid();
HIVData.showMembrane('lipid',1);

%% fit lipid
HIVData.getMembranePos('lipid');
HIVData.showMembrane('lipid');
%% get membrane position 
membrane = 'lamina';
HIVData.getMembranePos(membrane);
HIVData.showMembrane(membrane,1);
% get second membrane position
membrane = 'NUP';
HIVData.getMembranePos(membrane);
HIVData.showMembrane(membrane,2);

%%
HIVData.showAllMembranes;

%% Get Distance between membranes
distance = HIVData.getMembraneToMembraneDistance('lamina','NUP');

%% Get distance between particles and membrane

a = HIVData.getHIVToMembraneDistance('NUP');




