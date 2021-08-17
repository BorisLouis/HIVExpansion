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
path = 'D:\Documents\Unif\PhD\2021-Data\07 - Jul\Aline HIV\DATA Boris\Seperate channels post ExM NUP Lamin\Cell1'; %give empty brackets [], to open file selection
ext = '.tif'; %expected extension of the movie(s);
info.runMethod = 'load'; %'load'or 'run', if load is chosen it will try to load previously calculated data(e.g localized particles)
info.fitMethod = 'phasor'; %'Gauss' or 'phasor'
info.zMethod   = 'Intensity';
membranes = {'lamina','NUP'};
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
HIVData.segmentLamina();
HIVData.showMembrane('lamina');
%% Segmentation of NUP
HIVData.segmentNUP();
HIVData.showMembrane(membranes{2});

%% get membrane Positions
HIVData.getMembranePos(membranes{1});
HIVData.getMembranePos(membranes{2});

%% Show membrane and fit
HIVData.showMembrane(membranes{1},1);
HIVData.showMembrane(membranes{2},1);

%% Plot membrane on top of each other
HIVData.showAllMembranes;

%% Get distance between particles and Lamina
%HIVData.getHIVToMembraneDistance(membranes{1});

%% Get distance between membranes
distance = HIVData.getMembraneToMembraneDistance(membranes{1},membranes{2});

