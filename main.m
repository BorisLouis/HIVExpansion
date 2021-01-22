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
path = 'D:\Documents\Unif\PhD\2021-Data\19\LaminHIV'; %give empty brackets [], to open file selection
ext = '.tif'; %expected extension of the movie(s);
info.runMethod = 'run'; %'load'or 'run', if load is chosen it will try to load previously calculated data(e.g localized particles)
info.fitMethod = 'phasor'; %'Gauss' or 'phasor'
info.zMethod   = 'Intensity';
%if it exist. run will always re-run the analysis and erase previous data.

%% Data Loading

HIVData = Core.HIVLocMovie(ext,info,path);

HIVData.getExtraInfo();

%% Particle detection
%rough detection of particle
detectParam.delta = 6; % ROI around the molecule for testing
detectParam.chi2  = 60;% Threshold for detections ([24-80],24 for single molecules, up to 80 for brighter objects
HIVData.findCandidatePos(detectParam);

HIVData.showCandidate(1,5);

%% Localization
%Fitting for accurate localization of the detected particle
HIVData.SRLocalizeCandidate();

%% consolidation of particle position
%check that particle are detected in more than one z-slice
HIVData.consolidatePlanes();

%% Super-resolve in 3D
HIVData.superResolve();