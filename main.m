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
ext = '.tif'; %expected extension of the movie(s);
info = []; %placeholder for future information potentially needed


%% Data Loading

HIVData = Core.HIVMovie(info);