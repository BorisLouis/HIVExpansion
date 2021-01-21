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
path = 'D:\Documents\Unif\PhD\2021-Data\19\Lamina-NUP\cell1'; %give empty brackets [], to open file selection
ext = '.tif'; %expected extension of the movie(s);
info = []; %placeholder for future information potentially needed


%% Data Loading

HIVData = Core.HIVMovie(ext,info,path);

HIVData.getExtraInfo()

%% Data Processing





