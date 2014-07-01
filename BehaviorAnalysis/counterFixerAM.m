function [ fixedFiles ] = counterFixerAM( dataPath )
%Given fName = file directory for running wheel locomotion data, go through
%the entire folder and, if file is a .mat file, search for ds.tCounterOn
%and ds.tCounterOff and insert a itiCounter variable that is
%itiCounter=tCounterOn-tCounterOff.

%%   Detailed explanation goes here
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
% Data Path Formatting
disp('Loading data from MWorks Data folder...')
dirListing = dir(dataPath);
fileNames= {dirListing.name};
nFiles= length(fileNames);


end

