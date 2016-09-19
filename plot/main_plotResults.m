clear all
close all
clc;

%% Initialization
param.threshold = 0.5;

param.CLE.showAll = 0;
param.CLE.showAvg = 1;
param.CLE.showVar = 0;

param.CLE.saveAll = 0;
param.CLE.saveAvg = 1;
param.CLE.saveVar = 0;

param.VOC.showAll = 0;
param.VOC.showAvg = 1;

param.VOC.saveAll = 0;
param.VOC.saveAvg = 1;

param.resultsDir = 'D:\Ali\Framework Matlab\results-920722-v0.62';
% % param.resultsDir = '/home/ali/Documents/Thesis/codes/Matlab Framework/results-920510-v0.59';
% param.resultsDir = 'D:\Ali\Framework Matlab\Tracking - matlab framework v0.64\results';
addpath(genpath(param.resultsDir));

%% Add general and MAT files to path
addpath(genpath(fullfile(pwd,'..','bin')))
 
addpath(genpath(fullfile(param.resultsDir, 'results_OAB')))
addpath(genpath(fullfile(param.resultsDir, 'results_MIL')))
addpath(genpath(fullfile(param.resultsDir,'results_Frag')))
addpath(genpath(fullfile(param.resultsDir, 'results_IVT')))
addpath(genpath(fullfile(param.resultsDir, 'results_APG')))
addpath(genpath(fullfile(param.resultsDir, 'results_MTT')))

addpath(genpath(fullfile(param.resultsDir, 'results_PJS-M_gamma0.001')))
addpath(genpath(fullfile(param.resultsDir, 'results_PJS-S')))

%% Create folder of figures
if ~exist(fullfile(param.resultsDir,'Figures CLE'),'dir')
    mkdir(fullfile(param.resultsDir,'Figures CLE'))
end
if ~exist(fullfile(param.resultsDir,'Figures VOC'),'dir')
    mkdir(fullfile(param.resultsDir,'Figures VOC'))
end

%% Set Tracker and Dataset Names
% allDataSetNames ={'sylv','singer1','singer2','skating1','shaking','board','car4','car11','twinnings','faceocc','david_indoor','faceocc2','OneLeaveShopReenter2cor','PETS01D1Human1', 'cliffbar','dollar','woman','tiger1','bird1','bird2','trellis70','lemming','stone','liquor'};

% allDataSetNames = {'board','crossing','david','dollar','faceocc1','faceocc2','singer1','skating1','stone','sylv','trellis','walking2'};
% allTrackerNames = {'OAB','MIL','Frag','IVT','APG','MTT','PBGM','PBGS'};
% allTrackerSettings = {'','','','','','','_gamma0.001',''};

% allDataSetNames = {'Trellis','board','david','sylv','skating1','singer1','stone','dollar','FaceOcc1'};%,'tiger1'}; %'coke11','faceocc','faceocc2',
allDataSetNames = {'Trellis','board','david','singer1','skating1','FaceOcc2'};     % 
allTrackerNames = {'OAB','MIL','Frag','IVT','APG','MTT','TWSR'};
allTrackerSettings = {'','','','','','',''};
% allTrackerColors = {'r','g','b','c','m',[1 0.6 0],'k',[0.6 0.2 0.7],[0.6 0.1 0.3]};
% allTrackerColors = {[170 20 20]/255, [255 120 0]/255, [255 255 0]/255, [0 255 0]/255, [0 255 255]/255, [0 0 255]/255, [255 0 255]/255, [90 10 80]/255};
allTrackerColors = {[170 20 20]/255, [255 120 0]/255, [255 255 0]/255, [0 255 0]/255, [0 255 255]/255, [0 0 255]/255, [255 0 255]/255, [90 10 80]/255};

rngRange = 0:4;
% rngRange = 0;
legendNames = allTrackerNames;

%% Plot and Save
allSuccessRates=zeros(length(allTrackerNames),length(allDataSetNames),length(rngRange));
for dn=1:length(allDataSetNames)
    % ======= Plot all CLE and Find SuccessRate ======= 
    allCLE=[];
    [allCLE,allSuccessRates,avgSuccessRates,gtInterv,frameNum] = plotAllCLE(allCLE,allSuccessRates,allDataSetNames,allTrackerNames,allTrackerSettings,allTrackerColors,legendNames,dn,rngRange,param);
    allCLEmean(:,dn,:) = mean(allCLE,2);
    allCLEstd(:,dn,:) = std(allCLE,[],2);
    
    % ======= Plot avg CLE =======
    [allCLEAvg,allCLEvar] = plotAvgCLE(allCLE,allDataSetNames,allTrackerNames,allTrackerColors,legendNames,dn,gtInterv,frameNum,param);
    
    % ======= Plot CLE var =======
    plotCLEvar(allCLE,allCLEAvg,allCLEvar,allDataSetNames,allTrackerNames,allTrackerColors,dn,param);
    
    % ======= Plot all VOC =======
    allVOC=[];
    [allVOC] = plotAllVOC(allVOC,allDataSetNames,allTrackerNames,allTrackerSettings,allTrackerColors,legendNames,dn,rngRange,param);
    
    % ======= Plot avg VOC =======
    [allVOCAvg,allVOCvar] = plotAvgVOC(allVOC,allDataSetNames,allTrackerNames,allTrackerColors,legendNames,dn,gtInterv,frameNum,param);
%     VOCmean = squeeze(mean(allVOC,2));
%     VOCstd = squeeze(std(allVOC,[],2));
    allVOCmean(:,dn,:) = mean(allVOC,2);
    allVOCstd(:,dn,:) = std(allVOC,[],2);

    % ====== Plot ROC =======
    plotBestROC(allVOCmean(:,dn,:),allDataSetNames,allTrackerNames,allTrackerSettings,allTrackerColors,dn,param);
    plotAvgROC(allVOCmean(:,dn,:),allDataSetNames,allTrackerNames,allTrackerSettings,allTrackerColors,dn,rngRange,param);
end

%% Save
save(['results','_','thr',num2str(param.threshold),'.mat'])

%% Remove path
% rmpath(genpath(param.resultsDir))
% close all