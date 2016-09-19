clc;
clear all;
close all

load('results_thr0.5.mat')

%% Initialization
grayLevel = [0.6 0.8];
roundDigit = 0;
numFormat = '%2.0f';

% % allDataSetNames ={'board','david','dollar','faceocc','faceocc2','PETs','singer1','skating1','stone','sylv','trellis'};
% allDataSetNames ={'board','david','dollar','faceocc','faceocc2','PETs','singer1','skating1','stone','sylv'};
% allTrackerNames = {'OAB','MIL','Frag','IVT','APG','MTT','PBGM','PBGS'};

% allDataSetNames = {'sylv','singer1','singer2','skating1','shaking','board','car4','car11','twinnings','faceocc','david','faceocc2','PETS2','PETS1', 'cliffbar','dollar','woman','tiger1','bird1','bird2','trellis70','lemming','stone','liquor'};
% allTrackerNames = {'OAB','MIL','Frag','IVT','APG','MTT','PBGM','PBGS'};

% allDataSetNames = {'board','crossing','david','dollar','faceocc1','faceocc2','singer1','skating1','stone','sylv','trellis','walking2'};
% allTrackerNames = {'OAB','MIL','Frag','IVT','APG','MTT','PBGM','PBGS'};

% allDataSetNames = {'Trellis','board','david','sylv','skating1','singer1','stone','dollar','FaceOcc1'};%,'tiger1'}; %'coke11','faceocc','faceocc2',
allDataSetNames = {'Trellis','board','david','singer1','skating1','FaceOcc2'};
allTrackerNames = {'OAB','MIL','Frag','IVT','APG','MTT','TWSR'};

savePathDir = 'D:\Ali\Framework Matlab\Tracking - matlab framework v0.64\plot';

%% Add Average to Table
allDataSetNames = cat(2,allDataSetNames,{'Average'});

%% Average SuccessRate
avgSuccessRate = 100*mean(allSuccessRates,3)'; 
avgSuccessRate = [avgSuccessRate; mean(avgSuccessRate)];
filename = fullfile(savePathDir,'avgSuccessRate.tex');
mat2LatexTable(avgSuccessRate,allTrackerNames,allDataSetNames,grayLevel,numFormat,roundDigit,filename,'max');

%% My Best SuccessRate
myBestSuccessRate = 100*allSuccessRates(:,:,1);
myBestSuccessRate(end-1:end,:) = 100*max(allSuccessRates(end-1:end,:,:),[],3);
myBestSuccessRate = myBestSuccessRate';
myBestSuccessRate = [myBestSuccessRate; mean(myBestSuccessRate)];
filename = fullfile(savePathDir,'myBestSuccessRate.tex');
mat2LatexTable(myBestSuccessRate,allTrackerNames,allDataSetNames,grayLevel,numFormat,roundDigit,filename,'max');

%% All Best SuccessRates
bestSuccessRate = 100*max(allSuccessRates,[],3);
bestSuccessRate = bestSuccessRate';
bestSuccessRate = [bestSuccessRate; mean(bestSuccessRate)];
filename = fullfile(savePathDir,'bestSuccessRate.tex');
mat2LatexTable(bestSuccessRate,allTrackerNames,allDataSetNames,grayLevel,numFormat,roundDigit,filename,'max');

%% Mean Mean VOC
meanMeanVOC=mean(allVOCmean,3);
meanMeanVOC = meanMeanVOC';
meanMeanVOC = [meanMeanVOC; mean(meanMeanVOC)];
filename = fullfile(savePathDir,'meanMeanVOC.tex');
mat2LatexTable(meanMeanVOC,allTrackerNames,allDataSetNames,grayLevel,'%1.2f',-2,filename,'max');

%% Mean Mean CLE
meanMeanCLE = mean(allCLEmean,3);
meanMeanCLE = meanMeanCLE';
meanMeanCLE = [meanMeanCLE; mean(meanMeanCLE)];
filename = fullfile(savePathDir,'meanMeanCLE.tex');
mat2LatexTable(meanMeanCLE,allTrackerNames,allDataSetNames,grayLevel,numFormat,roundDigit,filename,'min');
