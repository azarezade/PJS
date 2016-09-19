clear all
fclose('all');
clc;
%% Selected Frame Nums
% frameNum_board = [121 150 200 226 300 401 466 496 631 700];
% frameNum_board = [121 200 300 401 466 696];
% frameNum_crossing = [32 41 58 87 93 109];
% frameNum_david = [50 124 171 229 270 422];
% frameNum_dollar = [50 141 171 211 271 301];
% frameNum_faceocc1 = [87 282 448 536 650 743];
% frameNum_faceocc2 = [173 336 495 611 694 773];
% frameNum_FaceOcc2 = [173 336 495 611 694 773];
% frameNum_singer1 = [42 84 117 228 299 346];
% frameNum_skating1 = [30 73 133 178 265 341];
% frameNum_stone = [141 261 300 406 436 566];
% frameNum_sylv = [25 105 234 344 604 725];
% frameNum_trellis = [42 111 212 307 423 560];
% frameNum_walking2 = [92 166 225 312 442 499];
frameNum_board = [121 300 466 696];
frameNum_david = [50 171 270 422];
frameNum_FaceOcc2 = [173 495 694 773];
frameNum_singer1 = [42 117 299 346];
frameNum_skating1 = [30 133 178 341];
frameNum_stone = [141 261 300 406 436 566];
frameNum_sylv = [25 105 234 344 604 725];
frameNum_trellis = [42 212 423 560];
frameNum_walking2 = [92 166 225 312 442 499];

%% Initialization
% % % % % allDataSetNames = {'Trellis','board','david','singer1','skating1','FaceOcc2'};
allTrackerNames = {'OAB','MIL','Frag','IVT','APG','MTT','TWSR'};
% allTrackerNames = {'OAB','MIL','Frag','IVT','APG','MTT','PBGM','PBGS'};
allTrackerSettings = {'','','','','','','',''};
% allTrackerColors = {'r','g','b','c','m','y','k',[0.6 0.2 0.7],[0.6 0.1 0.3]};
% allTrackerColors = {[255 0 0]/255, [255 120 0]/255, [255 255 0]/255, [0 255 0]/255, [0 255 255]/255, [0 0 255]/255, [255 0 255]/255, [90 10 80]/255};
allTrackerColors = {[170 20 20]/255, [255 120 0]/255, [255 255 0]/255, [0 255 0]/255, [0 255 255]/255, [0 0 255]/255, [255 0 255]/255, [90 10 80]/255};

allRngs = [0 0 0 0 0 0 0 0];

% Initialize dataset info
    settings.dataSet.name = 'trellis';
    sampleFramesNum = frameNum_trellis;
    dim1 = 1;
    dim2 = 4;
    
    settings.dataSet.type = 'Standard';
    settings.dataSet.path = fullfile(pwd, '..','..','Standard DataSet');
%     settings.dataSet.type = 'CVLab'; 
%     settings.dataSet.path = fullfile(pwd, '..','..','..','CVLab');
    
    addpath(genpath(fullfile(pwd,'..','bin')))
    [ frameNum, imgNames, imgSize ] = initDataSetInfo(settings.dataSet);
    
% MAT results directory
    MATresultsDir = 'D:\Ali\Framework Matlab\results-920722-v0.62';

%% Add MAT Files to Patch
addpath(genpath(fullfile(MATresultsDir, 'results_OAB')))
addpath(genpath(fullfile(MATresultsDir, 'results_MIL')))
addpath(genpath(fullfile(MATresultsDir, 'results_Frag')))
addpath(genpath(fullfile(MATresultsDir, 'results_IVT')))
addpath(genpath(fullfile(MATresultsDir, 'results_APG')))
addpath(genpath(fullfile(MATresultsDir, 'results_MTT')))

addpath(genpath(fullfile(MATresultsDir, 'results_PBGM_gamma0.001')))
addpath(genpath(fullfile(MATresultsDir, 'results_PBGS')))

%% Plot
h = figure;
for i=1:length(sampleFramesNum)
    img = loadFrame(sampleFramesNum(i),imgNames,settings);
    subaxis(dim1,dim2,i,'Spacing', 0.01, 'Padding', 0, 'Margin', 0);
    imshow(img); 
    for tn = 1:length(allTrackerNames)
        if tn==5 && (i==4)%(tn==5 && i==4) || (tn==5 && i==6)
%         if (i>=2 && tn==4)||(i>=4 && tn==6) 
%         if (i>=6 && tn==6) 
%         if (i>=2 && tn==4)||(i>=3 && tn==6) 
            continue
        end
        hold on
        load([settings.dataSet.name, '_', allTrackerNames{tn}, '_', 'rng', num2str(allRngs(tn)), allTrackerSettings{tn}, '.mat'],'resultsCorners');
        cornersX = resultsCorners(sampleFramesNum(i),1:2:end);
        cornersY = resultsCorners(sampleFramesNum(i),2:2:end);
        plot([cornersX,cornersX(1)],[cornersY cornersY(1)],'Color',allTrackerColors{tn});
    end
    text(10,10,sprintf('#%04d',sampleFramesNum(i)),'Color','r','FontSize',5);
	hold off
    axis tight
    axis off
end

%% Save 
[s1,s2]=size(img);
s=4*s2/s1;
set(h, 'PaperPosition', [0 0 s 1]); %Position plot at left hand corner with width 5 and height 5.
set(h, 'PaperSize', [s 1]); %Set the paper to have width 5 and height 5.
saveas(h, ['sample' '_' settings.dataSet.name], 'pdf') %Save figure

%% Save Legend

% h2 = figure;
% for tn=1:length(allTrackerNames)
%     hold on
%     load([settings.dataSet.name, '_', allTrackerNames{tn}, '_', 'rng', num2str(allRngs(tn)), allTrackerSettings{tn}, '.mat'], 'centerDistError', 'VOCMeasure', 'successRate');
%     plot(centerDistError,'Color',allTrackerColors{tn});
% end
% legend(allTrackerNames,'orientation','horizontal')
