if ~(exist('batchTest','var') && batchTest == true)
	close all;
    clear;
    clc;

    rng('default');
    randSeed=3;
    rng(randSeed);
    trackerName = 'APG';

    addpath(fullfile(pwd,'bin'),'-BEGIN');
    addpath(fullfile(pwd,'trackers',trackerName),'-BEGIN');

    clc;
    settings = set_param();
end

%% Initialization
% Initialize Data Set Information
    [ frameNum, imgNames, imgSize ] = initDataSetInfo(settings.dataSet);
disp(frameNum);
% Read Ground Truth
    [gtCenters, gtCorners,gtInterv] = readGroundTruth(settings.dataSet, frameNum);

% Make/Load Grayscale Frames
%     frame = loadFrames( settings.dataSet, frameNum, imgNames, imgSize);
%    frame=frame(:,:,12:end);
%     frameNum=frameNum-12;
%     gtCenters=gtCenters((12/gtInterv):end,:);
%     gtCorners=gtCorners((12/gtInterv):end,:);
% Initialize tracker state
    [ initGeoParam, initAffParam] = initState(gtCenters,gtCorners,settings.objParam);

%times  = 1; %operate times; to avoid overwriting previous saved tracking result in the .mat format
%title = 'car4';
%res_path='results\';

%% parameter setting for each sequence
%switch title
%    case 'car4'
%        fprefix		= '.\car4\';
%        fext		= 'jpg';    %Image format of the sequence
%        numzeros	= 4;		%number of digits for the frame index
%        start_frame = 12;		% first frame index to be tracked
%        nframes		= 600;		% number of frames to be tracked
        %Initialization for the first frame.
        %Each column is a point indicating a corner of the target in the first image.
        %The 1st row is the y coordinate and the 2nd row is for x.
        %Let [p1 p2 p3] be the three points, they are used to determine the affine parameters of the target, as following
        %    p1(65,55)-----------p3(170,53)
        %         | 				|
        %         |     target      |
        %         | 				|
        %   p2(64,140)--------------
%        init_pos= [55,140,53;
%                   65,64,170];
%              % size of template
%end
init_pos = [gtCorners(1,2),gtCorners(1,8),gtCorners(1,4);gtCorners(1,1),gtCorners(1,7),gtCorners(1,3)];
disp(init_pos);
%prepare the file name for each image
%s_frames = cell(nframes,1);
%nz	= strcat('%0',num2str(numzeros),'d'); %number of zeros in the name of image
%for t=1:nframes
%    image_no	= start_frame + (t-1);
%    id=sprintf(nz,image_no);
%    s_frames{t} = strcat(fprefix,id,'.',fext);
%end

%prepare the path for saving tracking results
%res_path=[res_path title '\'];
%if ~exist(res_path,'dir')
%    mkdir(res_path);
%end
%% parameters setting for tracking
%sz_T =[12,15];
sz_T = settings.objParam.size;
para.lambda = [0.2,0.001,10]; % lambda 1, lambda 2 for a_T and a_I respectively, lambda 3 for the L2 norm parameter
% set para.lambda = [a,a,0]; then this the old model
para.angle_threshold = 40;
para.Lip	= 8;
para.Maxit	= 5;
para.nT		= settings.dictObjNum;%number of templates for the sparse representation
para.rel_std_afnv = [settings.pfParam.affsig(3),settings.pfParam.affsig(5),settings.pfParam.affsig(4),settings.pfParam.affsig(6),...
    settings.pfParam.affsig(1),settings.pfParam.affsig(2)];%diviation of the sampling of particle filter
para.n_sample	= settings.pfParam.numsample;		%number of particles
para.sz_T		= sz_T;
%to be changed
para.init_pos	= init_pos;
para.bDebug		= 0;		%debugging indicator
bShowSaveImage	= 1;       %indicator for result image show and save after tracking finished
%para.s_debug_path = res_path;

%% main function for tracking
[tracking_res,output]  = L1TrackingBPR_APGup(imgNames, para,frameNum, settings, gtCenters, gtCorners, gtInterv);

disp(['fps: ' num2str(frameNum/sum(output.time))]);
FPS=frameNum/sum(output.time);
resultsCorners = zeros(frameNum,8);
results=round(aff2image(tracking_res, sz_T));
for i=1:size(tracking_res,2)
rect= results(:,i);
inp	= reshape(rect,2,4);

topleft_r = inp(1,1);
topleft_c = inp(2,1);
botleft_r = inp(1,2);
botleft_c = inp(2,2);
topright_r = inp(1,3);
topright_c = inp(2,3);
botright_r = inp(1,4);
botright_c = inp(2,4);
resultsCorners(i,:)=[topleft_r,topleft_c,topright_r,topright_c,botright_r,botright_c,botleft_r,botleft_c];
end
res = resultsCorners;
res(:,1)=resultsCorners(:,2);
res(:,2)=resultsCorners(:,1);
res(:,3)=resultsCorners(:,4);
res(:,4)=resultsCorners(:,3);
res(:,5)=resultsCorners(:,6);
res(:,6)=resultsCorners(:,5);
res(:,7)=resultsCorners(:,8);
res(:,8)=resultsCorners(:,7);
resultsCorners = res;

%%
centerDistError = computeCenterDistError(resultsCorners, gtCorners, gtInterv);
VOCMeasure = computeVOCMeasure(resultsCorners, gtCorners, gtInterv);
successRate = computeSuccessRate(VOCMeasure,settings.VOCThreshold);

%% Save MAT file
if ~(exist('batchTest','var') && batchTest == true)
    MATFileName = [settings.dataSet.name,'_',trackerName,'_rng',num2str(randSeed) '.mat'];
    save(fullfile('results',MATFileName),'centerDistError', 'VOCMeasure', 'successRate','FPS'...
         ,'resultsCorners','gtCorners','gtInterv','settings','randSeed');
end


%% Plot results

%% Remove from path
if ~(exist('batchTest','var') && batchTest == true)
    rmpath(genpath(fullfile(pwd,'toolbox')),'-BEGIN');
    rmpath(fullfile(pwd,'bin'),'-BEGIN');
    rmpath(fullfile(pwd,'APG'),'-BEGIN');
end
%% Output tracking results

%save([res_path title '_L1_APG_' num2str(times) '.mat'], 'tracking_res','sz_T','output');

%if ~para.bDebug&bShowSaveImage
%    for t = 1:nframes
%        img_color	= imread(s_frames{t});
%        img_color	= double(img_color);
%        imshow(uint8(img_color));
%        text(5,10,num2str(t+start_frame),'FontSize',18,'Color','r');
%        color = [1 0 0];
%        map_afnv	= tracking_res(:,t)';
%        drawAffine(map_afnv, sz_T, color, 2);%draw tracking result on the figure
%        drawnow
        %save tracking result image
%        s_res	= s_frames{t}(1:end-4);

%s_res	= fliplr(strtok(fliplr(s_res),'/'));
 %       s_res	= fliplr(strtok(fliplr(s_res),'\'));
 %       s_res	= [res_path s_res '_L1_APG.jpg'];
 %       saveas(gcf,s_res)
 %   end
%end
