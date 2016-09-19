if ~(exist('batchTest','var') && batchTest == true)
    clear;
    close all;
    clc;

    rng('default');
    randSeed=0;
    rng(randSeed);
    trackerName = 'MTT';

    addpath(fullfile(pwd,'bin'),'-BEGIN');
    addpath(fullfile(pwd,'trackers',trackerName),'-BEGIN');

    clc;
    settings = set_param();

end

% video_name = 'car11';
% video_path = fullfile('.\data\',video_name);
% m_start_frame = 1;  %starting frame number
% nframes		= 393; %393;	 %number of frames to be tracked
% Imgext		= 'png';				%image format
% numzeros	= 4;	%number of digits for the frame index
% all_images	= cell(nframes,1);
% nz			= strcat('%0',num2str(numzeros),'d'); %number of zeros in the name of image
% for t=1:nframes
%     image_no	= m_start_frame + (t-1);
%     fid			= sprintf(nz, image_no);
%     all_images{t}	= strcat(video_path,'\',fid,'.',Imgext);
% end
 % Initialize Data Set Information
 [ frameNum, imgNames, imgSize ] = initDataSetInfo(settings.dataSet);
 nframes  = frameNum;

 for t=1:nframes

	if isequal(settings.dataSet.type,'CVLab')
        all_images{t}= fullfile(settings.dataSet.path,settings.dataSet.name,'img',imgNames{t});
    elseif isequal(settings.dataSet.type,'Standard')
        all_images{t}= fullfile(settings.dataSet.path,settings.dataSet.name,'imgs',imgNames{t});
    end
    imgNames{t} = all_images{t};
end

% Read Ground Truth
    [gtCenters, gtCorners,gtInterv] = readGroundTruth(settings.dataSet, frameNum);

%% initialize bounding box
% m_boundingbox = [gtCorners(1,1),gtCorners(1,2),gtCorners(1,3)-gtCorners(1,1),gtCorners(1,6)-gtCorners(1,2)];  % [left-top-x, left-top-y, width, height];

% init_pos	= SelectTarget(all_images{1});  % automatically get bounding box
% init_pos =  [p1 p2 p3];
% 			  p1-------------------p3
% 				\					\
% 				 \       target      \
% 				  \                   \
% 				  p2-------------------\
init_pos = [gtCorners(1,2)  gtCorners(1,8)  gtCorners(1,4);
            gtCorners(1,1)  gtCorners(1,7)  gtCorners(1,3)];

opt.init_pos = double(init_pos);  %  initialization bounding box
%
% width = m_boundingbox(3);
% height = m_boundingbox(4);

%% 	set object size including height and width based on the initialization
%if min( 0.5*[height width]) < 25
 %   sz_T = 1.0 * [height width];
 %   if height > 80
 %       sz_T =  [ 0.5 *height width];
 %   end
%else
%    sz_T = 0.5 * [height width];
%end
%sz_T = ceil(sz_T);
%if min(sz_T>32)
%    sz_T = [32 32];
%end
sz_T = settings.objParam.size;
%    begin = 1;
%    firstFrame = frame(:,:,begin);
%    initGraphics;
resultsCorners = zeros(frameNum,8);
%% MTT tracking Parameters
opt.n_sample = settings.pfParam.numsample;		% number of particles   400
opt.sz_T= sz_T;         % object size
opt.tracker_type = 'L21';  opt.lambda = 0.01; % three different trackers: L21, L11, L01(denote L\infinity 1);
% opt.tracker_type = 'L11';  opt.lambda = 0.005;
% opt.tracker_type = 'L01';  opt.lambda = 0.2;
opt.eta  = 0.01;
opt.obj_fun_th = 1e-3;
opt.iter_maxi = 100; % lambda, eta,obj_fun_th, and iter_maxi are parameters for Accelerated Proximal Gradient (APG) Optimization. Please refer to our paper for details.
opt.rel_std_afnv =  [settings.pfParam.affsig(3),settings.pfParam.affsig(5),settings.pfParam.affsig(4),settings.pfParam.affsig(6),...
    settings.pfParam.affsig(1),settings.pfParam.affsig(2)];
%0.005,0.0005,0.0005,0.005,4,4]; % affine parameters for particle sampling

opt.m_theta = 0.6;  % [0 1] decide object template update
opt.show_optimization = false; % show optimization results to help tue eta and lambda for APG optimization.
opt.show_time = true; % show optimization speed

%% Run MTT tracking. To get better results for different videos, we can change sz_T, rel_std_afnv,  and m_theta.
% m_track: tracking result;
% m_speed: the time (second) per frame
%[tracking_res,m_speed] = MTT_Tracking(all_images,opt);
tic;
[tracking_res,m_speed] = MTT_Tracking(all_images,opt);
FPS=toc;
FPS=FPS/frameNum;
%% Save tracking results
%all_results_path = '.\MTT_Results\';
% if ~exist([all_results_path video_name])
%     mkdir([all_results_path video_name]);
% end
all_rect = [];
prevCLE = 0;
prevVOC = 1;
for t = 1:nframes
     img_color	= imread(imgNames{t});
     img_color	= double(img_color);
     imshow(uint8(img_color));
     text(5,10,num2str(t),'FontSize',18,'Color','r');
     color = [1 0 0];
    map_afnv	= tracking_res(:,t)';
    rect=drawAffine(map_afnv, sz_T, color, 2);
     resultsCorners(t,:) = [rect(2,1), rect(1,1), rect(2,1),rect(1,2),rect(2,3),rect(1,2),rect(2,3),rect(1,1)];
%    all_rect =[all_rect; rect(2,1) rect(1,1) rect(2,3)-rect(2,1) rect(1,2)-rect(1,1)];

%     s_res	= all_images{t}(1:end-4);
%     s_res	= fliplr(strtok(fliplr(s_res),'/'));
%     s_res	= fliplr(strtok(fliplr(s_res),'\'));
%     s_res	= [s_res '_MTT.png'];
%     f = getframe(gcf);
%     imwrite(uint8(f.cdata), [all_results_path video_name '\' s_res]);

%     % Compute current CLE
%     if mod(t-1,gtInterv)==0
%         CLE = computeCenterDistError(resultsCorners(t,:), gtCorners(t,:), gtInterv);
%         VOC = computeVOCMeasure(resultsCorners(t,:), gtCorners(t,:), gtInterv);
%         prevCLE = CLE;
%         prevVOC= VOC;
%     else
%         VOC = prevVOC;
%         CLE = prevCLE;
%     end
%     fprintf('CLE = %3.1f \n VOC = %3.1f \n',CLE,VOC);
end

centerDistError = computeCenterDistError(resultsCorners, gtCorners, gtInterv);
VOCMeasure = computeVOCMeasure(resultsCorners, gtCorners, gtInterv);
successRate = computeSuccessRate(VOCMeasure, 0.6);

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
    rmpath(fullfile(pwd,'MTT'),'-BEGIN');
    rmpath(genpath(fullfile(pwd,'MTT','MTT_Toolbox')),'-BEGIN');
end
