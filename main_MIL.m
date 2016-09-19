%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is matlab version for MIL tracker, which was proposed by Boris
% Babenko (see http://vision.ucsd.edu/~bbabenko/project_miltrack.shtml for
% more details). I rewrote it for understand it easily. The copyright
% belongs to Boris. Please use it for adacemic purpose. And it is not
% guaranteed to be safe, thus you should use this code at own risk
%
% For bug report, please mail to whluo.china@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~(exist('batchTest','var') && batchTest == true)
	fclose('all');
    close all;
    clear;

    rng('default');
    randSeed=0;
    rng(randSeed);
    trackerName = 'MIL';

    addpath(fullfile(pwd,'bin'),'-BEGIN');
    addpath(fullfile(pwd,'trackers',trackerName),'-BEGIN');

    clc;
    settings = set_param();
end

%% Framework Parameters
% Initialize Data Set Information
    [ frameNum, imgNames, imgSize ] = initDataSetInfo(settings.dataSet);

% Read Ground Truth
    [gtCenters, gtCorners,gtInterv] = readGroundTruth(settings.dataSet, frameNum);

% Initialize tracker state
    [ initGeoParam, initAffParam] = initState(gtCenters,gtCorners,settings.objParam);

% Make Patch Indices
    [patchIndex,settings.objParam]= makePatchIndex(settings.objParam);

% Show First Frame
    firstFrame = loadFrame( 1, imgNames, settings);

resultsCorners = zeros(frameNum,8);


%% MIL Tracker Initialization
inistate=[gtCorners(1,1:2) gtCorners(1,3)-gtCorners(1,1) gtCorners(1,8)-gtCorners(1,2)];
states = zeros(frameNum,4);
states(1,:) = round(inistate);

%set parameters
trparams.posradtrain = 4;
trparams.negnumtrain = 65;
trparams.posmaxtrain = 1000000;

trparams.srchwinsz = 25;
trparams.initpostrainrad = 3;
trparams.initnegnumtrain = 65;
trparams.states = inistate;

%parameters of strong classfier
clfparams.numsel = 50;
clfparams.numfeat = 250;
clfparams.lrate = 0.85;  %learning rate
clfparams.selectedftr = zeros(50,5);

%sep feature parameters
ftrparams.width = 75;
ftrparams.height = 97;
ftrparams.minnumrect = 2;
ftrparams.maxnumrect = 6;

ftrparams.weight = 0;
ftrparams.numrects = 0;
ftrparams.rsums = 0;
ftrparams.maxsum = 0;

ftrs = generateharrftrs(ftrparams);

%% Main Loop
reverseStr = '';
fprintf('\n');
for f=1:frameNum
    tic;
    I0 = loadFrame( f, imgNames, settings);

    c = 0;
    [m n c] = size(I0);
    if c==3
        I1 = I0;
        I0=rgb2gray(I0);
    else
        I1 = repmat(I0,[1,1,3]);%for drawing result
    end

    integral_I0 = integralImage(I0);
    %%%%%%%%initialization
    if f == 1
        [clfle] = iniclfweakle(clfparams);
        possamples = sampleimage(trparams,states(f,:),1,I0);
        negsamples = sampleimage(trparams,states(f,:),0,I0);
        [posftrval,negftrval]=computeharr(possamples,negsamples,ftrs,I0,integral_I0);
        [clf] = updateclfweakle(posftrval,negftrval,clfle,clfparams,0);
        [posres,negres] = getlearnres(clf,posftrval,negftrval);
        selectors = getselectors(posres,negres,clfparams.numsel);
    %%%%%%%%%%run time
    else
        detectx = sampleimage(trparams,states(f-1,:),4,I0);
        prob = classify(detectx,selectors,clf,I0,ftrs,integral_I0);
        [maxprob,idx] = max(prob);
        states(f,:) = detectx(idx,:);
        newim = draw_rect_to_im(I1,states(f,:));
        figure(1)
        imshow(newim);
        hold on
%         imwrite(newim,sprintf('%s_%04d.jpg',[OutputDir 'TrackResult'],f));

        tracknegsamples = sampleimage(trparams,states(f,:),2,I0);
        trackpossamples = sampleimage(trparams,states(f,:),3,I0);
        [trackposftrval,tracknegftrval]=computeharr(trackpossamples,tracknegsamples,ftrs,I0,integral_I0);
        [clf]=updateclfweakle(trackposftrval,tracknegftrval,clf,clfparams,1);
        [trackposres,tracknegres] = getlearnres(clf,trackposftrval,tracknegftrval);
        selectors = getselectors(trackposres,tracknegres,clfparams.numsel);
        clear detectx
        clear prob
        clear idx
        clear tracknegsamples
        clear trackpossamples
        clear trackposftrval
        clear tracknegftrval
        clear trackposres
        clear tracknegres
    end

    trackTime = toc;
	FPS(f) = 1/trackTime;
    msg = sprintf(' frame# = %u \n FPS = %2.2f \n progress = %2.2f%s \n ----------- \n elapsed time = %us',f,FPS(f),f*100/frameNum,'%%',uint32(etime(clock,t0)));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg)-1);
end

%% Find corrners
for i=1:size(states,1)
    resultsCorners(i,:) = paramLoc2Corner([states(i,:),0]);
end
%% Computing Evaluation Measures
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
%
% %% Remove from path
% if ~(exist('batchTest','var') && batchTest == true)
%     rmpath(genpath(fullfile(pwd,'..','toolbox')));
%     rmpath(fullfile(pwd,'bin'));
%     rmpath(fullfile(pwd,trackerName));
% end
