if ~(exist('batchTest','var') && batchTest == true)
	close all;
    clear;

    rng('default');
    randSeed=0;
    rng(randSeed);
    trackerName = 'IVT';

    addpath(fullfile(pwd,'bin'),'-BEGIN');
    addpath(fullfile(pwd,'trackers',trackerName),'-BEGIN');

    clc;
    settings = set_param();
end

%% Initialization
% Initialize Data Set Information
    [ frameNum, imgNames, imgSize ] = initDataSetInfo(settings.dataSet);

% Read Ground Truth
    [gtCenters, gtCorners,gtInterv] = readGroundTruth(settings.dataSet, frameNum);

% % Make/Load Grayscale Frames
%     frame = loadFrames( settings.dataSet, frameNum, imgNames, imgSize);

% Initialize tracker state
    [ initGeoParam, initAffParam] = initState(gtCenters,gtCorners,settings.objParam);

% Make Patch Indices
    [patchIndex,settings.objParam]= makePatchIndex(settings.objParam);

% % Show First Frame
%     firstFrame = frame(:,:,1);
%     firstCorner = gtCorners(1,:);
%     shape = initFigures(firstFrame,firstCorner,settings.videoParam);

% Show First Frame
    firstFrame = loadFrame( 1, imgNames, settings);
    firstCorner = gtCorners(1,:);

% Initializeing Variables
    tmpl.mean = warpimg(firstFrame, initAffParam, settings.objParam.size);
    tmpl.basis = [];
    tmpl.eigval = [];
    tmpl.numsample = 0;
    tmpl.reseig = 0;

    param = [];
    param.est = initGeoParam;

    param.wimg = tmpl.mean;
    pts = [];

% track the sequence from frame 2 onward
    wimgs = [];
    resultsCorners = zeros(frameNum,8);
    prevParticles = repmat(initGeoParam', 1, settings.pfParam.numsample);
    FPS = zeros(frameNum,1);
    trackerr = zeros(1,frameNum);
    meanerr = zeros(1,frameNum);

reverseStr = [];
fprintf('\n');
prevCLE = 0;
prevVOC = 1;
allCLE = prevCLE;
allVOC = prevVOC;
for f = 1:frameNum

    if (settings.verifyMode == 1) &&  f == settings.dataSet.lastFrame
        verify;
        return;
    end

    tic;
%         currentFrame = frame(:,:,f);
    currentFrame = loadFrame( f, imgNames, settings);

    % Sampling
        if f > 1
            currentParticles = sampling(bestParticle, settings.pfParam);
        else
            currentParticles = prevParticles;
        end


%             sampling like original code
%             if f > 1
%                 currentParticles = IVTsampling(prevParticles, settings.pfParam.numsample, settings.pfParam.affsig, likelihood);
%             else
%                 currentParticles = prevParticles + randn(6,settings.pfParam.numsample).*repmat(settings.pfParam.affsig(:),[1,settings.pfParam.numsample]);
%             end

        param.param = currentParticles;

    % Compute Likelihood
        param = estwarp_condens(currentFrame, tmpl, param, settings.pfParam, settings.omParam);

    % Resampling
        likelihood = param.conf;
        %prevParticles = resampling(currentParticles, likelihood);
        prevParticles = currentParticles;

    % Update Dictionary
        updateDictionary;

    % Compute Corners And Center of the Answer
        bestParticle = param.est;
        bestCorners = paramAff2Corner(paramGeom2Aff(bestParticle), settings.objParam.size);
        resultsCorners(f,:) = bestCorners;

    % Show Tracker
    if settings.videoParam.show
        if f>20
            objDict = tmpl.basis;
%             [objDict, objDictNorm] = normalizeMat(tmpl.basis);
            objDictNorm = ones(1,10);
%             objDictNorm = normalizeMat(tmpl.basis);
            showTracker(currentFrame,bestCorners,f,objDict,objDictNorm,settings)
        end
    end

    % Evaluate Measures
%     if mod(f-1,gtInterv)==0
%         CLE = computeCenterDistError(resultsCorners(f,:), gtCorners(1+((f-1)/gtInterv),:), gtInterv);
%         VOC = computeVOCMeasure(resultsCorners(f,:), gtCorners(1+((f-1)/gtInterv),:), gtInterv);
%         prevCLE = CLE;
%         prevVOC= VOC;
%         allCLE = [allCLE CLE];
%         allVOC = [allVOC VOC];
%     else
%         VOC = prevVOC;
%         CLE = prevCLE;
%     end

    trackTime = toc;
    FPS(f) = 1/trackTime;

    % Display tracker statistics
%     msg = sprintf(' frame = %u \n FPS = %2.2f \n CLE = %2.0f\t\t\tavgCLE = %2.1f \n VOC = %1.2f\t\t\tavgVOC = %1.2f  \n progress = %2.2f%s \n ----------- \n elapsed time = %us',f,FPS(f),CLE,mean(allCLE),VOC,mean(allVOC),f*100/frameNum,'%%',uint32(etime(clock,t0)));
%     msg = sprintf(' frame = %u \n FPS = %2.2f \n CLE = %2.0f \n VOC = %1.3f \n progress = %2.2f%s \n ----------- \n elapsed time = %us',f,FPS(f),CLE,VOC,f*100/frameNum,'%%',uint32(etime(clock,t0)));
    msg = sprintf(' frame = %u \n FPS = %2.2f \n progress = %2.2f%s \n ----------- \n elapsed time = %us',f,FPS(f),f*100/frameNum,'%%',uint32(etime(clock,t0)));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg)-1);
end
fprintf('\n\n');

%% Computing Evaluation Measures
centerDistError = computeCenterDistError(resultsCorners, gtCorners, gtInterv);
VOCMeasure = computeVOCMeasure(resultsCorners, gtCorners, gtInterv);
% centerDistError = mean(allCLE);
% VOCMeasure = mean(allVOC);
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
    rmpath(fullfile(pwd,'IVT'),'-BEGIN');
end
