if ~(exist('batchTest','var') && batchTest == true)
	fclose('all');
    close all;
    clear;

    rng('default');
    randSeed=0;
    rng(randSeed);
    trackerName = 'PJSM';

%     addpath(genpath(fullfile(pwd,'..','toolbox')),'-BEGIN');
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
    
% Initialize tracker state
    [ initGeoParam, initAffParam] = initState(gtCenters,gtCorners,settings.objParam);

% Make Patch Indices
    [patchIndex,settings.objParam]= makePatchIndex(settings.objParam);

% Show First Frame
    firstFrame = loadFrame( 1, imgNames, settings);
    firstCorner = gtCorners(1,:);
    shape = initFigures(firstFrame,firstCorner,settings.videoParam);

% Save Video
if (settings.videoParam.save)
        paramNames=['_gamma',num2str(settings.SC_mFoccus_Param.gamma),'_upd',num2str(settings.updateInterv)];
        videoName= [settings.dataSet.name '_',trackerName,'_rng',num2str(randSeed),paramNames,'.avi'];
        settings.writerMain = VideoWriter(fullfile('results',videoName));
        settings.writerMain.FrameRate = 5;
        open(settings.writerMain);
end

%% Initial Tracking and make objDict
resultsCorners = zeros(frameNum,8);
objDict = [];
objDictNorm = [];

% Initial tracking
 [objDict,objDictNorm,resultsCorners,bestParticle ] = ...
            initTracker(objDict,objDictNorm,resultsCorners,initGeoParam,shape,imgNames,frameNum,patchIndex,settings);
% [objDict,objDictNorm,resultsCorners,bestParticle ] = ...
%             initTracker_old(objDict,objDictNorm,resultsCorners,initGeoParam,shape,imgNames,frameNum,settings);
    
% Make patchDict
    patchDict = reshape(objDict(patchIndex,:),settings.objParam.patchWidth*settings.objParam.patchHeight,[]);
    [patchDict, patchDictNorm] = normalizeMat(patchDict);
    
% make patchBuffers
     for pn=1:settings.objParam.patchNum
        patchBuffer{pn} = patchDict(:,end - (settings.omParam.batchsize*settings.objParam.patchNum) + pn : settings.objParam.patchNum : end);
        patchBufferNorm{pn} = patchDictNorm(end - (settings.omParam.batchsize*settings.objParam.patchNum) + pn : settings.objParam.patchNum : end);
     end    
    
%% Initialize Best Candidate Buffer
% vector of index of group which start from indx 0
    listGroups = int32([(1:settings.groupSize+1:(settings.groupSize+1)*settings.objParam.patchNum*settings.pfParam.numsample)-1, ...
                          settings.pfParam.numsample*settings.objParam.patchNum*(settings.groupSize+1)]);
    
% Initialize group buffer
% Initialize group buffer
    bestBuffer = objDict(:,end-settings.groupSize+1:end);
    groupBuffer = updateBuffer(bestBuffer,settings.groupSize,settings.objParam,patchIndex);

% Initialize candidPatchGroup
%    candidPatchGroup = zeros(settings.objParam.patchWidth*settings.objParam.patchHeight, settings.objParam.patchNum*settings.pfParam.numsample*(settings.groupSize+1));
  

poolInd = repmat(eye(settings.objParam.patchNum), settings.dictObjNum, settings.pfParam.numsample);

%% Tracker Main Loop
alpha = ones(2,settings.objParam.patchNum);
beta = ones(2,settings.objParam.patchNum);

alpha(1,:) = 3;
alpha(2,:) = 7;

beta(1,:) = 7;
beta(2,:) = 3;

wasOccluded = zeros(settings.objParam.patchNum,1);
FPS = zeros(frameNum,1);

initPoint = ones(size(patchDict,2),settings.pfParam.numsample * size(groupBuffer,2));
%initPoint = ones(160,38400);

reverseStr = ''; 
fprintf('\n');
prevCLE = 0;
prevVOC = 1;
allCLE = prevCLE;
allVOC = prevVOC;
for f = 2:frameNum
    
    if (settings.verifyMode == 1) &&  f == settings.dataSet.lastFrame
        verify;
        return;
    end
    
    tic;
    currentFrame = loadFrame( f, imgNames, settings);

    % sampling
     currentParticles = sampling(bestParticle, settings.pfParam);
     
    % Crop candidates
    [ candidates, candidates_norm ] = cropCandidate(currentFrame,currentParticles,settings.objParam);

    % Sparse Coding
	[patchCoef, initPoint] = Coding( candidates, groupBuffer, patchDict, listGroups, patchIndex, initPoint, settings);
    
    % Alignment Pooling
    likelihood = Pooling(patchCoef,settings.objParam, patchDict, candidates, patchIndex, poolInd);

	for p_i = 1:size(currentParticles,2)
        particle = currentParticles(:,p_i);
        corners = paramAff2Corner(paramGeom2Aff(particle), settings.objParam.size);
        for c_i = 1:2:7
            if corners(c_i) < 1 || corners(c_i) > size(currentFrame,2) ...
                    || corners(c_i+1) < 0 || corners(c_i+1) > size(currentFrame,1)
                likelihood(p_i) = -Inf;
                break;
            end
        end
    end
	
    % Find the Best Candidate
    [~, bestParticleInd] = max(likelihood);
    bestParticle = currentParticles(:, bestParticleInd);     
    
    % Dictionary Update
    bestCandidDataNrom=candidates_norm(bestParticleInd);
    bestCandidData=candidates(:,bestParticleInd);
    [ patchDict,patchDictNorm,objDict,objDictNorm, alpha, beta, wasOccluded] ...
       = updateDict(bestCandidData,bestCandidDataNrom,patchDict,patchDictNorm,objDict,objDictNorm,patchIndex,f,settings, alpha, beta, wasOccluded);
    
    % Update Buffer
    if settings.groupSize > 0
        bestBuffer(:,1) = [];
        bestBuffer = [bestBuffer, bestCandidData];
    end
    groupBuffer = updateBuffer(bestBuffer,settings.groupSize,settings.objParam,patchIndex);

    % Compute corners of the answer
    bestCorners = paramAff2Corner(paramGeom2Aff(bestParticle), settings.objParam.size);
    resultsCorners(f,:) = bestCorners;

   % Show Results   
    showTracker(currentFrame,bestCorners,f,objDict,objDictNorm,settings)
    
    % Save Results
    settings = saveVideo(settings,f,frameNum);

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
    MATFileName = [settings.dataSet.name,'_',trackerName,'_rng',num2str(randSeed),...
        '_gamma',num2str(settings.SC_mFoccus_Param.gamma),'_upd',num2str(settings.updateInterv),'.mat'];
    save(fullfile('results',MATFileName),'centerDistError', 'VOCMeasure', 'successRate','FPS'...
         ,'resultsCorners','gtCorners','gtInterv','settings','randSeed');
end

%% Plot results

%% Remove from path
if ~(exist('batchTest','var') && batchTest == true)
    rmpath(genpath(fullfile(pwd,'..','toolbox')));
    rmpath(fullfile(pwd,'bin'));
    rmpath(fullfile(pwd,trackerName));
end