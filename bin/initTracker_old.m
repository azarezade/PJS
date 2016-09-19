function [objDict,objDictNorm,resultsCorners,bestParticle ] = ...
            initTracker_old(objDict,objDictNorm,resultsCorners,initGeoParam,shape,imgNames,frameNum,settings)
objSize = settings.objParam.size;
begin = 1;
% simple tracking and collecting tracking results as templates
%% simple tracking
    prevParticles = repmat(initGeoParam',1,settings.pfParam.numsample);
    for f = begin : begin+settings.dictObjNum-1  
%         tic;
        currentFrame = loadFrame( f, imgNames, settings);
        % Sampling
            if f > 1
                %currentParticles = sampling(prevParticles, settings.pfParam.numsample, settings.pfParam.affsig);
%                 currentParticles = sampling(repmat(prevParticles(:,bestParticleInd), 1, settings.pfParam.numsample), settings.pfParam.numsample, settings.pfParam.affsig);
                currentParticles = sampling(bestParticle, settings.pfParam);
            else
                currentParticles = prevParticles;
            end

        candidates = warpimg(currentFrame, paramGeom2Aff(currentParticles), objSize);
        candi_data = reshape(candidates, objSize(1)*objSize(2), settings.pfParam.numsample); 
        candi_data = candi_data.*(candi_data>0);  

        if f > 1
            % use knn function of the vlfeat open source library
            candidate_kdTree = vl_kdtreebuild(candi_data);   
            [bestParticleInd, distances] = vl_kdtreequery( candidate_kdTree, candi_data, objDict(:,end), 'NumNeighbors', 1);        
        else
            bestParticleInd = 1;
        end
        bestParticle = currentParticles(:,bestParticleInd);
        %likelihood = zeros(1,settings.pfParam.numsample);
        %likelihood(bestParticleInd) = 1;

        objDict = [objDict, candi_data(:,bestParticleInd)];
        objDictNorm = [objDictNorm, 1];

        % Compute Corners And Centers of the Answer
            bestCorners = paramAff2Corner(paramGeom2Aff(bestParticle), objSize);
            resultsCorners(f,:) = bestCorners;


        % Show Results
%         showTracker(currentFrame,bestCorners,f,objDict,objDictNorm,settings)

        % Resampling
            %%%%%%%%%%%%%%%% bayad resampling ezafe konim????? %%%%%%%%%%%
            %prevParticles = resampling(currentParticles, likelihood);
            prevParticles = currentParticles;

%         trackTime = toc;
%         FPS(f) = 1/trackTime;
%         disp(['frame#: ' num2str(f) ',  FPS = ' num2str(FPS(f))]);
    end
[objDict, objDictNorm] = normalizeMat(objDict);
bestParticle = currentParticles(:,bestParticleInd);
end