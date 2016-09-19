function [objDict,objDictNorm,resultsCorners,bestParticle ] = ...
            initTracker(objDict,objDictNorm,resultsCorners,initGeoParam,shape,imgNames,frameNum,patchIndex,settings)
objSize = settings.objParam.size;
% simple tracking and collecting tracking results as templates
%% simple tracking

firstFrame = loadFrame( 1, imgNames, settings);
move = [0 0 0 0 0 0;
    2 0 0 0 0 0;
    1 0 0 0 0 0;
    0 1 0 0 0 0;
    -1 0 0 0 0 0;
    0 -1 0 0 0 0;
    1 1 0 0 0 0;
    -1 -1 0 0 0 0;
    1 -1 0 0 0 0;
    -1 1 0 0 0 0];

for i=1:settings.dictObjNum
    particle = initGeoParam' + move(i,:)';
    candidate = warpimg(firstFrame, paramGeom2Aff(particle), objSize);
    candidate = candidate.*(candidate>0); 
    candi_data = candidate(:);
    objDict = [objDict, candi_data];
end

bestParticle = initGeoParam';
bestCorners = paramAff2Corner(paramGeom2Aff(bestParticle), objSize);
resultsCorners = bestCorners;

[objDict, objDictNorm] = normalizeMat(objDict);

% showTracker(firstFrame,bestCorners,1,shape,objDict,objDictNorm,settings);
showTracker(firstFrame,bestCorners,1,objDict,objDictNorm,settings)

end