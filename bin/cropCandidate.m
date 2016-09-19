function [ candidates, candidates_norm ] = cropCandidate(currentFrame,currentParticles,objParam)

objSize = objParam.size;
candidates = warpimg(currentFrame, paramGeom2Aff(currentParticles), objSize); 
candidates = candidates.*(candidates>0); 
[candidates,candidates_norm] = normalizeMat(reshape(candidates,objSize(1)*objSize(2), size(currentParticles,2)));

end