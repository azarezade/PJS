% Written by Ahmad Khajenezhad <khajenezhad@ce.sharif.edu>, Ali Zarezade <zarezade@ce.sharif.edu>
 % Copyright 2012 by Ali Soltani-Farani, Ali Zarezade, Ahmad Khajenezhad
function likelihood = Pooling(patchCoef,objParam, patchDict, candidates, patchIndex, poolInd)

patchNum = objParam.patchNum;
coef = patchCoef .* poolInd;

candidatesPatches = reshape(candidates(patchIndex,:),objParam.patchHeight*objParam.patchWidth,[]);
candidatesPatches = normalizeMat(candidatesPatches);


dists = -sum((candidatesPatches - patchDict * coef).^2);

similarity = sum(reshape(dists, patchNum, []));

likelihood = similarity;
% likelihood = exp(similarity);

end
