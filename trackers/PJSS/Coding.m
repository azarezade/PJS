function [ patchCoef, initPoint ] = Coding(candidates,groupBuffer,patchDict,listGroups, patchIndex, initPoint, settings)

objParam = settings.objParam;
groupSize = settings.groupSize;

% ------
    candidatesPatches = reshape(candidates(patchIndex,:),objParam.patchHeight*objParam.patchWidth,[]);
    candidatesPatches = normalizeMat(candidatesPatches);
    
% ------

    candidPatchGroup = repmat(groupBuffer,[1,size(candidates,2)]);
    candidPatchGroup(:,1:groupSize+1:end)=candidatesPatches;
    
% sparse coding
    patchCoef = mexSOMP(candidPatchGroup, patchDict, listGroups, settings.SC_SOMP_Param);
    patchCoef = patchCoef(:,1:groupSize+1:end);
%      size(find(isnan(patchCoef)))
end