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

    listGroups(length(listGroups)+1) = size(candidPatchGroup,2);
    initPoint = MFoccusL2L1(candidPatchGroup, patchDict, double(listGroups'), initPoint, settings.SC_mFoccus_Param ); 
    patchCoef = initPoint(:,1:groupSize+1:end);
    
%     size(find(isnan(patchCoef)))

%     imshow(initPoint)
%     x = (sum(abs(patchCoef))./sqrt(sum(patchCoef.^2)));
    
end