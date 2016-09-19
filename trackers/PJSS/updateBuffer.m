function [bestCandidPatchBuffer] = updateBuffer(bestBuffer,groupSize,objParam,patchIndex)
patchWidth = objParam.patchWidth;
patchHeight = objParam.patchHeight;
patchNum = objParam.patchNum;
    
    bestCandidPatchBuffer = bestBuffer(patchIndex, :);
    bestCandidPatchBuffer = reshape(bestCandidPatchBuffer,patchWidth*patchHeight, patchNum*groupSize);
	bestCandidPatchBuffer= normalizeMat(bestCandidPatchBuffer);
    %***
    bestCandidPatchBufferAugmented = [zeros(patchWidth*patchHeight,patchNum) bestCandidPatchBuffer];
    tmp=zeros(patchWidth*patchHeight,patchNum*(groupSize+1));
    for i=1:patchNum
        tmp(:,(groupSize+1)*(i-1)+1:(groupSize+1)*i) = bestCandidPatchBufferAugmented(:,i:patchNum:end);
    end
    bestCandidPatchBuffer=tmp;
    
end