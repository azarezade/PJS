function [patchIndex,objParam]= makePatchIndex(objParam)

objSize = objParam.size;
patchHeight = objParam.patchHeight;
patchWidth = objParam.patchWidth;
patchWidthOverlap = objParam.patchWidthOverlap;
patchHeightOverlap = objParam.patchHeightOverlap;

% Computing patchNumH and patchNumW
    patchNumH = floor( (objSize(1) - patchHeightOverlap)/(patchHeight - patchHeightOverlap) );
    if patchNumH * patchHeight - (patchNumH - 1) * patchHeightOverlap ~= objSize(1)
        error('There is some kind of mismatch between patchHeight, patchHeightOverlap and ObjectHeight.');
    end
    patchNumW = floor( (objSize(2) - patchWidthOverlap)/(patchWidth - patchWidthOverlap) );
    if patchNumW * patchWidth - (patchNumW - 1) * patchWidthOverlap ~= objSize(2)
        error('There is some kind of mismatch between patchWidth, patchWidthOverlap and ObjectWidth.');
    end
    patchNum = patchNumH * patchNumW;

    % Make Patch Indices
    patchInd = zeros(patchHeight * patchWidth, patchNum);
    ind = reshape(1:objSize(1)*objSize(2), objSize(1), objSize(2));

    counter = 0;
    for j=1:patchWidth-patchWidthOverlap:objSize(2)-patchWidth+1
        for i=1:patchHeight-patchHeightOverlap:objSize(1)-patchHeight+1
            counter = counter + 1;
            patchInd(:,counter) = reshape(ind(i:i+patchHeight-1,j:j+patchWidth-1), patchHeight*patchWidth, 1);
        end
    end
    clear ind;
    

    patchIndex = reshape(patchInd,[],1);
    objParam.patchNum = patchNum;
    objParam.patchNumH = patchNumH;
    objParam.patchNumW = patchNumW;
end