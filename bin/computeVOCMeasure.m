function VOC = computeVOCMeasure(resultsCorners, gtCorners, gtInterv)
    ind = 1:gtInterv:size(resultsCorners,1);
    resolution = 100;
    VOC = zeros(size(ind,2),1);
    interSection = zeros(size(ind,2),1);
    union = zeros(size(ind,2),1);
    for i = 1:size(ind,2)
        j=ind(i);
        resCor = reshape(resultsCorners(j,:),2,4)';
        gtCor = reshape(gtCorners(i,:),2,4)';
        
        interSection(i) = areaIntersection(resCor, gtCor,  resolution);
        union(i) = areaIntersection(resCor, resCor, resolution) + areaIntersection(gtCor, gtCor, resolution) - interSection(i);
        VOC(i) = interSection(i) / union(i);
    end
    
end