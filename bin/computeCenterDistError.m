function centerDistError = computeCenterDistError(resultsCorners, gtCorners, gtInterv)
    resultsCenters=[resultsCorners(:,1)+resultsCorners(:,5) resultsCorners(:,2)+resultsCorners(:,6)]/2;
    gtCenters=[gtCorners(:,1)+gtCorners(:,5) gtCorners(:,2)+gtCorners(:,6)]/2;
    ind = 1:gtInterv:size(resultsCenters,1);
    difference = (resultsCenters(ind, :) - gtCenters);
    centerDistError = sqrt(sum( difference.^2, 2));
end