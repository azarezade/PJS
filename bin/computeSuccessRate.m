function successRate = computeSuccessRate(VOCMesure,threshold)
    successRate = sum(VOCMesure > threshold)/size(VOCMesure,1);
end
