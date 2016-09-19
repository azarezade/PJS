function newParticles = resampling(currentParticles, likelihood)
    numOfParticles = length(likelihood);
	pdf = likelihood / sum(likelihood);
    newParticlesInd = randsample(1:numOfParticles,numOfParticles,true,pdf);
    newParticles = currentParticles(:,newParticlesInd);
end