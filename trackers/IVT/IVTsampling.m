function [nextParticles] = IVTsampling(prevParticles, numSample, affsig, likelihood)
    cumconf = cumsum(likelihood);
    idx = floor(sum(repmat(rand(1,numSample),[numSample,1]) > repmat(cumconf,[1,numSample])))+1;
    nextParticles = prevParticles(:,idx);
    nextParticles = nextParticles + randn(6,numSample).*repmat(affsig(:),[1,numSample]);
end
