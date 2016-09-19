% function [nextParticles] = sampling(prevParticles, pfParam, bestParticleInd)
%     prevParticles = repmat(prevParticles(:,bestParticleInd), 1, pfParam.numsample);
%     randomnum = randn(6,pfParam.numsample);
%     nextParticles = ( prevParticles + (randomnum).*repmat(pfParam.affsig(:),[1,pfParam.numsample]) );
% end

function [nextParticles] = sampling(bestParticle, pfParam)
    prevParticles = repmat(bestParticle, 1, pfParam.numsample);
    randomnum = randn(6,pfParam.numsample);
    nextParticles = ( prevParticles + (randomnum).*repmat(pfParam.affsig(:),[1,pfParam.numsample]) );
end
