function [ patchDict,patchDictNorm,objDict,objDictNorm, alpha, beta, wasOccluded] ...
     = updateDict(bestCandidData,bestCandidDataNrom,patchDict,patchDictNorm,objDict,objDictNorm,patchIndex,frameNum,settings, alpha, beta, wasOccluded)
 
    patchWidth = settings.objParam.patchWidth;
    patchHeight = settings.objParam.patchHeight;
    patchNum = settings.objParam.patchNum;

    
    bestCandidPatches = reshape(bestCandidData(patchIndex),patchHeight*patchWidth,[]);

    [bestCandidPatches, bestCandidPatchesNorm] = normalizeMat(bestCandidPatches);
    
    bestObjectPatchCoef = mexLasso(bestCandidPatches, patchDict, settings.SC_Lasso_Param);
    
    isOccluded = ones(settings.objParam.patchNum,1);
    
    alpha = settings.bern_aging_factor * alpha;
    beta = settings.bern_aging_factor * beta;
    
    sigmaInv = 1 / settings.objectModelLikelihoodSigma;%0.01;
    
    for i=1:patchNum
        
        ind_1 = zeros(size(bestObjectPatchCoef(:,i)));
        ind_1(i:settings.objParam.patchNum:end,1) = 1;
        ind_2 = 1 - ind_1;
        
        subPatchCoef_1 = ind_1 .* bestObjectPatchCoef(:,i);
        subPatchCoef_2 = ind_2 .* bestObjectPatchCoef(:,i);
        
        err_1 = norm(bestCandidPatches(:,i) -  (patchDict * subPatchCoef_1))^2;
        err_2 = norm(bestCandidPatches(:,i) -  (patchDict * subPatchCoef_2))^2;
        
%        fprintf('%1.2f\t%1.2f\t%1.2f\n',log(occProb(i)),log(1 - occProb(i)) - err_1, log(occProb(i)) - err_2);
        

        WOIndex = wasOccluded(i) + 1;
        occProb = alpha(WOIndex,i) / (beta(WOIndex,i) + alpha(WOIndex,i));

        
        if sigmaInv * log(1 - occProb) - err_1 > sigmaInv * log(occProb) - err_2
            isOccluded(i) = 0;
            beta(WOIndex,i) = beta(WOIndex,i) + 1;
        else
            alpha(WOIndex,i) = alpha(WOIndex,i) + 1;
        end
        
%        fprintf('%1.8f\t%f\t%f\t%f\t%d\n', log(1-occProb), -err_1, log(occProb), -err_2, isOccluded(i));
        
    end
    
%    '---------------------------------------------'
%    pause;
    
    wasOccluded = isOccluded;
    
%    disp(['frame: ' num2str(frameNum)]);
%    reshape(isOccluded,4,4)
    
    % Update Object Dictionary (Replace a Word Randomly) and Patch
    % Dictionary
    if mod(frameNum,settings.updateInterv) == 0 || frameNum <= settings.dictObjNum
        if frameNum <= settings.dictObjNum
            remIdx = 2;
        else
            random_weight = [0,(2).^(1:settings.dictObjNum-1)];
            random_weight = cumsum(random_weight/sum(random_weight));
            random_num = rand(1,1);
            for remIdx=2:settings.dictObjNum-1
                if random_num >= random_weight(remIdx-1) && random_num<random_weight(remIdx)
                    break;
                end
            end
            if random_num>=random_weight(settings.dictObjNum-1)
                remIdx = settings.dictObjNum;
            end
        end
        
        newPatches = zeros(settings.objParam.patchHeight*settings.objParam.patchWidth, patchNum);
        newPatchesNorm = zeros(1,patchNum);
        for i=1:patchNum
            if ~isOccluded(i) || frameNum <= settings.dictObjNum
                newPatches(:,i) = bestCandidPatches(:,i);
                newPatchesNorm(i) = bestCandidPatchesNorm(i);
            else
                newPatches(:,i) = patchDict(:, (remIdx-1)*patchNum+i);
                newPatchesNorm(i) = patchDictNorm((remIdx-1)*patchNum+i);
            end
        end
        
        newPatches = normalizeMat(newPatches);
        
        patchDict(:,(remIdx-1)*patchNum+1:remIdx*patchNum) = [];
        patchDictNorm((remIdx-1)*patchNum+1:remIdx*patchNum) = [];
        patchDict = [patchDict, newPatches];
        patchDictNorm = [patchDictNorm, newPatchesNorm];

        if (settings.videoParam.show)
        
            objDict(:,remIdx)=[];

            recon = [];
            cnt = 0;
            for i=1:settings.objParam.patchNumW
                tmp = [];
                for j=1:settings.objParam.patchNumH
                    cnt = cnt + 1;
                    tmp = [tmp; reshape(newPatches(:,cnt)*newPatchesNorm(cnt), patchHeight, patchWidth)];
                end
                recon = [recon, tmp];
            end

            objDict(:,settings.dictObjNum) = normalizeMat(recon(:));
            objDictNorm(settings.dictObjNum) = bestCandidDataNrom;
        end
        
        %pause;
    end
end

%++++++++++++++++
% Show eigenbasis
%     figure(3);
%     set(figure(3),'Name',['eigenbasis - frame #' num2str(f)],'NumberTitle','off');
%     imshow(reshape(tmpl.basis,[objSize(1), objSize(2)*size(objDict,2)]),[]);
%     drawnow;
%++++++++++++++++
