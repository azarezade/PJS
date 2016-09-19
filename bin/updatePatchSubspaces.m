function [tmpl, patchDict, patchDictNorm, objDict, objDictNorm] = updatePatchSubspaces(isOccluded, tmpl, bestCandidPatches, bestCandidPatchesNorm, frameNum, patchDict, patchDictNorm, objDict, objDictNorm, bestCandidDataNrom, settings)
    
    patchWidth = settings.objParam.patchWidth;
    patchHeight = settings.objParam.patchHeight;
    patchNum = settings.objParam.patchNum;

    for i=1:patchNum
        if isOccluded(i) ~= 1
            tmpl{i}.warpimg = [tmpl{i}.warpimg,bestCandidPatches(:,i)];
            tmpl{i}.warpimgNorm = [tmpl{i}.warpimgNorm, bestCandidPatchesNorm(i)];
        end
    end
    
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
        
        %%%%

        newPatches = zeros(settings.objParam.patchHeight*settings.objParam.patchWidth, patchNum);
        newPatchesNorm = zeros(1,patchNum);
        
        for i=1:patchNum
            
            if size(tmpl{i}.warpimg,2) > 0
                [tmpl{i}.basis, tmpl{i}.eigval, tmpl{i}.mean, tmpl{i}.numsample] = sklm(tmpl{i}.warpimg, tmpl{i}.basis, tmpl{i}.eigval, tmpl{i}.mean, tmpl{i}.numsample, settings.omParam.ff);
                if  (size(tmpl{i}.basis,2) > settings.omParam.maxbasis)
                    tmpl{i}.basis  = tmpl{i}.basis(:,1:settings.omParam.maxbasis);   
                    tmpl{i}.eigval = tmpl{i}.eigval(1:settings.omParam.maxbasis);    
                end
            end
            
            if isOccluded(i) ~= 1
                updatingPatch = bestCandidPatches(:,i);
                updatingPatchNorm = bestCandidPatchesNorm(i);
            else
                if size(tmpl{i}.warpimg,2) > 0
                    updatingPatch = mean(tmpl{i}.warpimg,2);
                    updatingPatchNorm = mean(tmpl{i}.warpimgNorm);
                else
                    updatingPatch = patchDict(:, (remIdx-1)*patchNum+i);
                    updatingPatchNorm = patchDictNorm((remIdx-1)*patchNum+i);
                end
            end
            tmpl{i}.warpimg = [];
            tmpl{i}.warpimgNorm = [];
            
            reconCoef = mexLasso((updatingPatch-tmpl{i}.mean), [tmpl{i}.basis, eye(size(tmpl{i}.basis,1)) ], settings.eigParam.SC_Lasso_Param);
            recon = tmpl{i}.basis*reconCoef(1:size(tmpl{i}.basis,2))+tmpl{i}.mean;
            newPatches(:,i) = normalizeMat(recon);
            newPatchesNorm(i) = updatingPatchNorm;
        end
        
        %%%%
        
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

