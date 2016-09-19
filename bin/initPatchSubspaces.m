function [ tmpl ] = initPatchSubspaces(patchDict,patchDictNorm,settings)
    for i=1:settings.objParam.patchNum
        
        objs = patchDict(:,i:settings.objParam.patchNum:end);
        
        tmpl{i}.mean = mean(objs,2);
        tmpl{i}.basis = [];                                        
        tmpl{i}.eigval = [];                                      
        tmpl{i}.numsample = 0;                                     
        tmpl{i}.warpimg = [];
        [tmpl{i}.basis, tmpl{i}.eigval, tmpl{i}.mean, tmpl{i}.numsample] = sklm(objs, tmpl{i}.basis, tmpl{i}.eigval, tmpl{i}.mean, tmpl{i}.numsample);
        tmpl{i}.warpimgNorm = patchDictNorm(i:settings.objParam.patchNum:end);
    end
end

