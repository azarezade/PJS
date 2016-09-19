    currentFrame = loadFrame( f, imgNames, settings);
    gtInd = (f-1)/gtInterv + 1;

    w = norm(gtCorners(gtInd,3:4) - gtCorners(gtInd,1:2));
    h = norm(gtCorners(gtInd,7:8) - gtCorners(gtInd,1:2));

    diversity = [40,40]; 
    stpSize = 2;
%     stpSize = 1;
    firstW = gtCorners(gtInd,1) - diversity(1) * stpSize;
    firstH = gtCorners(gtInd,2) - diversity(2) * stpSize;
    lastW = gtCorners(gtInd,1) + diversity(1) * stpSize;
    lastH = gtCorners(gtInd,2) + diversity(2) * stpSize;

    %%%%%%

    particles = zeros(6,(2*diversity(1)+1) * (2*diversity(2)+1));
    cnt = 0;
    for i=firstW:stpSize:lastW
        for j=firstH:stpSize:lastH
            corners = [i,j,i+w-1,j,i+w-1,j+h-1,i,j+h-1];
            center = [i+floor(w/2)-1, j+floor(h/2)-1];
            loc = paramCornersAndCenter2Loc(corners,center);
            geom = paramLoc2Geom(loc,settings.objParam.size(1));

            cnt = cnt + 1;
            particles(:,cnt) = geom(:);
        end
    end

    %%%%%%%%%%%%%

%         bestParticle = paramLoc2Geom(paramCornersAndCenter2Loc(gtCorners(gtInd,:), gtCenters(gtInd,:)),settings.objParam.size(2))';
%         newPfParam = struct('numsample', 4000, 'affsig', [4, 4, .01, .01, .001,.000]);
%         particles = sampling(bestParticle,newPfParam);

    %%%%%%%%%%%%5

%         bestParticle = paramLoc2Geom(paramCornersAndCenter2Loc(gtCorners(gtInd,:), gtCenters(gtInd,:)),settings.objParam.size(2))';
%         particles = sampling(bestParticle,settings.pfParam);

    [ candidates, candidates_norm ] = cropCandidate(currentFrame, particles,settings.objParam);
    if isequal(trackerName,'PBGS') || isequal(trackerName,'PBGM')
        initPoint = ones(size(patchDict,2),size(candidates,2) * size(groupBuffer,2));
        newListGroups = int32((1:settings.groupSize+1:(settings.groupSize+1)*settings.objParam.patchNum*size(candidates,2))-1);
    end

    algs = {'APG'}; % {'IVT','APG','L1','LST','MTT_L2L0','PBGS','PBGM'};
    for algI = 1:length(algs)
        algName = algs{algI};

        addpath(fullfile(pwd,'trackers',algName),'-BEGIN');

        if isequal(algName,'IVT')
            param.param = currentParticles;
            % Compute Likelihood
            param = estwarp_condens(currentFrame, tmpl, param, settings.pfParam, settings.omParam);
            likelihood = param.conf;
        end

        if isequal(algName,'PBGM')
            poolInd = repmat(eye(settings.objParam.patchNum), settings.dictObjNum, size(candidates,2));
            patchCoef = Coding(candidates, groupBuffer, patchDict, newListGroups, patchIndex, initPoint, settings);
            likelihood = exp(Pooling(patchCoef,settings.objParam, patchDict, candidates, patchIndex, poolInd));
        end

        if isequal(algName,'PBGS')
            poolInd = repmat(eye(settings.objParam.patchNum), settings.dictObjNum, size(candidates,2));
            patchCoef = Coding(candidates, groupBuffer, patchDict, newListGroups, patchIndex, [], settings);
            likelihood = exp(Pooling(patchCoef,settings.objParam, patchDict, candidates, patchIndex, poolInd));
        end

        if isequal(algName,'LST')
            patchCoef = Coding(candidates,patchDict, settings.SC_Lasso_Param,settings.objParam,patchIndex);
            likelihood = Pooling(patchCoef,settings.objParam,settings.dictObjNum);
        end

        if isequal(algName,'MTT_L2L0')
            coef = Coding(candidates,T, settings);
            likelihood = Pooling(candidates,coef,T);
        end

        if isequal(algName,'L1')
            coef = Coding(candidates,objDict, settings);
            likelihood = Pooling(candidates,coef,objDict);
        end

        if isequal(algName,'APG')
            APG_step;
            likelihood = p';
        end
        
%         for p_i = 1:size(particles,2)
%             particle = particles(:,p_i);
%             corners = paramAff2Corner(paramGeom2Aff(particle), settings.objParam.size);
%             for c_i = 1:2:7
%                 if corners(c_i) < 1 || corners(c_i) > size(currentFrame,2) ...
%                         || corners(c_i+1) < 0 || corners(c_i+1) > size(currentFrame,1)
%                     likelihood(p_i) = -inf;
%                     break;
%                 end
%             end
%         end


        fig = figure;
        set(fig,'Name',algName,'NumberTitle','off');
        colormap('jet')
        colors=colormap;
        dinfo = [particles(1,:); particles(2,:); likelihood];
        maxl=max(dinfo(3,:));
        minl=min(dinfo(3,:));
        imshow(currentFrame);
        hold on            
        for i=1:size(dinfo,2)
            cl = colors(floor((length(colors)-1)*(dinfo(3,i)-minl)/(maxl-minl)+1),:);
            plot(dinfo(1,i),dinfo(2,i),'s','MarkerFacecolor',cl, 'MarkerEdgecolor',cl);
        end
    end
    