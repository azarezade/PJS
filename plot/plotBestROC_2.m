function [ ] = plotBestROC_2(allVOCmean,allDataSetNames,allTrackerNames,allTrackerSettings,allTrackerColors,dn,param)
%     allVOCmean = squeeze(allVOCmean);
%     [~, indx] = max(allVOCmean,[],2);
    figure(4);
    for tn=1:length(allTrackerNames)
        % load data set
        
        ROCData = [];
        for thr=0:0.1:1
            
            if isequal(allTrackerNames{tn},'OAB') || isequal(allTrackerNames{tn},'Frag') || isequal(allTrackerNames{tn},'MIL')
                load([allDataSetNames{dn}, '_', allTrackerNames{tn}, '_', 'rng', num2str(0), allTrackerSettings{tn}, '.mat']);
                bestSuccessRate = computeSuccessRate(VOCMeasure,thr);
            else
                bestSuccessRate = -Inf;
                for rng=0:9
                    load([allDataSetNames{dn}, '_', allTrackerNames{tn}, '_', 'rng', num2str(rng), allTrackerSettings{tn}, '.mat']);
                    bestSuccessRate = max(bestSuccessRate,computeSuccessRate(VOCMeasure,thr));
                end
            end
            
            ROCData = [ROCData, bestSuccessRate];
        end
                
        
        plot(0:0.1:1,ROCData,'Color',allTrackerColors{tn})
        hold on
    end
    title(allDataSetNames{dn})
    xlabel('Threshold')
    ylabel('Success rate (Avg)')
    legend(allTrackerNames,'location','SouthWest')
    hold off
    
    figName = [allDataSetNames{dn}, '_ROC'];
%     print(gcf,'-djpeg','-r150',fullfile(param.resultsDir,'Figures VOC',figName))
%     print(gcf,'-depsc','-r75',fullfile(param.resultsDir,'Figures VOC',figName))
    set(gcf,'PaperPositionMode','auto')
    print(gcf,'-depsc','-r50',fullfile(param.resultsDir,'Figures VOC',figName))

end
