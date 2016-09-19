function [ ] = plotAvgROC(allVOCmean,allDataSetNames,allTrackerNames,allTrackerSettings,allTrackerColors,dn,rngRange,param)
%     allVOCmean = squeeze(allVOCmean);
%     [~, indx] = max(allVOCmean,[],2);
    figure(5);
    for tn=1:length(allTrackerNames)
        % load data set
        
        ROCData = [];
        for thr=0:0.1:1
            
            if isequal(allTrackerNames{tn},'OAB') || isequal(allTrackerNames{tn},'Frag') || isequal(allTrackerNames{tn},'MIL')
                load([allDataSetNames{dn}, '_', allTrackerNames{tn}, '_', 'rng', num2str(0), allTrackerSettings{tn}, '.mat']);
                avgSuccessRate = computeSuccessRate(VOCMeasure,thr);
            else
                sumSuccessRate = 0;
                for rng=rngRange
                    load([allDataSetNames{dn}, '_', allTrackerNames{tn}, '_', 'rng', num2str(rng), allTrackerSettings{tn}, '.mat']);
                    sumSuccessRate = sumSuccessRate + computeSuccessRate(VOCMeasure,thr);
                end
                avgSuccessRate = sumSuccessRate / 10;
            end
            
            ROCData = [ROCData, avgSuccessRate];
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
