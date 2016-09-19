function [ ] = plotBestROC(allVOCmean,allDataSetNames,allTrackerNames,allTrackerSettings,allTrackerColors,dn,param)
    allVOCmean = squeeze(allVOCmean);
    [~, indx] = max(allVOCmean,[],2);
    figure(4);
    for tn=1:length(allTrackerNames)
        % load data set
        if isequal(allTrackerNames{tn},'OAB') || isequal(allTrackerNames{tn},'Frag') || isequal(allTrackerNames{tn},'MIL')
            load([allDataSetNames{dn}, '_', allTrackerNames{tn}, '_', 'rng', num2str(0), allTrackerSettings{tn}, '.mat']);
        else
            load([allDataSetNames{dn}, '_', allTrackerNames{tn}, '_', 'rng', num2str(indx(tn)-1), allTrackerSettings{tn}, '.mat']);
        end
        
        ROCData = [];
        for thr=0:0.1:1
            ROCData = [ROCData, computeSuccessRate(VOCMeasure,thr)];
        end
        plot(0:0.1:1,ROCData,'Color',allTrackerColors{tn})
        hold on
    end
    title(allDataSetNames{dn})
    xlabel('Threshold')
    ylabel('Success rate (Best)')
    legend(allTrackerNames,'location','SouthWest')
    hold off
    
    figName = [allDataSetNames{dn}, '_ROC'];
%     print(gcf,'-djpeg','-r150',fullfile(param.resultsDir,'Figures VOC',figName))
%     print(gcf,'-depsc','-r75',fullfile(param.resultsDir,'Figures VOC',figName))
    set(gcf,'PaperPositionMode','auto')
    print(gcf,'-depsc','-r50',fullfile(param.resultsDir,'Figures VOC',figName))

end
