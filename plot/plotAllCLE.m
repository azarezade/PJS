function [allCLE,allSuccessRates,avgSuccessRates,gtInterv,frameNum] = plotAllCLE(allCLE,allSuccessRates,allDataSetNames,allTrackerNames,allTrackerSettings,allTrackerColors,legendNames,dn,rngRange,param)
for rs=1:length(rngRange)
    % Set figure properties
    if (param.CLE.showAll)
        figure(1);
%         set(figure(1), 'Position', [200 200 1000 500], 'PaperPosition', [0 0 10 5])
    end
    % plot all plots for each rng
    CLE = [];
    for tn=1:length(allTrackerNames)
        % load data set
        if isequal(allTrackerNames{tn},'OAB') || isequal(allTrackerNames{tn},'Frag') || isequal(allTrackerNames{tn},'MIL')
            load([allDataSetNames{dn}, '_', allTrackerNames{tn}, '_', 'rng', num2str(0), allTrackerSettings{tn}, '.mat']);
        else
            load([allDataSetNames{dn}, '_', allTrackerNames{tn}, '_', 'rng', num2str(rngRange(rs)), allTrackerSettings{tn}, '.mat']);
        end
        
        CLE = [CLE; centerDistError'];
        successRate = computeSuccessRate(VOCMeasure,param.threshold);
        allSuccessRates(tn,dn,rs) = successRate;
        
        frameNum = length(resultsCorners);
        % plot figure
        if (param.CLE.showAll)
            plot(1:gtInterv:frameNum,centerDistError,'color',allTrackerColors{tn})
            hold on;
        end
    end

    % Show figure
    if (param.CLE.showAll)
        xlabel('Frame Number')
        ylabel('Center Location Error (CLE)')
%         title([allDataSetNames{dn},' - ','rng=' num2str(rngRange(rs))])
        title(allDataSetNames{dn})
        legend(legendNames,'Location','NorthWest')
        hold off
    end

    % Save figure
    if (param.CLE.saveAll)
        figName = [allDataSetNames{dn}, '_', 'CLE', '_', 'rng', num2str(rngRange(rs))];
        print(gcf,'-djpeg','-r150',fullfile(param.resultsDir,'Figures CLE',figName))
%         print(gcf,'-depsc',fullfile(param.resultsDir,'Figures CLE',figName))
        set(gcf,'PaperPositionMode','auto')
        print(gcf,'-depsc','-r50',fullfile(param.resultsDir,'Figures CLE',figName))
    end

    % add to allCLE
    allCLE = cat(3,allCLE,CLE);

end
avgSuccessRates = mean(allSuccessRates,3);
end