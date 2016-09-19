function [allVOC] = plotAllVOC(allVOC,allDataSetNames,allTrackerNames,allTrackerSettings,allTrackerColors,legendNames,dn,rngRange,param)
for rs=1:length(rngRange)
    if (param.VOC.showAll)
        % Set figure properties
        figure(1);
%         set(figure(1), 'Position', [200 200 1000 500], 'PaperPosition', [0 0 10 5])
    end
    % plot all plots for each rng
    VOC = [];
    for tn=1:length(allTrackerNames)
        % load data set
        if isequal(allTrackerNames{tn},'OAB') || isequal(allTrackerNames{tn},'Frag') || isequal(allTrackerNames{tn},'MIL')
            load([allDataSetNames{dn}, '_', allTrackerNames{tn}, '_', 'rng', num2str(0), allTrackerSettings{tn}, '.mat']);
        else
            load([allDataSetNames{dn}, '_', allTrackerNames{tn}, '_', 'rng', num2str(rngRange(rs)), allTrackerSettings{tn}, '.mat']);
        end
 
        VOC = [VOC; VOCMeasure'];

        % plot figure
        if (param.VOC.showAll)
            plot(1:gtInterv:length(resultsCorners),VOCMeasure,'color',allTrackerColors{tn})
            hold on;
        end
    end

    % Show figure
    if (param.VOC.showAll)
        xlabel('Frame Number')
        ylabel('Overlap ratio (VOC)')
%         title([allDataSetNames{dn},' - ','rng=' num2str(rngRange(rs))])
        title(allDataSetNames{dn})
        legend(legendNames,'Location','SouthWest')
        hold off
    end

    % Save figure
    if (param.VOC.saveAll)
        figName = [allDataSetNames{dn}, '_VOC', '_', 'rng', num2str(rngRange(rs))];
        print(gcf,'-djpeg','-r150',fullfile(param.resultsDir,'Figures VOC',figName))
%         print(gcf,'-depsc',fullfile(param.resultsDir,'Figures VOC',figName))
        set(gcf,'PaperPositionMode','auto')
        print(gcf,'-depsc','-r50',fullfile(param.resultsDir,'Figures VOC',figName))
    end

    % add to allVOC
    allVOC = cat(3,allVOC,VOC);
end
end