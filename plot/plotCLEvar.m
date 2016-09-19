function [ ] = plotCLEvar(allCLE,allCLEAvg,allCLEvar,allDataSetNames,allTrackerNames,allTrackerColors,dn,param)
    if (param.CLE.showVar)
%         allLegendNames = cell(1,length(allTrackerNames));
%         for tn=1:length(allTrackerNames)
%             allLegendNames{tn} = [allTrackerNames{tn} ' avgcle = ' num2str(mean(allCLEAvg(tn,:)))  ' varcle = ' num2str(std(allCLEAvg(tn,:)))];
%         end

        % set figure properties
        figure(3);
        set(figure(3), 'Position', [200 200 1000 500], 'PaperPosition', [0 0 8 3])
        % plot CLE var
        for tn=1:length(allTrackerNames)
            plot(1:gtInterv:length(resultsCorners),allCLEAvg(tn,:),'color',allTrackerColors{tn})
            hold on
        end 
        xlabel('Frame Number')
        ylabel('Center Location Error (CLE)')
        title([allDataSetNames{dn},' - ','average rng 0 to 9'])
%         legend(legendNames,'Location','NorthWest')
        hold on
        for tn=1:length(allTrackerNames)
            jbfill(1:size(allCLE,2), allCLEAvg(tn,:)+allCLEvar(tn,:), allCLEAvg(tn,:)-allCLEvar(tn,:), allTrackerColors{tn},allTrackerColors{tn},0,0.15); 
            hold on
        end
        hold off
        % save figure
        if (param.CLE.saveVar)
            figName = [allDataSetNames{dn}, '_','CLEvar_', 'rngAvg0to9_var'];
            print(gcf,'-djpeg','-r150',fullfile(param.resultsDir,'Figures CLE',figName))
            print(gcf,'-depsc',fullfile(param.resultsDir,'Figures CLE',figName))
        end
    end
end