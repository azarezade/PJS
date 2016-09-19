function [allCLEAvg,allCLEvar] = plotAvgCLE(allCLE,allDataSetNames,allTrackerNames,allTrackerColors,legendNames,dn,gtInterv,frameNum,param)
    allCLEAvg = mean(allCLE,3);
    allCLEvar = std(allCLE,0,3);

    % write legends
    allLegendNames = cell(1,length(allTrackerNames));
    for tn=1:length(allTrackerNames)
        allLegendNames{tn} = [allTrackerNames{tn} ' avg = ' num2str(mean(allCLEAvg(tn,:)))  ' std = ' num2str(std(allCLEAvg(tn,:)))];
    end
    % plot the averaged CLE    
    if (param.CLE.showAvg)
        % set figure properties
        figure(2);
%         set(figure(2), 'Position', [200 200 1000 500], 'PaperPosition', [0 0 10 5])
        
        % plot avg cle
        for tn=1:length(allTrackerNames)
            % plot current figure
            plot(1:gtInterv:frameNum,allCLEAvg(tn,:),'color',allTrackerColors{tn})
            hold on;
        end
        xlabel('frame number')
        ylabel('center location error (CLE)')
%         title([allDataSetNames{dn},' - ','average rng 0 to 9'])
        title(allDataSetNames{dn})
%         legend(allLegendNames,'location','northwest')
        legend(legendNames,'location','northwest')
        hold off
        % save figure
        if (param.CLE.saveAvg)
            figName = [allDataSetNames{dn}, '_', 'CLE_', 'rngavg0to9'];
            print(gcf,'-djpeg','-r150',fullfile(param.resultsDir,'Figures CLE',figName))
%             print(gcf,'-depsc',fullfile(param.resultsDir,'Figures CLE',figName))
            set(gcf,'PaperPositionMode','auto')
            print(gcf,'-depsc','-r50',fullfile(param.resultsDir,'Figures CLE',figName))
        end
    end
end