function [allVOCAvg,allVOCvar] = plotAvgVOC(allVOC,allDataSetNames,allTrackerNames,allTrackerColors,legendNames,dn,gtInterv,frameNum,param)
    allVOCAvg = mean(allVOC,3);
    allVOCvar = std(allVOC,0,3);
    % write legends
    allLegendNames = cell(1,length(allTrackerNames));
    for tn=1:length(allTrackerNames)
%         allLegendNames{tn} = [allTrackerNames{tn} ' avgVOC = ' num2str(mean(allVOCAvg(tn,:)))  ' stdVOC = ' num2str(std(allVOCAvg(tn,:)))];
        allLegendNames{tn} = [allTrackerNames{tn} ' avg = ' sprintf('%1.3f \n', mean(allVOCAvg(tn,:)))  ' std = ' sprintf('%1.3f \n', std(allVOCAvg(tn,:)))];
    end
    % plot the averaged VOC    
    if (param.VOC.showAvg)
        % set figure properties
        figure(2);
%         set(figure(2), 'Position', [200 200 1000 500], 'PaperPosition', [0 0 10 5])

        % plot avg VOC
        for tn=1:length(allTrackerNames)
            % plot current figure
            plot(1:gtInterv:frameNum,allVOCAvg(tn,:),'color',allTrackerColors{tn},'LineWidth',2)
            axis([0 frameNum+2 0 1]);
            hold on;
        end
        xlabel('frame number')
        ylabel('Overlap ratio (VOC)')
        title(allDataSetNames{dn})
%         title([allDataSetNames{dn},' - ','average rng 0 to 9'])
%         legend(allLegendNames,'location','SouthWest')
        legend(legendNames,'location','SouthWest')
        hold off
        % save figure
        if (param.VOC.saveAvg)
            figName = [allDataSetNames{dn}, '_', 'VOC_', 'rngavg0to9'];
set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold');
set(gca,'fontsize',13,'fontWeight','bold');                        
            print(gcf,'-djpeg','-r150',fullfile(param.resultsDir,'Figures VOC',figName))
%             print(gcf,'-depsc',fullfile(param.resultsDir,'Figures VOC',figName))
            set(gcf,'PaperPositionMode','auto')
            print(gcf,'-depsc','-r50',fullfile(param.resultsDir,'Figures VOC',figName))
        end
    end
end
