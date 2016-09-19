function [settings]=saveVideo(settings,f,frameNum)

% saveTracker
if (settings.videoParam.save)
    figure(1);
    frame1 = getframe;
    writeVideo(settings.writerMain,frame1);
    if (f==frameNum)
        close(settings.writerMain);
    end
end
end