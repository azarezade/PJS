function [] = playGroundTruth()
settings = setSettings();
% Data Set Information
    [ frameNum, imgNames, ~ ] = initDataSetInfo(settings.dataSet);

% Read Ground Truth
    [~, gtCorners,gtInterv] = readGroundTruth(settings.dataSet, frameNum);
    
% Plot ground truth squre on data set
for i=1:(floor(frameNum/gtInterv))
    img = loadFrame(i,imgNames,settings);
    
    cornersX = gtCorners(i,1:2:end);
    cornersY = gtCorners(i,2:2:end);
    
    figure(1);
    set(figure(1),'Name',['Main Window - frame # ' num2str(i)],'NumberTitle','off');
    imshow(img,[]);
	hold on;
    plot([cornersX,cornersX(1)],[cornersY cornersY(1)],'y-');
    drawnow
    pause(gtInterv/30)

end