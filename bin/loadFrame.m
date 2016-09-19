function frame = loadFrame( f, imgNames, settings)
    dataSetPath=settings.dataSet.path;
    dataSetName=settings.dataSet.name;
    
    if isequal(settings.dataSet.type,'CVLab')
        currentFrame = imread(fullfile(dataSetPath,dataSetName,'img',imgNames{f}));
    elseif isequal(settings.dataSet.type,'Standard')
            currentFrame = imread(fullfile(dataSetPath,dataSetName,'imgs',imgNames{f}));
    end
    
    if size(currentFrame,3)==3
        grayFrame = rgb2gray(currentFrame);
    else
        grayFrame = currentFrame;
    end
    
	frame = double(grayFrame)/255;
end