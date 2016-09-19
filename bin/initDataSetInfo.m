function [ frameNum, imgNames, imgSize ] = initDataSetInfo(dataSet)
% Get dir from image folder "dataSetName"  that is in current directory (pwd)
% dataSet must be in the same dicrectory of the main tracker m-file

    dataSetPath = dataSet.path;
    dataSetName = dataSet.name;
    if isequal(dataSet.type,'CVLab')
        imgsAddress = fullfile(dataSetPath,dataSetName,'img');
    elseif isequal(dataSet.type,'Standard')
        imgsAddress = fullfile(dataSetPath,dataSetName,'imgs');
    end
    dataSetDir =dir(imgsAddress);

    % Read images name
    fileNames = dataSetDir(~[dataSetDir.isdir]);    % keep file names only
    imgNames = {fileNames.name}.';                  % all file names in cell array
    imgNames = setdiff(imgNames,'Thumbs.db');
    totalFrameNum = numel(imgNames);                    % number of images

    % Find image size
    if isequal(dataSet.type,'CVLab')
        imgSize = size(imread(fullfile(dataSetPath,dataSetName,'img',imgNames{1})));
    elseif isequal(dataSet.type,'Standard')
        imgSize = size(imread(fullfile(dataSetPath,dataSetName,'imgs',imgNames{1})));
    end
    imgSize = imgSize(1:2);

    % Find number of frames
    if isfield(dataSet,'lastFrame')
        frameNum = min(totalFrameNum, dataSet.lastFrame);
    else
        frameNum = totalFrameNum;
    end
end

