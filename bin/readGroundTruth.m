function [ gtCenters, gtCorners, gtInterv] = readGroundTruth(dataSet,frameNum)

dataSetName=dataSet.name;
dataSetPath=dataSet.path;

if isequal(dataSet.type,'CVLab')
    % read gt file
    gtInterv = 1;
    gtFileAddress = fullfile(dataSetPath, dataSetName, 'groundtruth_rect.txt');
elseif isequal(dataSet.type,'Standard')
    % read gt file
    gtInrevAddress = fullfile(dataSetPath, dataSetName, [dataSetName '_gtInterv.txt']);
    fid1 = fopen(gtInrevAddress);
    gtInterv = str2num(fgetl(fid1));
    gtFileAddress = fullfile(dataSetPath, dataSetName, [dataSetName '_gt.txt']);
end

fid2 = fopen(gtFileAddress);
gt = textscan(fid2, '%f %f %f %f', 'delimiter',',');
gtTopLeftCorner_x = gt{1}(1:end);
gtTopLeftCorner_y = gt{2}(1:end);
gtWidth =gt{3}(1:end);
gtHeight =gt{4}(1:end);

% Find dimention of other corners
gtTopRightCorner_y = gtTopLeftCorner_y;
gtTopRightCorner_x = gtTopLeftCorner_x + gtWidth;

gtBottomLeftCorner_y = gtTopLeftCorner_y + gtHeight;
gtBottomLeftCorner_x = gtTopLeftCorner_x ;

gtBottomRightCorner_y = gtTopLeftCorner_y + gtHeight;
gtBottomRightCorner_x = gtTopLeftCorner_x + gtWidth;

% All corneres of all frames in clockwise direction 
gtCorners = [gtTopLeftCorner_x, gtTopLeftCorner_y,...
                gtTopRightCorner_x, gtTopRightCorner_y,...
                gtBottomRightCorner_x gtBottomRightCorner_y,...
                gtBottomLeftCorner_x, gtBottomLeftCorner_y];

% find center of gt squares
gtCenter_x = gtTopLeftCorner_x + (gtWidth/2);
gtCenter_y= gtTopLeftCorner_y + (gtHeight/2);

gtCenters = [gtCenter_x gtCenter_y];

tmp = floor( (frameNum + gtInterv -1) / gtInterv );

gtCenters = gtCenters(1:tmp,:);
gtCorners = gtCorners(1:tmp,:);

gtFirstFrameFileAddress = fullfile(dataSetPath, dataSetName, [dataSetName '_gtFirstFrame.txt']);
if exist( gtFirstFrameFileAddress, 'file')
    fid3 = fopen(gtFirstFrameFileAddress);
    gt1st = textscan(fid3, '%f %f %f %f %f', 'delimiter',',');
    gt1stCenter_x = gt1st{1}(1:end);
    gt1stCenter_y = gt1st{2}(1:end);
    gt1stWidth = gt1st{3}(1:end);
    gt1stHeight = gt1st{4}(1:end);
    gt1stRotate = gt1st{5}(1:end);
    
    gtCorners(1,:) = paramLoc2Corner([gt1stCenter_x, gt1stCenter_y, gt1stWidth, gt1stHeight, gt1stRotate]);
    gtCenters(1,:) = mean(reshape(gtCorners(1,:),2,[]),2)';
end
end