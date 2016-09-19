function [] = showTracker(img,corners,f,objDict,objDictNorm,settings)
if (settings.videoParam.show)
    imgTracker=repmat(img,[1,1,3]);
    
    % Multiply object dictionry atom by their norms
    objDict = objDict .* repmat(objDictNorm, size(objDict,1), 1);
    imageDict = reshape(objDict, settings.objParam.size(1),[]);
    
    % Glue tracker and dictionary
    scale=size(imgTracker,2)/size(imageDict,2);
    imageDict = imresize(imageDict,scale);
    imageDict = repmat(imageDict,[1,1,3]);
    img = [imgTracker; imageDict];
    
    % Resize final image and corners 
    scale2 = 1;
    img = imresize(img,scale2);
    corners = scale2 * corners;
    cornersX = corners(1:2:end);
    cornersY = corners(2:2:end);
    
    % Draw
    figure(1);
    set(figure(1),'Name',['Main Window - frame # ' num2str(f)],'NumberTitle','off');
    imshow(img,[]);
	hold on;
    plot([cornersX,cornersX(1)],[cornersY cornersY(1)],'y-');
    drawnow
end
end