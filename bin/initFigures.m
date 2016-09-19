function [shape] = initFigures(firstFrame,firstCorner,videoParam)
%  Initialize Graphics            
shape = vision.ShapeInserter('Shape','Polygons',...
    'Fill',false,...
    'BorderColorSource', 'Property',...
    'BorderColor', 'Custom',...
    'CustomBorderColor',[255 255 0],...
    'Antialiasing',false);

if (videoParam.show == true)
    % draw initial track window
        % showFirstFrameWithGt;
            image = zeros(size(firstFrame,1), size(firstFrame,2), 3);
            image(:,:,1) = firstFrame;
            image(:,:,2) = firstFrame;
            image(:,:,3) = firstFrame;
            image = step(shape, image, floor(firstCorner));
            figure(1);
            set(figure(1),'Name',['Main Window - frame # 1' ],'NumberTitle','off');
            imshow(image);
            drawnow
            % pause(gtInterv/24)


    % Show Dictionary
        % showInitialDictionary;
%     figure(2);
%     set(figure(2), 'position', [500 50 800 300], 'paperposition', [0 0 8 3])
%     set(figure(2),'Name','Doctionary - frame #1','NumberTitle','off');
%     %imshow(reshape(objDict(:,1:dictObjNum) .* repmat(objDictNorm(1:dictObjNum),size(objDict,1),1), objSize(1), objSize(2) * dictObjNum));
%     drawnow;

%     disp('Press any key to coninue.');
%     pause;
end

end