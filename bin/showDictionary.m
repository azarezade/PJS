function [] = showDictionary(objDict,objDictNorm,f,settings)

objDict = objDict .* repmat(objDictNorm, size(objDict,1), 1);
image = reshape(objDict, settings.objParam.size(1),[]);
image = imresize(image,3);
figure(2);
set(figure(2),'Name',['Dictionary Window - frame # ' num2str(f)],'NumberTitle','off');
imshow(image);
drawnow

end