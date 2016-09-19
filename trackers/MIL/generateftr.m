function [ftr] = generateftr(ftrpara)

%   ---------------------------------> x
%   |
%   |
%   |
%   |
%   |
%   |
%   |
%   \/
%   y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to generate a harr-like feature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

width = ftrpara.width;
height = ftrpara.height;

numrects = randi([ftrpara.minnumrect,ftrpara.maxnumrect],1,1);
weights = zeros(numrects,1);
rects = zeros(numrects,4);
for i = 1:numrects
    weights(i,1) = rand(1)*2 -1;
    rects(i,1) = randi([0,width-2]); % x of top left
    rects(i,2) = randi([0,height-2]);% y of top left 
    rects(i,3) = randi([1,width-rects(i,1)-1]); %width
    rects(i,4) = randi([1,height-rects(i,2)-1]);%height
    rsums = abs(weights(i,1)*(rects(i,3)+1)*(rects(i,4)+1)*255);
end

ftr.rects = rects;
ftr.weights = weights;
ftr.rsums = rsums;