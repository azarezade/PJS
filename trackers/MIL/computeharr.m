function [posval,negval]=computeharr(possamples,negsamples,ftrs,image,integral_image)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to calculate value of harr-like features over negative and positive samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ftrnum = size(ftrs,2);
posnum = size(possamples,1);
negnum = size(negsamples,1);
posval = zeros(posnum,ftrnum);
negval = zeros(negnum,ftrnum);


for i = 1:ftrnum
    for j = 1:posnum
        posval(j,i) = computeftr(ftrs{1,i},image,possamples(j,:),integral_image);
    end
    for k = 1:negnum
        negval(k,i) = computeftr(ftrs{1,i},image,negsamples(k,:),integral_image);
    end
end

 