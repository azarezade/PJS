function [posres,negres]=getlearnres(clf,posval,negval)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to get likihood of features over samples 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num = size(clf,2);

numpos = size(posval,1);
numneg = size(negval,1);

for i = 1:num
    for j = 1:numpos
        posres(j,i) = computeres(clf{1,i},posval(j,i));
    end
    for k = 1:numneg
        negres(k,i) = computeres(clf{1,i},negval(k,i));
    end
end
