function probsum=classify(samples,selectors,clf,image,ftrs,integral_image)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to estimate candidates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

samnum = size(samples,1);
res = zeros(samnum,1);

selnum = size(selectors,1);

val = zeros(samnum,selnum);

for i = 1:selnum
    for j = 1:samnum
        val(j,i) = computeftr(ftrs{1,selectors(i,1)},image,samples(j,:),integral_image);
        prob(j,i) = computeres(clf{1,selectors(i,1)},val(j,i));
    end
end

probsum = sum(prob,2);