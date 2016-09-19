function selectors = getselectors(posval,negval,selnum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get selectors to construct strong clasifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[possamnum,ftrnum] = size(posval);
[negsamnum] = size(negval,1);

%initialize H
Hpos = zeros(possamnum,1);
Hneg = zeros(negsamnum,1);

selectors = zeros(selnum,1);

for m= 1:selnum
    poslikl = ones(ftrnum,1);
    neglikl = zeros(ftrnum,1);
    likelihood = zeros(ftrnum,1);
    for i = 1:ftrnum
        lll = 1.0;
        for j = 1:possamnum
            lll = lll*(1-sigmoid(Hpos(j,1)+posval(j,i)));
        end
        poslikl(i,1) = -log(1-lll+0.00001);
    
        lll = 0;
        for k = 1:negsamnum
            lll = lll-log(0.00001+1-sigmoid(Hneg(k,1)+negval(k,i)));
        end
        neglikl(i,1) = lll;
    end

    likelihood = poslikl/possamnum + neglikl/negsamnum;

    [newlikl,id] = sort(likelihood,1,'ascend');


     
      
     for i = 1:size(id,1)
         if sum(sum(selectors==id(i,1),1),2)==0
            selectors(m,1) = id(i,1);
            break;
         end
     end


    %update H = H + h_m
    Hpos = Hpos + posval(:,selectors(m,1));
    Hneg = Hneg + negval(:,selectors(m,1));

end
    

