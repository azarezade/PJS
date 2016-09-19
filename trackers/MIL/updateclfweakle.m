function newclfle = updateclfweakle(posval,negval,oldclf,params,trainflag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to update weak classifiers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num = size(oldclf,2);%num of weak classifiers

numpos = size(posval,1);
numneg = size(negval,1);

posw = zeros(numpos,1);
negw = zeros(numneg,1);

learnrate = params.lrate;
if trainflag == 0 % the first frame
    for i = 1:num
        newclfle{1,i}.mu0 = mean(negval(:,i));
        newclfle{1,i}.sig0 = var(negval(:,i))+0.000000001;
        newclfle{1,i}.mu1 = mean(posval(:,i));
        newclfle{1,i}.sig1 = var(posval(:,i))+0.000000001;
    end
end
if trainflag == 1 %run-time stage
    for i = 1:num
        newclfle{1,i}.mu0 = mean(negval(:,i))*(1-learnrate) + oldclf{1,i}.mu0*learnrate;
        newclfle{1,i}.sig0 = mean((negval(:,i)-repmat(oldclf{1,i}.mu0,[numneg,1])).^2)*(1-learnrate) + oldclf{1,i}.sig0*learnrate;
        newclfle{1,i}.mu1 = mean(posval(:,i))*(1-learnrate) + oldclf{1,i}.mu1*learnrate;
        newclfle{1,i}.sig1 = mean((posval(:,i)-repmat(oldclf{1,i}.mu1,[numpos,1])).^2)*(1-learnrate) + oldclf{1,i}.sig1*learnrate;
    end
end

for i = 1:num
    newclfle{1,i}.q = (newclfle{1,i}.mu1 - newclfle{1,i}.mu0)/2;
    newclfle{1,i}.s = sign(newclfle{1,i}.mu1 - newclfle{1,i}.mu0);
    newclfle{1,i}.n0 = 1/(newclfle{1,i}.sig0^0.5);
    newclfle{1,i}.n1 = 1/(newclfle{1,i}.sig1^0.5);
    newclfle{1,i}.e0 = -1/(newclfle{1,i}.sig0*2+eps);
    newclfle{1,i}.e1 = -1/(newclfle{1,i}.sig1*2+eps);
end