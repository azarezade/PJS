function clfle = iniclfweakle(params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialization of strong classifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num = params.numfeat;


    weakle.mu0 = 0;
    weakle.mu1 = 0;
    weakle.sig0 = 1;
    weakle.sig1 = 1;
    weakle.q = -431602080;
    weakle.s = -842150451;
    weakle.n0 = -431602080;
    weakle.n1 = -431602080;
    weakle.e0 = -431602080;
    weakle.e1 = -431602080;
    
for i = 1:num
    clfle{1,i} = weakle;
end

%clfle{:}=repmat(weakle,[num,1]);