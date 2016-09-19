n_sample=size(particles,2);

Y = candidates;

% [Y, Y_inrange] = crop_candidates(im2double(currentFrame), paramGeom2Aff(particles)',settings.objParam.size);
% if(sum(Y_inrange==0) == n_sample)
%     sprintf('Target is out of the frame!\n');
% end

[Y,Y_crop_mean,Y_crop_std] = whitening(Y);	 % zero-mean-unit-variance
[Y, Y_crop_norm] = normalizeTemplates(Y); %norm one

%-L1-LS for each candidate target
eta_max	= -inf;
q   = zeros(n_sample,1); % minimal error bound initialization

% first stage L2-norm bounding    
for j = 1:n_sample
%     if Y_inrange(j)==0 || sum(abs(Y(:,j)))==0
%         continue;
%     end

    % L2 norm bounding
    q(j) = norm(Y(:,j)-Temp1*Y(:,j));
    q(j) = exp(-alpha*q(j)^2);
end
%  sort samples according to descend order of q
[q,indq] = sort(q,'descend');    

% second stage
p	= zeros(n_sample,1); % observation likelihood initialization
n = 1;
tau = 0;
while (n<n_sample)&&(q(n)>=tau)        

    [c] = APGLASSOup(Temp'*Y(:,indq(n)),Dict,para);

    D_s = (Y(:,indq(n)) - [A(:,1:nT) fixT]*[c(1:nT); c(end)]).^2;%reconstruction error
    p(indq(n)) = exp(-alpha*(sum(D_s))); % probability w.r.t samples
    tau = tau + p(indq(n))/(2*n_sample-1);%update the threshold

    if(sum(c(1:nT))<0) %remove the inverse intensity patterns
        continue;
    elseif(p(indq(n))>eta_max)
        id_max	= indq(n);
        c_max	= c;
        eta_max = p(indq(n));
%         Min_Err(f) = sum(D_s);
    end
    n = n+1;
end