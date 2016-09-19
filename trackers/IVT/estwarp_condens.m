function param = estwarp_condens(frm, tmpl, param, pfParam, omParam)
% function param = estwarp_condens(frm, tmpl, param, pfParam, omParam)
%

%% Copyright (C) Jongwoo Lim and David Ross.
%% All rights reserved.

n = pfParam.numsample;
sz = size(tmpl.mean);
N = sz(1)*sz(2);


wimgs = warpimg(frm, paramGeom2Aff(param.param), sz);
diff = repmat(tmpl.mean(:),[1,n]) - reshape(wimgs,[N,n]);
coefdiff = 0;
if (size(tmpl.basis,2) > 0)
  coef = tmpl.basis'*diff;
  diff = diff - tmpl.basis*coef;
  if (isfield(param,'coef'))
    coefdiff = (abs(coef)-abs(param.coef))*tmpl.reseig./repmat(tmpl.eigval,[1,n]);
  else
    coefdiff = coef .* tmpl.reseig ./ repmat(tmpl.eigval,[1,n]);
  end
  param.coef = coef;
end
if (~isfield(omParam,'errfunc'))  omParam.errfunc = [];  end
switch (omParam.errfunc)
  case 'robust';
    param.conf = exp(-sum(diff.^2./(diff.^2+omParam.rsig.^2))./omParam.condenssig)';
  case 'ppca';
    param.conf = exp(-(sum(diff.^2) + sum(coefdiff.^2))./omParam.condenssig)';
  otherwise;
    param.conf = exp(-sum(diff.^2)./omParam.condenssig)';
end
param.conf = param.conf ./ sum(param.conf);
[maxprob,maxidx] = max(param.conf);
param.est = param.param(:,maxidx);
param.wimg = wimgs(:,:,maxidx);
param.err = reshape(diff(:,maxidx), sz);
param.recon = param.wimg + param.err;
