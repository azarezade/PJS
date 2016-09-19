wimgs = [wimgs, param.wimg(:)];
if (size(wimgs,2) >= settings.omParam.batchsize)
    if (isfield(param,'coef'))
      ncoef = size(param.coef,2);
      recon = repmat(tmpl.mean(:),[1,ncoef]) + tmpl.basis * param.coef;
      [tmpl.basis, tmpl.eigval, tmpl.mean, tmpl.numsample] = ...
        sklm(wimgs, tmpl.basis, tmpl.eigval, tmpl.mean, tmpl.numsample, settings.omParam.ff);
      param.coef = tmpl.basis'*(recon - repmat(tmpl.mean(:),[1,ncoef]));
    else
      [tmpl.basis, tmpl.eigval, tmpl.mean, tmpl.numsample] = ...
        sklm(wimgs, tmpl.basis, tmpl.eigval, tmpl.mean, tmpl.numsample, settings.omParam.ff);
    end
    wimgs = [];

    if (size(tmpl.basis,2) > settings.omParam.maxbasis)
      %tmpl.reseig = omParam.ff^2 * tmpl.reseig + sum(tmpl.eigval(tmpl.maxbasis+1:end).^2);
      tmpl.reseig = settings.omParam.ff * tmpl.reseig + sum(tmpl.eigval(settings.omParam.maxbasis+1:end));
      tmpl.basis  = tmpl.basis(:,1:settings.omParam.maxbasis);
      tmpl.eigval = tmpl.eigval(1:settings.omParam.maxbasis);
      if (isfield(param,'coef'))
        param.coef = param.coef(1:settings.omParam.maxbasis,:);
      end
    end
end