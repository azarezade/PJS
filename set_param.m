function settings = set_param()

% Initialize Parameters
    settings.dataSet.name = 'Trellis';
    settings.dataSet.type = 'Standard';
    settings.dataSet.path = fullfile(pwd, '..','Standard DataSet');
%     settings.dataSet.lastFrame= 31;
    settings.verifyMode = 0;

% L1 minimization parameters (Lasso in SPAMS library)
    settings.SC_Lasso_Param.mode = 2;
    settings.SC_Lasso_Param.lambda = 0.01;
    settings.SC_Lasso_Param.pos = 'true';
    
% Weighted L1 parameters
    settings.SC_WL1_Param.lambda = 0.0001;
    settings.SC_WL1_Param.alpha = 0.01;
    settings.SC_WL1_Param.betta = 0.01;
    settings.SC_WL1_Param.tol = 5e-2;
    settings.SC_WL1_Param.max_iter = 100;
    settings.SC_WL1_Param.L = 64;
    settings.SC_WL1_Param.numThreads = 1;
    
% SOMP minimization parameters 
    settings.SC_SOMP_Param.L = 4;

% MTT SOMP minimization parameters 
    settings.MTT_SC_SOMP_Param.L = 10;
    
% mFoccus parametes
    settings.SC_mFoccus_Param.gamma = 0.001;
    settings.SC_mFoccus_Param.tol = 1e-3;
    settings.SC_mFoccus_Param.max_iter = 1000;
    settings.SC_mFoccus_Param.nthreads = 4;
    
% OMP minimization parameters 
    settings.OMP_Param.L = 5;
    
% Sigma of Likelihood function of Object Model (used in updateDict to
% balance the ratio of likelihood to the prior)
settings.objectModelLikelihoodSigma = 100;
    
% Patch Kernel
    settings.patchKernel{8} = [
        0.5   0.7   0.7   0.5;
        0.7    1     1    0.7;
        0.7    1     1    0.7;
        0.5   0.7   0.7   0.5];
    
    settings.patchKernel{16} = [
        0.5   0.7   0.5;
        0.7    1    0.7;
        0.5   0.7   0.5];
    
% Size of buffer
    settings.groupSize = 4;    
    
% Eigen Basis L1 minimization parameters
    settings.eigParam.SC_Lasso_Param.mode = 2;
    settings.eigParam.SC_Lasso_Param.lambda = 0.01;
    settings.eigParam.SC_Lasso_Param.pos = 'ture';
    settings.eigParam.updateInterv = 5;

% Update Interval and Bernoulli Aging Factor
    settings.updateInterv = 5;
    settings.bern_aging_factor = 0.1;
        
% Number of dictionary atoms 
    settings.dictObjNum = 10;
    settings.dictBGNum = 1;
%     settings.BGDictUpdateInterv = 5;
    
% Object settings
    settings.objParam.size = [32, 32];
    
% Patch settings
    settings.objParam.patchWidth = 8;%16;
    settings.objParam.patchHeight = 8;%16
    settings.objParam.patchWidthOverlap = 0;%8
    settings.objParam.patchHeightOverlap = 0;%8
    settings.objParam.patchNum = [];  % it will be initialize in "makePatchIndex" function
    settings.objParam.backgroundPathNum = 1;

% Particle filter parameters
    settings.pfParam = struct('numsample', 600, 'affsig', [6,6,.02,.002,.002,0]);

% Initialize Observation Model Parameters
    settings.omParam = struct( 'condenssig',0.2, 'ff',0.95, 'batchsize', 5, 'maxbasis', 10, 'errfunc', 'L2' );
    
% Video settings
    settings.videoParam.save = 0;
    settings.videoParam.show = 1;

% Success Rate Threshold
    settings.VOCThreshold = 0.6;
    
% Make results directory    
    if ~exist('results','dir')
        mkdir('results')
    end
end