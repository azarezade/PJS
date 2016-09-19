 % MFoccusL2L1.CPP Group sparse coding
 % 
 % This code solves the below optimization problem for each group of signals
 %
 % min_P 0.5*||X - D*P||^2 + gamma*SUM_1:K ||P(k,:)||_2
 %
 % Syntax:      P = MFoccusL2L1(X, D, G, P0, param)
 %              P = MFoccusL2L1(X, D, G, P0)
 %
 % Inputs:
 %          D:          Dictionary
 %          X:          Signals
 %          P0:         Initial Guess
 %          G:          Array of zero-based indexes into start of each group
 %          param:      Struct array of following variables
 %              gamma:      Regularization parameter (default = 1)
 %              tol:        Tolerance (default = 10^-2)
 %              max_iter:   Maximum number of iterations (default 100)
 %              nthreads:   Number of threads to use (default No. of CPUs)
 %
 % Output:
 %          P:          Group Sparse Coefficients
 %
 % Written by Ali Soltani-Farani <a_soltani@ce.sharif.edu>
 % Copyright 2012 by Ali Soltani-Farani                                 

 
 % MATLAB Code for a single group
 % function P = MFOCCUS_L2L1_N(X, D, G, P, gamma, tol, max_iter )
 %
 %   gI = gamma*eye(size(D,1),size(D,1));
 %   
 %   delta = 1;
 %   iter = 1;
 %   while delta > tol && iter <= max_iter
 %       L = sum(P.^2,2).^0.25;
 %       Dhat = bsxfun(@times,D,L');
 %       newP = bsxfun( @times, Dhat'*((gI+Dhat*Dhat')\X), L );
 %       delta = sum(sum((P-newP).^2))/sum(sum(P.^2));
 %       P = newP;
 %       iter = iter + 1;
 %   end
 %   
 % end                         