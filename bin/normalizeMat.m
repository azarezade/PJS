function [A, A_norm] = normalizeMat(A)
    A_norm = sqrt(sum(A.*A));
    A = A./(repmat(A_norm,[size(A,1) 1])+eps);
end