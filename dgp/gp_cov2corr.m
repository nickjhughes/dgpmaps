
function A = gp_cov2corr(A)
%GP_COV2CORR Convert a covariance matrix into a correlation matrix.
%
% A = gp_cov2corr(A)

stdA = repmat(sqrt(diag(A)),1,size(A,1));
A = A./stdA./(stdA');
