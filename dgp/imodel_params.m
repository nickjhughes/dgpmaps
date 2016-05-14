
function params = imodel_params(op, a0)
%IMODEL_PARAMS Estimate hyperparameters using the given vector averaged map.
%
% params = imodel_params(op, a0)
%
% Estimates parameters for the independent (i.e., OP only ) model.
% Estimates parameters by maximising the correlation between the empirical
% and estimated prior kernel functions. params is a struct with fields
% alpha and sigma. a0 is a vector of the initial guesses for each
% hyperparameter, in the same order.

% Initialise
kernel_types;
n = size(op,1);
op_x = real(op);
op_y = imag(op);

% Calculate covariances
xc_op_x = xcorr2(op_x);
xc_op_x = xc_op_x/max(xc_op_x(:));
xc_op_y = xcorr2(op_y);
xc_op_y = xc_op_y/max(xc_op_y(:));

% Fit
if nargin < 3
    % alpha, sigma
    a0 = [1, 2];
end
options = optimset('Display', 'Iter');
a = fminsearch(@(a)(ifit(a, n, xc_op_x, xc_op_y)), a0, options);
params = struct('alpha', a(1), 'sigma', a(2));

end


function err = ifit(a, n, xc_op_x, xc_op_y)

% Initialise
kernel_types;
[xx, yy] = meshgrid(-n:n, -n:n);
x = [yy(:) xx(:)]';
N = 400;

% Parameters
params_G = struct('alpha',a(1), 'sigma',a(2));
kern_G = kernels.dog;

% OP_x = (G*G)
k = reshape(kernel_conv(N, kern_G, kern_G, params_G, params_G, x, []), 2*n+1, 2*n+1);
k = k(2:end-1,2:end-1);
err_op_x = 1 - corr(k(:), xc_op_x(:));

% OP_y = (G*G)
k = reshape(kernel_conv(N, kern_G, kern_G, params_G, params_G, x, []), 2*n+1, 2*n+1);
k = k(2:end-1,2:end-1);
err_op_y = 1 - corr(k(:), xc_op_y(:));

% Total error
err = err_op_x + err_op_y;

end
