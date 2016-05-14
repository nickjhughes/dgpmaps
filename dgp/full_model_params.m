
function params = full_model_params(od, op, a0)
%FULL_MODEL_PARAMS Estimate hyperparameters using the given vector averaged maps.
%
% params = full_model_params(od, op, a0)
%
% Estimates parameters by maximising the correlation between the empirical
% and estimated prior kernel functions. params is a struct with fields
% alpha, sigma, alphaprime, sigmaprime, beta, gamma, rho, and lambda.
% a0 is a vector of the initial guesses for each hyperparameter, in the
% same order.

% Initialise
kernel_types;
n = size(od,1);
op_x = real(op);
op_y = imag(op);

% Calculate covariances
xc_od = xcorr2(od);
xc_od = xc_od/max(xc_od(:));
xc_op_x = xcorr2(op_x);
xc_op_x = xc_op_x/max(xc_op_x(:));
xc_op_y = xcorr2(op_y);
xc_op_y = xc_op_y/max(xc_op_y(:));

% Fit
if nargin < 3
    % alpha, sigma, alphaprime, sigmaprime, beta, gamma, rho, lambda
    a0 = [1, 2, 1, 2, 0.2, 0.75, 2, 1.2];
end
options = optimset('Display', 'Iter');
a = fminsearch(@(a)(fullfit(a, n, xc_od, xc_op_x, xc_op_y)), a0, options);
params = struct('alpha', a(1), ...
                'sigma', a(2), ...
                'alphaprime', a(3), ...
                'sigmaprime', a(4), ...
                'beta', a(5), ...
                'gamma', a(6), ...
                'rho', a(7), ...
                'lambda', a(8));

end


function err = fullfit(a, n, xc_od, xc_op_x, xc_op_y)

% Initialise
kernel_types;
[xx, yy] = meshgrid(-n:n, -n:n);
x = [yy(:) xx(:)]';
N = 400;

% Parameters
params_G = struct('alpha',a(1), 'sigma',a(2));
kern_G = kernels.dog;
params_Gprime = struct('alpha',a(3), 'sigma',a(4));
kern_Gprime = kernels.dog;
params_H = struct('beta',a(5), 'gamma',a(6), 'rho',a(7), 'lambda',a(8));
params_V = params_H;
kern_H = kernels.gaborh;
kern_V = kernels.gaborv;

% OD = (G*G)
k = reshape(kernel_conv(N, kern_G, kern_G, params_G, params_G, x, []), 2*n+1, 2*n+1);
k = k(2:end-1,2:end-1);
err_od = 1 - corr(k(:), xc_od(:));

% OP_x = (G*H)*(G*H) + (G'*G')
k1 = reshape(kernel_conv(N, [kern_G, kern_H], [kern_G, kern_H], ...
    {params_G, params_H}, {params_G, params_H}, x, []), 2*n+1, 2*n+1);
k1 = k1(2:end-1,2:end-1);
k2 = reshape(kernel_conv(N, kern_Gprime, kern_Gprime, ...
    params_Gprime, params_Gprime, x, []), 2*n+1, 2*n+1);
k2 = k2(2:end-1,2:end-1);
k = k1 + k2;
err_op_x = 1 - corr(k(:), xc_op_x(:));

% OP_y = (G*V)*(G*V) + (G'*G')
k1 = reshape(kernel_conv(N, [kern_G, kern_V], [kern_G, kern_V], ...
    {params_G, params_V}, {params_G, params_V}, x, []), 2*n+1, 2*n+1);
k1 = k1(2:end-1,2:end-1);
k2 = reshape(kernel_conv(N, kern_Gprime, kern_Gprime, ...
    params_Gprime, params_Gprime, x, []), 2*n+1, 2*n+1);
k2 = k2(2:end-1,2:end-1);
k = k1 + k2;
err_op_y = 1 - corr(k(:), xc_op_y(:));

% Total error
err = err_od + err_op_x + err_op_y;

end
