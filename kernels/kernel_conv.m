
function k = kernel_conv(n, kernel1, kernel2, params1, params2, x, y)
%KERNEL_CONV The numerical convolution of two kernels.
%
% k = kernel_conv(n, kernel1, kernel2, params1, params2, x, y)
%
% Evaluates the kernel function which results from the convolution of the
% two given kernels. n is the size of the grid on which to numerically
% evaluate the convolution. The kernel parameters are given as integers,
% following the definitions in kernel_types.m. Each kernel can also be
% given as a vector of two kernels itself, which are evaluated recursively,
% in which case the corresponding parameters should be given as a cell
% array of two vectors of parameters.
% 
% Requires convolve2:
% http://www.mathworks.com/matlabcentral/fileexchange/22619-fast-2-d-convolution

% Initialisation
kernel_types;
if isempty(y)
    y = zeros(size(x,1),size(x,2));
end
N = size(x,2);

% Convolve kernels recursively
[px, py] = meshgrid(-n:n, -n:n);
if length(kernel1) > 1
    K1 = reshape(kernel_conv(n, kernel1(1), kernel1(2:end), params1{1}, params1{2:end}, [px(:),py(:)]', []), size(px,1), size(px,2));
else
    switch kernel1
        case kernels.dog
            K1 = reshape(kernel_dog(params1.alpha, params1.sigma, [px(:),py(:)]', []), size(px,1), size(px,2));
        case kernels.gaborh
            K1 = reshape(kernel_gaborh(params1.beta, params1.gamma, params1.rho, params1.lambda, [px(:),py(:)]', []), size(px,1), size(px,2));
        case kernels.gaborv
            K1 = reshape(kernel_gaborv(params1.beta, params1.gamma, params1.rho, params1.lambda, [px(:),py(:)]', []), size(px,1), size(px,2));
        otherwise
            error('kernel_conv does not support that kernel.');
    end
end
if length(kernel2) > 1
    K2 = reshape(kernel_conv(n, kernel2(1), kernel2(2:end), params2{1}, params2{2:end}, [px(:),py(:)]', []), size(px,1), size(px,2));
else
    switch kernel2
        case kernels.dog
            K2 = reshape(kernel_dog(params2.alpha, params2.sigma, [px(:),py(:)]', []), size(px,1), size(px,2));
        case kernels.gaborh
            K2 = reshape(kernel_gaborh(params2.beta, params2.gamma, params2.rho, params2.lambda, [px(:),py(:)]', []), size(px,1), size(px,2));
        case kernels.gaborv
            K2 = reshape(kernel_gaborv(params2.beta, params2.gamma, params2.rho, params2.lambda, [px(:),py(:)]', []), size(px,1), size(px,2));
        otherwise
            error('kernel_conv does not support that kernel.');
    end
end
KK = convolve2(K1, K2, 'same');

% Look up covariances
k = nan(1,N);
for j = 1:N
    dx = x(:,j) - y(:,j);
    d = n + dx + 1;
    k(j) = KK(d(2), d(1));
end
