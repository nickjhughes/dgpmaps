
function k = kernel_dog(alpha, sigma, x, y)
%KERNEL_DOG A difference-of-Gaussians kernel function.
%
% k = kernel_dog(alpha, sigma, x, y)
%
% Evaluates a difference-of-Gaussians kernel function with parameters
% alpha and sigma. x is a set of coordinates to evaluate the kernel at,
% and y is an optional second set of coordinates, which if not given the
% origin we be used for.

if isempty(y)
    y = zeros(size(x,1),size(x,2));
end

tausqr = sum((x-y).^2,1);
k = (alpha^2)*(exp(-tausqr/(4*sigma^2)) + ...
         (1/4)*exp(-tausqr/(16*sigma^2)) - ...
         (4/5)*exp(-tausqr/(10*sigma^2)));
