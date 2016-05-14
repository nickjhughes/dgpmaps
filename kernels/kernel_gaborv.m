
function k = kernel_gaborv(beta, gamma, rho, lambda, x, y)
%KERNEL_GABORV A vertical Gabor kernel function.
%
% k = kernel_gaborv(beta, gamma, rho, lambda, x, y)
%
% Evaluates a vertical Gabor kernel function with parameters beta, gamma,
% rho and lambda. x is a set of coordinates to evaluate the kernel at,
% and y is an optional second set of coordinates, which if not given the
% origin we be used for.

if isempty(y)
    y = zeros(size(x,1),size(x,2));
end

k = beta*exp(-((x(1,:)-y(1,:)).^2+gamma^2*(x(2,:)-y(2,:)).^2)/(2*rho^2)).*sin((2*pi*(x(1,:)-y(1,:)))/lambda);
