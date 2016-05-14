
function [od, op] = full_model_sample(n, seed, alpha, alphaprime, sigma, sigmaprime, beta, gamma, rho, lambda, ridge)
%FULL_MODEL_SAMPLE Sample a pair of OD and OP maps from the full model.
%
% [od, op] = full_model_sample(n, seed, alpha, alphaprime, sigma, sigmaprime, beta, gamma, rho, lambda, ridge)
%
% n is map size, seed is RNG seed, greek letters are hyperparameters,
% ridge is the value to add to the diagonal of the prior covariance matrix
% for numerical stability.

% Settings
if nargin < 1
    n = 16;
end
N1 = n^2;
N2 = N1;
N3 = N1;
N = N1 + N2 + N3;
[x, y] = meshgrid(1:n, 1:n);
x = [x(:) y(:)];
if nargin < 2
    seed = 1;
end

% Parameters
if nargin < 3
    alpha = 1;
end
if nargin < 4
    alphaprime = 1;
end
if nargin < 5
    sigma = 2;
end
if nargin < 6
    sigmaprime = 2;
end
if nargin < 7
    beta = 0.2;
end
if nargin < 8
    gamma = 0.75;
end
if nargin < 9
    rho = 2;
end
if nargin < 10
    lambda = 1.2;
end
if nargin < 11
    ridge = 1e-10;
end

% Construct covariance function lookups
[px, qy] = meshgrid(-n:n, -n:n);
kern_dog = reshape(kernel_dog(alpha,sigma,[px(:),qy(:)]',[]),size(px,1),size(px,2));
kern_dog_prime = reshape(kernel_dog(alphaprime,sigmaprime,[px(:),qy(:)]',[]),size(px,1),size(px,2));
kern_gabor_h = reshape(kernel_gaborh(beta,gamma,rho,lambda,[px(:),qy(:)]',[]),size(px,1),size(px,2));
kern_gabor_v = reshape(kernel_gaborv(beta,gamma,rho,lambda,[px(:),qy(:)]',[]),size(px,1),size(px,2));
kern_gabor_h_dog = conv2(kern_dog, kern_gabor_h, 'same');
kern_gabor_v_dog = conv2(kern_dog, kern_gabor_v, 'same');
C_11 = conv2(kern_dog, kern_dog, 'same');
C_12 = -conv2(kern_dog, kern_gabor_h_dog, 'same');
C_13 = -conv2(kern_dog, kern_gabor_v_dog, 'same');
C_21 = -conv2(kern_gabor_h_dog, kern_dog, 'same');
C_22 = -conv2(kern_gabor_h_dog, kern_gabor_h_dog, 'same');
C_23 = -conv2(kern_gabor_h_dog, kern_gabor_v_dog, 'same');
C_31 = -conv2(kern_gabor_v_dog, kern_dog, 'same');
C_32 = -conv2(kern_gabor_v_dog, kern_gabor_h_dog, 'same');
C_33 = -conv2(kern_gabor_v_dog, kern_gabor_v_dog, 'same');
C_22_ = conv2(kern_dog_prime, kern_dog_prime, 'same');
C_33_ = conv2(kern_dog_prime, kern_dog_prime, 'same');

% Construct prior covariance matrix
C11 = nan(N1,N1);
C12 = nan(N1,N2);
C13 = nan(N1,N3);
C21 = nan(N2,N1);
C22 = nan(N2,N2);
C23 = nan(N2,N3);
C31 = nan(N3,N1);
C32 = nan(N3,N2);
C33 = nan(N3,N3);
for j = 1:N1
    for k = 1:N1
        d = n+(x(j,:)-x(k,:));
        C11(j,k) = C_11(d(1)+1, d(2)+1);
    end
end
for j = 1:N1
    for k = 1:N2
        d = n+(x(j,:)-x(k,:));
        C12(j,k) = C_12(d(1)+1, d(2)+1);
    end
end
for j = 1:N1
    for k = 1:N3
        d = n+(x(j,:)-x(k,:));
        C13(j,k) = C_13(d(1)+1, d(2)+1);
    end
end

for j = 1:N2
    for k = 1:N1
        d = n+(x(j,:)-x(k,:));
        C21(j,k) = C_21(d(1)+1, d(2)+1);
    end
end
for j = 1:N2
    for k = 1:N2
        d = n+(x(j,:)-x(k,:));
        C22(j,k) = C_22(d(1)+1, d(2)+1) + C_22_(d(1)+1, d(2)+1);
    end
end
for j = 1:N2
    for k = 1:N3
        d = n+(x(j,:)-x(k,:));
        C23(j,k) = C_23(d(1)+1, d(2)+1);
    end
end

for j = 1:N3
    for k = 1:N1
        d = n+(x(j,:)-x(k,:));
        C31(j,k) = C_31(d(1)+1, d(2)+1);
    end
end
for j = 1:N3
    for k = 1:N2
        d = n+(x(j,:)-x(k,:));
        C32(j,k) = C_32(d(1)+1, d(2)+1);
    end
end
for j = 1:N3
    for k = 1:N3
        d = n+(x(j,:)-x(k,:));
        C33(j,k) = C_33(d(1)+1, d(2)+1) + C_33_(d(1)+1, d(2)+1);
    end
end
C = [C11, C12', C13'; C21, C22, C23'; C31, C32, C33];
C = C + ridge*eye(size(C));

% Sample functions
rng(seed);
Y = randn(1,N)*chol(C);
Y1 = Y(1:N1);
Y2 = Y(N1+1:N1+N2);
Y3 = Y(N1+N2+1:end);
od = reshape(Y1,n,n);
op = reshape(Y2,n,n) + 1i*reshape(Y3,n,n);
