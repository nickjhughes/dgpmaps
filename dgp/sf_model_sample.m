
function [od, sf, op] = sf_model_sample(n, seed, alpha, alphaprime, alphasf, sigma, sigmaprime, sigmasf, beta, gamma, rho, lambda, ridge)
%SF_MODEL_SAMPLE Sample a triplet of OD, OP, SF maps from the SF model.
%
% [od, sf, op] = sf_model_sample(n, seed, alpha, alphaprime, alphasf, sigma, sigmaprime, sigmasf, beta, gamma, rho, lambda, ridge)
%
% n is map size, seed is RNG seed, greek letters are hyperparameters,
% ridge is the value to add to the diagonal of the prior covariance matrix
% for numerical stability.

% Parameters
if nargin < 1 || isempty(n)
    n = 16;
end
if nargin < 2 || isempty(seed)
    seed = 1;
end
if nargin < 3 || isempty(alpha)
    alpha = 1;
end
if nargin < 4 || isempty(alphaprime)
    alphaprime = 0.25;
end
if nargin < 5 || isempty(alphasf)
    alphasf = 0.25;
end
if nargin < 6 || isempty(sigma)
    sigma = 2;
end
if nargin < 7 || isempty(sigmaprime)
    sigmaprime = 2;
end
if nargin < 8 || isempty(sigmasf)
    sigmasf = 2*sigma;
end
if nargin < 9 || isempty(beta)
    beta = 0.2;
end
if nargin < 10 || isempty(gamma)
    gamma = 0.75;
end
if nargin < 11 || isempty(rho)
    rho = 2;
end
if nargin < 12 || isempty(lambda)
    lambda = 1.2;
end
if nargin < 13 || isempty(ridge)
    ridge = 1e-6;
end

% Setup
Ni = n^2;
N = 4*Ni;
[x, y] = meshgrid(1:n, 1:n);
x = [x(:) y(:)];

% Construct covariance function lookups
[px, qy] = meshgrid(-n:n, -n:n);
kern_dog = reshape(kernel_dog(alpha,sigma,[px(:),qy(:)]',[]),size(px,1),size(px,2));
kern_dog_prime = reshape(kernel_dog(alphaprime,sigmaprime,[px(:),qy(:)]',[]),size(px,1),size(px,2));
kern_dog_sf = reshape(kernel_dog(alphasf,sigmasf,[px(:),qy(:)]',[]),size(px,1),size(px,2));
kern_gabor_h = reshape(kernel_gaborh(beta,gamma,rho,lambda,[px(:),qy(:)]',[]),size(px,1),size(px,2));
kern_gabor_v = reshape(kernel_gaborv(beta,gamma,rho,lambda,[px(:),qy(:)]',[]),size(px,1),size(px,2));
kern_gabor_h_dog = conv2(kern_dog, kern_gabor_h, 'same');
kern_gabor_v_dog = conv2(kern_dog, kern_gabor_v, 'same');
% OD
C_11 = conv2(kern_dog, kern_dog, 'same');
C_12 = conv2(kern_dog, kern_dog_sf, 'same');
C_13 = -conv2(kern_dog, kern_gabor_h_dog, 'same');
C_14 = -conv2(kern_dog, kern_gabor_v_dog, 'same');
% SF
C_21 = conv2(kern_dog_sf, kern_dog, 'same');
C_22 = conv2(kern_dog_sf, kern_dog_sf, 'same');
C_23 = -conv2(kern_dog_sf, kern_gabor_h_dog, 'same');
C_24 = -conv2(kern_dog_sf, kern_gabor_v_dog, 'same');
% real(OP)
C_31 = -conv2(kern_gabor_h_dog, kern_dog, 'same');
C_32 = -conv2(kern_gabor_h_dog, kern_dog_sf, 'same');
C_33 = -conv2(kern_gabor_h_dog, kern_gabor_h_dog, 'same');
C_33_ = conv2(kern_dog_prime, kern_dog_prime, 'same');
C_34 = -conv2(kern_gabor_h_dog, kern_gabor_v_dog, 'same');
% imag(OP)
C_41 = -conv2(kern_gabor_v_dog, kern_dog, 'same');
C_42 = -conv2(kern_gabor_v_dog, kern_dog_sf, 'same');
C_43 = -conv2(kern_gabor_v_dog, kern_gabor_h_dog, 'same');
C_44 = -conv2(kern_gabor_v_dog, kern_gabor_v_dog, 'same');
C_44_ = conv2(kern_dog_prime, kern_dog_prime, 'same');

% Construct prior covariance matrix
C11 = nan(Ni,Ni);
C12 = nan(Ni,Ni);
C13 = nan(Ni,Ni);
C14 = nan(Ni,Ni);
C21 = nan(Ni,Ni);
C22 = nan(Ni,Ni);
C23 = nan(Ni,Ni);
C24 = nan(Ni,Ni);
C31 = nan(Ni,Ni);
C32 = nan(Ni,Ni);
C33 = nan(Ni,Ni);
C34 = nan(Ni,Ni);
C41 = nan(Ni,Ni);
C42 = nan(Ni,Ni);
C43 = nan(Ni,Ni);
C44 = nan(Ni,Ni);
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C11(j,k) = C_11(d(1)+1, d(2)+1);
    end
end
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C12(j,k) = C_12(d(1)+1, d(2)+1);
    end
end
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C13(j,k) = C_13(d(1)+1, d(2)+1);
    end
end
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C14(j,k) = C_14(d(1)+1, d(2)+1);
    end
end

for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C21(j,k) = C_21(d(1)+1, d(2)+1);
    end
end
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C22(j,k) = C_22(d(1)+1, d(2)+1);
    end
end
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C23(j,k) = C_23(d(1)+1, d(2)+1);
    end
end
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C24(j,k) = C_24(d(1)+1, d(2)+1);
    end
end

for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C31(j,k) = C_31(d(1)+1, d(2)+1);
    end
end
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C32(j,k) = C_32(d(1)+1, d(2)+1);
    end
end
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C33(j,k) = C_33(d(1)+1, d(2)+1) + C_33_(d(1)+1, d(2)+1);
    end
end
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C34(j,k) = C_34(d(1)+1, d(2)+1);
    end
end

for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C41(j,k) = C_41(d(1)+1, d(2)+1);
    end
end
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C42(j,k) = C_42(d(1)+1, d(2)+1);
    end
end
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C43(j,k) = C_43(d(1)+1, d(2)+1);
    end
end
for j = 1:Ni
    for k = 1:Ni
        d = n+(x(j,:)-x(k,:));
        C44(j,k) = C_44(d(1)+1, d(2)+1) + C_44_(d(1)+1, d(2)+1);
    end
end

C = [C11, C12', C13', C14'; C21, C22, C23', C24'; C31, C32, C33, C34; C41, C42, C43, C44];
C = C + ridge*eye(size(C));

% Sample functions from prior
rng(seed);
Y = randn(1,N)*chol(C);
Y1 = Y(1:Ni);
Y2 = Y(Ni+1:Ni+Ni);
Y3 = Y(Ni+Ni+1:Ni+Ni+Ni);
Y4 = Y(Ni+Ni+Ni+1:end);

% Feature maps
od = reshape(Y1,n,n);
sf = exp(reshape(Y2,n,n));
op = reshape(Y3,n,n) + 1i*reshape(Y4,n,n);
