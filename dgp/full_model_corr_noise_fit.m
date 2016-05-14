
function [post_mu, noise_model] = full_model_corr_noise_fit(data, Y, ntrials, alpha, alphaprime, sigma, sigmaprime, beta, gamma, rho, lambda)
%FULL_MODEL_CORR_NOISE_FIT GP estimate of OP and OD maps from given data.
%
% [post_mu, noise_model] = full_model_corr_noise_fit(data, Y, ntrials, alpha, alphaprime, sigma, sigmaprime, beta, gamma, rho, lambda)
%
% Get the posterior mean of the GP estimate of the OD and OP maps from the
% given data array (see vector_avg.m), using the given hyperparameters. To
% estimate parameters, use full_model_params.m.

% Settings
n = size(data, 1);
N1 = n^2;
N2 = N1;
N3 = N1;
N = N1 + N2 + N3;
[x, y] = meshgrid(1:n, 1:n);
x = [x(:) y(:)];

% Parameters
if nargin < 4
    alpha = 1;
end
if nargin < 5
    alphaprime = 1;
end
if nargin < 6
    sigma = 2;
end
if nargin < 7
    sigmaprime = 2;
end
if nargin < 8
    beta = 0.2;
end
if nargin < 9
    gamma = 0.75;
end
if nargin < 10
    rho = 2;
end
if nargin < 11
    lambda = 1.2;
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

% Check for no noise condition
no_noise = false;
noise_iterations = 3;
min_noise = 1e-3;
noise_var_est = mean(reshape(var(data,[],5),n*n*size(data,3)*size(data,4),1));
if noise_var_est < min_noise
    fprintf('No noise.\n');
    no_noise = true;
    noise_iterations = 1;
    noise_model = struct();
    noise_cov = min_noise*eye(N1,N1);
end

% Learn noise, derive posterior loop
for jj = 1:noise_iterations
    if ~no_noise
        % Learn noise model
        if jj == 1
            noise_model = dgp_learn_noise(data, [], []);
        else
            noise_model = dgp_learn_noise(data, [], [], post_mu);
        end
        noise_cov = diag(noise_model.D) + noise_model.G*noise_model.G.';
    end
    
    % Construct observation prior = C + noise covariance
    Cs11 = C11 + (1/ntrials)*noise_cov;
    Cs12 = C12;
    Cs13 = C13;
    Cs21 = C21;
    Cs22 = C22 + (1/ntrials)*noise_cov;
    Cs23 = C23;
    Cs31 = C31;
    Cs32 = C32;
    Cs33 = C33 + (1/ntrials)*noise_cov;
    Cs = [Cs11, Cs12', Cs13'; Cs21, Cs22, Cs23'; Cs31, Cs32, Cs33];
    
    % Calculate joint posterior
    CsdivY = Cs\Y;
    post_mu = nan(N1,3);
    for l = 1:N1
        k = nan(N,3);
        for j = 1:N1
            k(j,1) = C11(l,j);
            k(j,2) = C21(l,j);
            k(j,3) = C31(l,j);
        end
        for j = N1+1:N1+N2
            k(j,1) = C12(j-N1,l);
            k(j,2) = C22(l,j-N1);
            k(j,3) = C32(l,j-N1);
        end
        for j = N1+N2+1:N1+N2+N3
            k(j,1) = C13(j-N1-N2,l);
            k(j,2) = C23(j-N1-N2,l);
            k(j,3) = C33(l,j-N1-N2);
        end
        post_mu(l,1) = k(:,1).'*CsdivY;
        post_mu(l,2) = k(:,2).'*CsdivY;
        post_mu(l,3) = k(:,3).'*CsdivY;
    end    
end
