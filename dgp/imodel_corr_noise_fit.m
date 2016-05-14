
function [post_mu, noise_model] = imodel_corr_noise_fit(data, Y, ntrials, alpha, sigma)
%IMODEL_CORR_NOISE_FIT GP estimate of OP map from given data.
%
% [post_mu, noise_model] = imodel_corr_noise_fit(data, Y, ntrials, alpha, sigma)
%
% Get the posterior mean of the GP estimate of the OP map from the given
% data array (see vector_avg.m), using the given hyperparameters. To
% estimate parameters, use imodel_params.m.

% Settings
n = size(data, 1);
N1 = n^2;
N2 = N1;
N = N1 + N2;
[x, y] = meshgrid(1:n, 1:n);
x = [x(:) y(:)];

% Parameters
if nargin < 4
    alpha = 1;
end
if nargin < 5
    sigma = 2;
end

% Construct covariance function lookups
[px, qy] = meshgrid(-n:n, -n:n);
kern_dog = reshape(kernel_dog(alpha,sigma,[px(:),qy(:)]',[]),size(px,1),size(px,2));
C_11 = conv2(kern_dog, kern_dog, 'same');
C_22 = C_11;

% Construct prior covariance matrix
C11 = nan(N1,N1);
C12 = zeros(N1,N2);
C21 = zeros(N2,N1);
C22 = nan(N2,N2);
for j = 1:N1
    for k = 1:N1
        d = n+(x(j,:)-x(k,:));
        C11(j,k) = C_11(d(1)+1, d(2)+1);
    end
end
for j = 1:N2
    for k = 1:N2
        d = n+(x(j,:)-x(k,:));
        C22(j,k) = C_22(d(1)+1, d(2)+1);
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
    Cs21 = C21;
    Cs22 = C22 + (1/ntrials)*noise_cov;
    Cs = [Cs11, Cs12; Cs21, Cs22];
    
    % Calculate posterior
    CsdivY = Cs\Y;
    post_mu = nan(N1,2);
    for l = 1:N1
        k = nan(N,2);
        for j = 1:N1
            k(j,1) = C11(l,j);
            k(j,2) = C21(l,j);
        end
        for j = N1+1:N1+N2
            k(j,1) = C12(j-N1,l);
            k(j,2) = C22(l,j-N1);
        end
        post_mu(l,1) = k(:,1).'*CsdivY;
        post_mu(l,2) = k(:,2).'*CsdivY;
    end    
end
