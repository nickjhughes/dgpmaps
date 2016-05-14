
function [G, D, details] = dgp_factor_analysis(y, options, G0, D0)
%DGP_FACTOR_ANALYSIS Factor analysis noise model learning.
%
% [G, D, details] = dgp_factor_analysis(y, options, G0, D0)
%
% Used to learn the noise model. See dgp_learn_noise.m for usage.
% Based on code writen by Macke et al. for the paper,
% "J. H. Macke, S. Gerwinn, L. White, M. Kaschube, and M. Bethge
% Gaussian process methods for estimating cortical maps
% NeuroImage, 2010",
% available at: https://bitbucket.org/mackelab/gp_maps/src
% 
% See also:
% dgp_learn_noise

% Initialisation
q = options.q;
max_iters = options.max_iters;
min_iters = 5;
precision = options.precision;
ridge = options.ridge;
scores = zeros(1,max_iters);
[n, N] = size(y); % n=pixels, N=trials

% Valdiate algorithm by comparing results to the true covariance of a small set of random pixels
test_size = min(200,n);
test_index = randperm(n,test_size);
y_test = y(test_index,:);
truecov = (y_test*y_test')/N;
truecorr = gp_cov2corr(truecov);
truecorr = truecorr(~eye(size(truecorr)));

% Initialise covariance
if nargin < 3
    G = randn(n,q)/sqrt(n);
else
    G = G0;
end

% Initialise variance
var_y = mean(y'.^2)' + ridge;
if nargin < 4
    D = var_y;
else
    D = D0;
end

details = struct();

for j = 1:max_iters
    Dinv = spdiags(D.^1,0,n,n);
    
    % Calculate score (which is marginal likelihood)
    scores(j) = -N*n/2*log(2*pi) - N/2*logdet(eye(q)+G(test_index,:)'*Dinv(test_index,test_index)*G(test_index,:)) - N/2*sum(log(D));
    scores(j) = scores(j) - 1/2*trace(y_test'*lr_left_divide(y_test, D(test_index), G(test_index,:), G(test_index,:)', eye(q)));
    scores(j) = -scores(j)/(N*n);
    
    % Compare to true covariance
    testcov = G(test_index,:)*G(test_index,:)' + diag(D(test_index));
    testcorr = gp_cov2corr(testcov);
    testcorr = testcorr(~eye(size(testcorr)));
    
    details.corr_corr(j) = corr(truecorr,testcorr);
    details.corr_mse(j) = mean((truecorr-testcorr).^2);
    details.cov_mse(j) = mean((truecov(~eye(size(truecov)))-testcov(~eye(size(testcov)))).^2);
    details.var_corr(j) = corr(diag(truecov),diag(testcov));
    details.var_mse(j) = mean((diag(truecov)-diag(testcov)).^2);
    a = polyfit(truecorr,testcorr,1);
    details.corr_scaling(j) = a(1);
    
    % Check stopping criteria (and that the losses aren't increasing)
    if j > min_iters && abs(scores(j)-scores(j-1)) < precision && ...
       abs(details.corr_scaling(j)-details.corr_scaling(j-1)) < precision && ...
       abs(details.cov_mse(j)-details.cov_mse(j-1)) < precision
        fprintf('Losses do not improve anymore, aborting at %g.\n', scores(j));
        break
    elseif j > 1 && scores(j) > scores(j-1) + 10*precision
        fprintf('Warning: Loss is actually increasing, from %g to %g.\n', scores(j-1), scores(j));
    end
    
    % Update D and G (as per equations in the paper)
    beta = (eye(q) + G'*Dinv*G)\(G'*Dinv);
    A = y*(y'*beta');
    B = beta*A;
    A = A/N;
    B = B/N;
    B = eye(size(B)) - beta*G + B;
    G = A/B;
    for k = 1:n
        D(k) = A(k,:)*G(k,:)';
    end
    D = var_y - D;
    if min(D) < ridge
        fprintf('Warning: Noise variances are small.\n');
    end
    D = max(D,ridge);
end

% Return other output
scores = scores(1:j);
details.truecov = truecov;
details.testcov = testcov;
details.var_y = var_y;
details.scores = scores;
