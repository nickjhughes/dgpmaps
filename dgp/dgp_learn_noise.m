
function output = dgp_learn_noise(data, eyes, orientations, post_mu)
%DGP_LEARN_NOISE Learn the noise model for the given data.
%
% output = dgp_learn_noise(data, eyes, orientations, post_mu)
%
% Uses factor analysis method to learn the low-rank noise model.
% Based on code writen by Macke et al. for the paper,
% "J. H. Macke, S. Gerwinn, L. White, M. Kaschube, and M. Bethge
% Gaussian process methods for estimating cortical maps
% NeuroImage, 2010",
% available at: https://bitbucket.org/mackelab/gp_maps/src
% 
% See also:
% dgp_factor_analysis

% Inputs
if isempty(eyes)
    eyes = [-1, 1];
end
if isempty(orientations)
    orientations = [0, 45, 90, 135];
end

% Settings
G_scalar = 0.9;
q = 2;

% Initialisation
output = struct();
n = size(data,1);
ntrials = size(data,5);
have_post_mean = (nargin == 4);

if have_post_mean
    % Normalise posterior mean
    norm_post_mu = nan(size(post_mu));
    for j = 1:size(post_mu,2)
        norm_post_mu(:,j) = (post_mu(:,j) - mean(post_mu(:,j)))/std(post_mu(:,j) - mean(post_mu(:,j)));
    end
end

% Calculate residuals
res = nan(n*n,length(eyes),length(orientations),ntrials);
for j = 1:length(eyes)
    for k = 1:length(orientations)
        if have_post_mean
            if size(post_mu,2) == 2
                stimuli_mean = (cosd(2*orientations(k))*norm_post_mu(:,1) + ...
                    sind(2*orientations(k))*norm_post_mu(:,2));
            elseif size(post_mu,2) == 3
                stimuli_mean = eyes(j)*norm_post_mu(:,1) + ...
                    (cosd(2*orientations(k))*norm_post_mu(:,2) + ...
                    sind(2*orientations(k))*norm_post_mu(:,3));
            end
        else
            stimuli_mean = reshape(mean(data(:,:,j,k,:),5), n*n, 1);
        end
        res(:,j,k,:) = reshape(data(:,:,j,k,:),n*n,ntrials) - repmat(stimuli_mean,1,ntrials);
    end
end
res = reshape(res, n*n, length(eyes)*length(orientations)*ntrials);

% Default noise options
options = struct('q',q, 'max_iters',500, 'precision',0.01, 'ridge',0.02);
[G, ~, fadetails] = dgp_factor_analysis(res, options);

G = G_scalar*G;
wvar = zeros(size(G,1),1);
for k = 1:size(G,1)
    wvar(k,1) = G(k,:)*G(k,:)';
end
output.D = fadetails.var_y - wvar;
output.D = output.D(:);

output.factor_analysis = fadetails;
output.G = G;
output.q = q;
