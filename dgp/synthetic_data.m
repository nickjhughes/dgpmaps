
function [data, od, op, all_noise] = synthetic_data(n, ntrials, seed, sigma_noise, sigma_corr_noise, eyes, orientations, od, op)
%SYNTHETIC_DATA Generate synthetic imaging data.
%
% [data, od, op, all_noise] = synthetic_data(n, ntrials, seed, sigma_noise, sigma_corr_noise, eyes, orientations, od, op)
%
% Generates synthetic experimental imaging data from either the given
% OP and OD map pair, or from a pair of maps sampled from the joint prior.
% n is map size, ntrials the number of imaging trials per stimulus
% condition, seed the RNG seed to use for added noise and sampling maps.

% Inputs
if nargin < 5
    sigma_corr_noise = 0;
end
if nargin < 6
    eyes = [-1, 1];
end
if nargin < 7
    orientations = [0, 45, 90, 135];
end
if nargin < 8
    sample_maps = true;
else
    sample_maps = false;
end

if sample_maps
    % Sample maps
    [od, op] = full_model_sample(n, seed);
else
    % Or just seed the RNG
    rng(seed);
end

% Normalise maps
norm_od = (od-mean(od(:)))/std(od(:)-mean(od(:)));
norm_op = (real(op)-mean(real(op(:))))/std(real(op(:))-mean(real(op(:)))) + ...
       1i*(imag(op)-mean(imag(op(:))))/std(imag(op(:))-mean(imag(op(:))));

% Generate synthetic imaging responses
sing_conds = nan(n, n, length(eyes), length(orientations));
for j = 1:length(eyes)
    for k = 1:length(orientations)
        sing_conds(:,:,j,k) = eyes(j)*norm_od + ...
            (cosd(2*orientations(k))*real(norm_op) + ...
            sind(2*orientations(k))*imag(norm_op));
    end
end

if sigma_corr_noise == 0
    % Add independent noise
    data = repmat(sing_conds, [1,1,1,1,ntrials]);
    all_noise = sqrt(sigma_noise)*randn(size(data));
    data = data + all_noise;
else
    % Add correlated noise
    noise_dim = 2;
    noise_filter = lowpass_filter(50, 0.5, 0.5);
    D = abs(conv2(randn(n,n), noise_filter, 'same'));
    D = D(:)/sqrt(mean(D(:).^2))*sigma_noise;
    G = nan(n*n, noise_dim);
    for j = 1:noise_dim
        l = conv2(randn(n,n), noise_filter, 'same');
        l = l(:)/sqrt(mean(l(:).^2))*sigma_corr_noise/sqrt(noise_dim);
        G(:,j) = l(:);
    end
    noise_std = reshape(sqrt(D), n, n);
    data = nan(n, n, length(eyes), length(orientations), ntrials);
    all_noise = nan(size(data));
    for j = 1:length(eyes)
        for k = 1:length(orientations)
            for l = 1:ntrials
                noise = randn(n, n);
                common_input = randn(noise_dim, 1);
                noise = noise.*noise_std + reshape(G*common_input, size(noise));
                all_noise(:,:,j,k,l) = noise;
                data(:,:,j,k,l) = sing_conds(:,:,j,k) + noise;
            end
        end
    end
    
end

end


function A = lowpass_filter(n, decay, p)
% Generate a low-pass filtering kernel.

[x, y] = meshgrid(-n:n, -n:n);
r = sqrt(x.^2 + y.^2);
A = exp(-(r/decay).^p);
A = A/sum(A(:));

end
