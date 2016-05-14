
%% Add required directories to path

addpath(strcat(pwd, '\dgp'));
addpath(strcat(pwd, '\kernels'));
addpath(strcat(pwd, '\mapstats'));
addpath(strcat(pwd, '\plotting'));
addpath(strcat(pwd, '\va'));

% 3rd party libraries - see README.md
addpath(strcat(pwd, '\convolve2'));
addpath(strcat(pwd, '\export_fig'));
addpath(strcat(pwd, '\intersections'));

%% Plotting settings

line_width = 2;
axis_width = 1;
font_size = 10;
figure_width = 8;
figure_height = (3/4)*figure_width;
cmap_op = hsv(16);
cmap_od = gray(17);
cmap_od = cmap_od(1:16,:);
cmap_sf = jet(16);
cmap_noise = gray(128);
cmap_filter = jet(16);

%% Schematic filters

% Noise
n = 50;
rng(1);
X = randn(n);
plot_od(X, 'schematic_X0.pdf', cmap_noise, ...
    figure_width, figure_height);

% G
n = 15;
dn = 1;
alpha = 1;
sigma = 2;
[x, y] = meshgrid(-n:dn:n, -n:dn:n);
G = nan(size(x,1),size(x,2));
for j = 1:size(x,1)
    for k = 1:size(x,2)
        G(j,k) = kernel_dog(alpha, sigma, ...
            [x(j,k), y(j,k)].', []);
    end
end
plot_filter(G, 'schematic_G.pdf', cmap_filter, ...
    figure_width, figure_height);

% G'
dn = 1;
n = 15;
alpha = 0.5;
sigma = 1.5;
[x, y] = meshgrid(-n:dn:n, -n:dn:n);
Gprime = nan(size(x,1),size(x,2));
for j = 1:size(x,1)
    for k = 1:size(x,2)
        Gprime(j,k) = kernel_dog(alpha, sigma, ...
            [x(j,k), y(j,k)].', []);
    end
end
plot_filter(Gprime, 'schematic_Gprime.pdf', cmap_filter, ...
    figure_width, figure_height);

% H
n = 15;
dn = 0.75;
beta = 0.2;
gamma = 0.75;
rho = 2;
lambda = 1.2;
[x, y] = meshgrid(-n:dn:n, -n:dn:n);
H = nan(size(x,1),size(x,2));
for j = 1:size(x,1)
    for k = 1:size(x,2)
        H(j,k) = -kernel_gaborh(beta, gamma, rho, lambda, ...
            [x(j,k), y(j,k)].', []);
    end
end
plot_filter(H, 'schematic_H.pdf', cmap_filter, ...
    figure_width, figure_height);

% V
n = 15;
dn = 0.75;
beta = 0.2;
gamma = 0.75;
rho = 2;
lambda = 1.2;
[x, y] = meshgrid(-n:dn:n, -n:dn:n);
V = nan(size(x,1),size(x,2));
for j = 1:size(x,1)
    for k = 1:size(x,2)
        V(j,k) = -kernel_gaborv(beta, gamma, rho, lambda, ...
            [x(j,k), y(j,k)].', []);
    end
end
plot_filter(V, 'schematic_V.pdf', cmap_filter, ...
    figure_width, figure_height);

%% Sampling maps from the prior

% Sample maps
n = 50;
[od, op] = full_model_sample(n);

% OD
plot_od(od, 'example_od.pdf', cmap_od, figure_width, figure_height);

% OP
plot_op(op, 'example_op.pdf', cmap_op, figure_width, figure_height);

%% Estimating maps with vector averaging

% Settings
n = 50;
ntrials = 10;
seed = 1;
sigma_noise = 6;

% Generate synthetic imaging data
[data, od, op] = synthetic_data(n, ntrials, seed, ...
    sigma_noise, 0.5*sigma_noise);
plot_od(od, 'true_od.pdf', cmap_od, figure_width, figure_height);
plot_op(op, 'true_op.pdf', cmap_op, figure_width, figure_height);

% Vector averaged maps
[va_od, va_op] = vector_avg(data);
plot_od(va_od, 'va_od.pdf', cmap_od, figure_width, figure_height);
plot_op(va_op, 'va_op.pdf', cmap_op, figure_width, figure_height);

% Optimally filtered vector averaged maps
[f_va_od, f_va_op] = optimal_filtering(od, op, va_od, va_op);
plot_od(f_va_od, 'f_va_od.pdf', cmap_od, figure_width, figure_height);
plot_op(f_va_op, 'f_va_op.pdf', cmap_op, figure_width, figure_height);

%% Estimating maps with full model

% Settings
n = 50;
ntrials = 10;
seed = 1;
sigma_noise = 6;

% Generate synthetic imaging data
data = synthetic_data(n, ntrials, seed, sigma_noise, 0.5*sigma_noise);

% Data
[va_od, va_op] = vector_avg(data);
Y = vertcat(va_od(:), real(va_op(:)), imag(va_op(:)));

% Estimate parameters
params = full_model_params(va_od, va_op);

% Derive posterior
[post_mu, noise_model] = full_model_corr_noise_fit(data, Y, ntrials, ...
    params.alpha, params.alphaprime, params.sigma, params.sigmaprime, ...
    params.beta, params.gamma, params.rho, params.lambda);
gp_od = reshape(post_mu(:,1), n, n);
gp_op = reshape(post_mu(:,2) + 1i*post_mu(:,3), n, n);
plot_od(gp_od, 'gp_od.pdf', cmap_od, figure_width, figure_height);
plot_op(gp_op, 'gp_op.pdf', cmap_op, figure_width, figure_height);

%% Estimating maps with independent model

% Settings
n = 50;
ntrials = 10;
seed = 1;
sigma_noise = 6;

% Generate synthetic imaging data
data = synthetic_data(n, ntrials, seed, sigma_noise, 0.5*sigma_noise);

% Data
[va_od, va_op] = vector_avg(data);
iY = vertcat(real(va_op(:)), imag(va_op(:)));

% Estimate parameters just for OP
iparams = imodel_params(va_op);

% Derive posterior just for OP
[ipost_mu, inoise_model] = imodel_corr_noise_fit(data, iY, ntrials, ...
    iparams.alpha, iparams.sigma);
igp_op = reshape(ipost_mu(:,1) + 1i*ipost_mu(:,2),n,n);
plot_op(igp_op, 'igp_op.pdf', cmap_op, figure_width, figure_height);

%% Sampling spatial frequency maps

% Sample maps
n = 50;
[od, sf, op] = sf_model_sample(n);

% OD
plot_od(od, 'spatial_frequency_example_od.pdf', cmap_od, ...
    figure_width, figure_height);

% SF
plot_sf(sf, 'spatial_frequency_example_sf.pdf', cmap_sf, ...
    figure_width, figure_height);

% OP
plot_op(op, 'spatial_frequency_example_op.pdf', cmap_op, ...
    figure_width, figure_height);
