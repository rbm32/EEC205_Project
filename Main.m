clc; clear;

addpath("Functions\");
%% --- Parameters ---
num_scans = 50;
num_time_samples = 50;
x = linspace(-3, 3, num_scans); % Antenna positions
z = linspace(0, 2, num_time_samples);

% Targets: [x, z, radius, reflectivity]
targets = [
    -1.5, .8, 0.5;
    1, .7, 0.5;
];

% Constants
c = 3e8;
t_max = 2e-8;
h = 0.25;
epsilon_r = 2;
alpha = 0.001;
t = linspace(0, t_max, num_time_samples);

% Gaussian pulse
pulse_sigma = 1e-9;        % Standard deviation (in same units as time)
pulse_duration = 6 * pulse_sigma;  % Total duration of pulse window (±3σ)

Ts = t(2) - t(1);
t_pulse = -pulse_duration/2 : Ts : pulse_duration/2;
pulse = exp(-t_pulse.^2 / (2 * pulse_sigma^2));

% Noise
noise_matrix = .002 * randn(num_time_samples, num_scans);

% --- Create Reflectivity Matrix g ---
g = create_reflectivity_matrix(x, z, targets);
clf;
imagesc(x, z, g);
colormap(gray); colorbar;
xlabel('x (m)');
ylabel('z (m)');
title('Reflectivity Matrix (g)');
set(gca, 'YDir', 'reverse'); % Depth increases downward

% --- Generate Artificial B-scan ---
[H, t, x_vec, X, Z] = generate_forward_model_matrix(x, z, num_time_samples, pulse, c, t_max, epsilon_r, h, alpha);
bscan = H * g(:) + noise_matrix(:);
bscan = reshape(bscan, num_time_samples, length(x_vec));  % Reshape into B-scan
plot_g_and_bscan(g, bscan, x, t, z, -12,"Subsurface Reflectivity", "Simulated B-Scan");

%% --- Standard Backprojection ---
tic;
reconstructed = BScanBackprojection(bscan, x, t, c, h, epsilon_r, num_scans);
ExecutionTime_BP = toc;

plot_gpr_results(g, bscan, reconstructed, x, t, num_scans, -3);
mse_BP = compute_mse(g, reconstructed)

%%
tic;
reconstructed_fast = BScanBackprojection_Fast(bscan, x, t, c, h, epsilon_r, num_scans);
ExecutionTime_BPF = toc;

plot_gpr_results(g, bscan, reconstructed_fast, x, t, num_scans, -12);
mse_BPF = compute_mse(g, reconstructed_fast)

%% --- Self-Correlation Backprojection ---
correlation_threshold = 0.9;
tic;
reconstructed_sbp = self_correlation_backprojection_vectorized(bscan, x, t, c, h, epsilon_r, num_scans, alpha, correlation_threshold);
executionTime_SBP = toc;

plot_gpr_results(g, bscan, reconstructed_sbp, x, t, num_scans, -3);
mse_SBP = compute_mse(g, reconstructed_sbp);
%% --- Fast Cross-Correlation Backprojection ---
correlation_threshold = 0.9;
tic;
reconstructed_fsbp = fast_bp_crosscorr(bscan, x, t, c, h, epsilon_r, num_scans, alpha, correlation_threshold);
executionTime_FSBP = toc;

plot_gpr_results(g, bscan, reconstructed_fsbp, x, t, num_scans, -3);
mse_FSBP = compute_mse(g, reconstructed_sbp);
%% --- Bayesian Reconstruction ---
lambda = 0.01;
tic;
reconstructed_bay = reconstruct_g_bayesian(H, bscan, num_scans, num_time_samples, lambda);
ExecutionTime_BAYES = toc;

mse_BR= compute_mse(g, reconstructed_bay);

plot_gpr_results(g, bscan, reconstructed_bay, x, t, num_scans, -6);
%% Real Data
addpath("Data");
raw = csvread("Data\ryanPulse_v1.csv");
meta = csvread("Data\ryanPulse_v2-meta_x_t_times.csv");
xraw = meta(1,:);
xraw = xraw(xraw~=0);
traw = meta(2,:);
%% preprocessing
bscan_BKGR = svd_background_removal(raw, 1);
imagesc(xraw, traw, bscan_BKGR);
x_range = [-inf inf];
t_range = [3  7] *1e-9;
[bscan_cropped, x_cropped, t_cropped] = crop_bscan(bscan_BKGR, xraw, traw, x_range, t_range);
imagesc(x_cropped, t_cropped, bscan_cropped);

time_factor = 5;
x_factor = 5;
[bscan_ds, t_ds, x_ds] = downsample_bscan(bscan_cropped, t_cropped, x_cropped, time_factor, x_factor);


alpha = 3
bscan_weighted = compensate_attenuation(bscan_ds, t_ds, c, epsilon_r, alpha)

bscan = bscan_weighted;
t = t_ds;
x = x_ds;
disp(numel(bscan));
imagesc(x, t, bscan);


%%
pulse_sigma = 1e-9;        % Standard deviation (in same units as time)
pulse_duration = 6 * pulse_sigma;  % 

num_scans = length(x);
num_time_samples = length(t);

Ts = t(2) - t(1);
t_pulse = -pulse_duration/2 : Ts : pulse_duration/2;
pulse = exp(-t_pulse.^2 / (2 * pulse_sigma^2));

z = linspace(0, .5, num_time_samples);

epsilon_r = 2.1
alpha = 0;
h = .2;
c = 3e8;
t_max = 2e-8;

%% System Matrix
[Hr, time, x_vec, X, Z] = generate_forward_model_matrix_fast(x, z, num_time_samples, pulse, c, t_max, epsilon_r, h, alpha);

%% --- Reconstruction using All Methods ---

lambda = 0.01;
tic;
reconstructed_bay = reconstruct_g_bayesian(Hr, bscan, num_scans, num_time_samples, lambda);
ExecutionTime_Bayes = toc;
plot_g_and_bscan(reconstructed_bay, bscan, x, t, z, -11);



