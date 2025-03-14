clear; clc;

% Define parameters
image_size = 100;
num_time_samples = 256;

% Define linear antenna positions (e.g., moving from -1 to +1)
x = linspace(-1, 1, 100);

% Define target locations and properties [x, z, radius, reflectivity]
targets = [
    0, .8, 0.05, 0.5;  % Target 1 at (0, 0.8) with reflectivity 0.8
    0.5, 1.2, 0.05, 0.8; % Target 2 at (0.5, 1.2) with reflectivity 0.5
];

% Define ground surface (flat line at normalized z = 1.0 with width 0.6, and reflectivity 0.2)
ground_z = .2;      % Ground level (normalized)
ground_width = .2;  % Width of the ground surface in terms of x
ground_reflectivity = 0.2; % Reflectivity of the ground

% Custom pulse example (e.g., Gaussian pulse)
t_pulse = linspace(-0.5, 0.5, 50);
pulse = exp(-t_pulse.^2 / (2 * 0.1^2));

% Create the reflectivity matrix g from the target locations and ground surface
g = create_reflectivity_matrix(image_size, targets, ground_z, ground_width, ground_reflectivity);

c = 3e8; % Speed of wave (e.g., 0.3 m/ns for EM waves in soil)
t_max = 10e-9; % Maximum time (ns)

noise_matrix = 0.5 * randn(num_time_samples, image_size);
%
% Generate the artificial B-scan data
[bscan, time, x] = generate_artificial_bscan(image_size, num_time_samples, x, g, pulse, c, t_max, noise_matrix);

% Perform Filtered Backprojection to estimate the targets
% Backprojection
reconstructed_image = filtered_backprojection(bscan, x, time, c);

% Convert to dB
reconstructed_image_dB = 20 * log10(abs(reconstructed_image) + eps); % Avoid log(0)

% Apply cutoff
cutoff_dB = 1; % Cutoff 20 dB below the max value
max_dB = max(reconstructed_image_dB(:));
reconstructed_image_dB(reconstructed_image_dB < (max_dB - cutoff_dB)) = NaN;

% Plot reconstructed image
nexttile;
imagesc(linspace(min(x), max(x), 100), linspace(0, max(time) * 0.3 / 2, 100), reconstructed_image_dB);
xlabel('Antenna Position (m)');
ylabel('Depth (m)');
title('Reconstructed Image (dB)');
colormap('jet'); % Set reconstructed image colormap to jet
colorbar;
axis tight;

tl = nexttile;
imagesc(x, time, abs(bscan));
xlabel('Antenna Position (m)');
ylabel('Time (ns)');
title('Artificial B-scan (Custom Pulse)');

colorbar;
axis tight;
colormap(tl, 'gray'); % Set B-scan colormap to gray

