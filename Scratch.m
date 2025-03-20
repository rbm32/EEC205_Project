clc; clear; close all;

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

%% --- Perform Kirchhoff Migration ---
epsilon_r = 1;
c = 3e8;

x_img = linspace(-1, 1, 50);  % Imaging grid x-coordinates
z_img = linspace(0, 2, 50);   % Imaging grid z-coordinates
BP_image = kirchhoff_migration(x, ground_z, x_img, z_img, epsilon_r, c, bscan);

%% --- Plot B-Scan Data ---
figure;
imagesc(x_ant, time, bscan);
colormap gray; colorbar;
xlabel('Antenna Position (m)');
ylabel('Time (s)');
title('Simulated B-Scan');

%% --- Plot GPR Image from Migration ---
figure;
imagesc(x_img, -z_img, BP_image'); % Flip z-axis for correct visualization
colormap jet; colorbar;
xlabel('Horizontal Distance (m)');
ylabel('Depth (m)');
title('GPR Image using Fast Kirchhoff Migration');
set(gca, 'YDir', 'normal');


function BP_image = kirchhoff_migration(x_ant, h, x_img, z_img, epsilon_r, c, signal_data)
    % Implements the fast back-projection algorithm with cross-correlation for GPR imaging.
    %
    % Inputs:
    %   x_ant      - Antenna positions (1xM array) (m)
    %   h          - Antenna height above ground (m)
    %   x_img      - Imaging grid x-coordinates (Nx1 array) (m)
    %   z_img      - Imaging grid z-coordinates (1xZ array) (m)
    %   epsilon_r  - Relative permittivity of the ground
    %   c          - Speed of light in vacuum (m/s)
    %   signal_data- GPR received signals (MxT matrix)
    %
    % Output:
    %   BP_image   - Reconstructed GPR image using Kirchhoff migration
    
    % Define speed of wave in ground medium
    v_ground = c / sqrt(epsilon_r);
    M = length(x_ant);  % Number of antenna positions
    N = length(x_img);  % Number of x-grid points
    Z = length(z_img);  % Number of z-grid points
    BP_image = zeros(N, Z);  % Initialize image

    % Lookup Table for Time Delays
    time_lookup = zeros(N, Z, M);

    % Precompute Time Delays Using Approximate Method
    for n = 1:N
        for z = 1:Z
            x0 = x_img(n);
            z0 = z_img(z);
            for m = 1:M
                xk = x_ant(m);

                % Compute Refraction Point xr using the improved approximation method
                x1 = xk + (x0 - xk) / sqrt(epsilon_r);
                if abs(xk - x0) < (z0 + h) * sqrt(epsilon_r / (epsilon_r - 1))
                    xr = x0 + (x1 - x0) / sqrt(epsilon_r);
                elseif xk >= x0 + (z0 + h) * sqrt(epsilon_r / (epsilon_r - 1))
                    xr = x0 + z0 / (sqrt(epsilon_r) - 1);
                else
                    xr = x0 - z0 / (sqrt(epsilon_r) - 1);
                end

                % Compute Time Delays
                d_air = sqrt(h^2 + (xr - xk)^2); % Air segment
                d_ground = sqrt(z0^2 + (x0 - xr)^2); % Ground segment
                
                % Total Time of Flight (ToF)
                time_lookup(n, z, m) = (2 * d_air / c) + (2 * d_ground / v_ground);
            end
        end
    end

    % Perform Kirchhoff Migration with Cross-Correlation
    for n = 1:N
        for z = 1:Z
            signal_sum = zeros(M, 1);
            for m = 1:M
                time_delay = time_lookup(n, z, m);
                t_index = round(time_delay * size(signal_data, 2) / max(time_lookup(:))); % Convert to index
                
                % Ensure index is within bounds
                if t_index > 0 && t_index <= size(signal_data, 2)
                    signal_sum(m) = signal_data(m, t_index);
                end
            end

            % Apply Cross-Correlation for Artifact Suppression
            for i = 1:M-1
                for j = i+1:M
                    BP_image(n, z) = BP_image(n, z) + signal_sum(i) * signal_sum(j);
                end
            end
        end
    end

    % Normalize Image
    BP_image = BP_image / max(abs(BP_image(:)));

    % Display Image
    figure;
    imagesc(x_img, -z_img, BP_image'); % Flip z-axis for correct visualization
    colormap jet; colorbar;
    xlabel('Horizontal Distance (m)');
    ylabel('Depth (m)');
    title('GPR Image using Fast Kirchhoff Migration');
    set(gca, 'YDir', 'normal');
end


