function g = mwbp_reconstruction(bscan, x_vec, time, z_vec, c, epsilon_r, h, max_iter, k_thresh, refine_factor, min_resolution)
% Multi-Scale Weighted Back-Projection Imaging for GPR
%
% Inputs:
%   bscan         - [time x antenna] GPR B-scan data
%   x_vec         - Antenna positions (x-axis)
%   time          - Time vector for each A-scan
%   z_vec         - Depth values (z-axis)
%   c             - Speed of light in vacuum
%   epsilon_r     - Relative permittivity
%   h             - Antenna height above ground
%   max_iter      - Max number of refinement iterations
%   k_thresh      - Threshold coefficient for PTR detection (e.g. 0.5)
%   refine_factor - Subdivision factor for grid refinement (e.g. 3)
%   min_resolution- Minimum cell size to stop iterating (e.g. 0.01 m)
%
% Output:
%   g - Reconstructed image (reflectivity)

% Initial coarse grid
[x_grid, z_grid] = meshgrid(x_vec, z_vec);
v = c / sqrt(epsilon_r); % propagation speed

% Initial image (same size as coarse grid)
g = zeros(size(x_grid));

% Initialize PTR to full image
ptr_mask = true(size(g));

% Loop through refinement levels
for iter = 1:max_iter
    fprintf('--- Iteration %d ---\n', iter);
    
    % Determine resolution
    dx = abs(x_vec(2) - x_vec(1));
    dz = abs(z_vec(2) - z_vec(1));
    fprintf('Grid size: dx=%.3f m, dz=%.3f m\n', dx, dz);
    
    % Weighted Back-Projection on current PTR
    g_iter = zeros(size(g));
    weight_map = zeros(size(g));
    
    for i = 1:size(x_grid,1)
        for j = 1:size(x_grid,2)
            if ~ptr_mask(i,j)
                continue;
            end

            x_p = x_grid(i,j);
            z_p = z_grid(i,j);

            sum_val = 0;
            signal_vals = zeros(1, length(x_vec));

            num_antennas = size(bscan, 2);
            for m = 1:num_antennas
                x_a = x_vec(m);
                % Compute time-of-flight
                r = sqrt((x_a - x_p)^2 + (z_p + h)^2);
                t = 2 * r / v;

                % Interpolate signal at this time
                signal_vals(m) = interp1(time, bscan(:,m), t, 'linear', 0);
            end

            % Weighted factor
            mean_val = mean(signal_vals);
            var_val = var(signal_vals);

            if var_val == 0
                weight = 1;
            else
                weight = mean_val / var_val;
            end

            weight_map(i,j) = weight;
            g_iter(i,j) = weight * sum(signal_vals);
        end
    end

    % Add to result
    g = g_iter;

    % PTR Detection
    g_abs = abs(g);
    max_val = max(g_abs(:));
    ptr_mask = g_abs >= k_thresh * max_val;

    % Check resolution stop condition
    if dx <= min_resolution && dz <= min_resolution
        fprintf('Reached resolution limit. Stopping.\n');
        break;
    end

    % Refine grid in PTR
    if iter < max_iter
        [x_vec_fine, z_vec_fine] = refine_grid(x_vec, z_vec, ptr_mask, refine_factor);

        % Update grid for next iteration
        [x_grid, z_grid] = meshgrid(x_vec_fine, z_vec_fine);
        x_vec = x_vec_fine;
        z_vec = z_vec_fine;
        g = zeros(size(x_grid));
        ptr_mask = true(size(g));  % Reset PTR for new grid
    end
end
end

function [x_refined, z_refined] = refine_grid(x_vec, z_vec, mask, factor)
% Refines the x and z grid where mask is true
% Simple implementation: globally refine for now

x_refined = linspace(min(x_vec), max(x_vec), length(x_vec) * factor);
z_refined = linspace(min(z_vec), max(z_vec), length(z_vec) * factor);
end
