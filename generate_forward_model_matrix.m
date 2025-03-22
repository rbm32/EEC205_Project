function [H, time, x_vec, X, Z] = generate_forward_model_matrix(x_vec, z_vec, num_time_samples, pulse, c, t_max, epsilon_r, h, alpha)
    % Optimized forward model matrix generator (parallel + vectorization)
    %
    % Outputs:
    %   H     - Forward model matrix (num_time_samples*num_antennas × numel(g))
    %   time  - Time vector
    %   x_vec - Antenna positions
    %   X, Z  - Meshgrid of x/z (for reference)

    % Generate spatial grid
    [X, Z] = meshgrid(x_vec, z_vec);
    X_flat = X(:);
    Z_flat = Z(:);

    % Time axis
    time = linspace(0, t_max, num_time_samples);
    pulse_len = length(pulse);

    % Sizes
    num_antennas = length(x_vec);
    num_pixels = numel(X);

    % Preallocate output per antenna
    H_cells = cell(num_antennas, 1);

    % Parallel loop over antennas
    parfor m = 1:num_antennas
        x_pos = x_vec(m);
        H_local = sparse(num_time_samples, num_pixels);  % Each antenna has a chunk of rows

        for p = 1:num_pixels
            x_b = X_flat(p);
            z_b = Z_flat(p);
            % Compute time-of-flight and refraction point
            [tof, x_r] = compute_tof(x_pos, h, x_b, z_b, epsilon_r, c, false);

            % Path distances
            d_air = 2 * sqrt(h^2 + (x_r - x_pos)^2);
            d_soil = 2 * sqrt(z_b^2 + (x_b - x_r)^2);
            R_total = d_air + d_soil;

            % Attenuation
            att_R = 1 / (R_total^2);
            att_alpha = exp(-alpha * d_soil);
            att = att_R * att_alpha;

            % Time index in B-scan
            [~, idx] = min(abs(time - tof));
            if idx > 0 && idx <= num_time_samples
                row_end = min(idx + pulse_len - 1, num_time_samples);
                pulse_crop = pulse(1:(row_end - idx + 1));

                % Fill in sparse matrix block
                H_local(idx:row_end, p) = att * pulse_crop(:);
            end
        end

        H_cells{m} = H_local;
    end

    % Assemble full H matrix from antenna-wise blocks
    H = spalloc(num_time_samples * num_antennas, num_pixels, num_antennas * num_pixels * pulse_len);
    for m = 1:num_antennas
        row_range = (m - 1) * num_time_samples + (1:num_time_samples);
        H(row_range, :) = H_cells{m};
    end
end
