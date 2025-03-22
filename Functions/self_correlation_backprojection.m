function reconstructed_image = self_correlation_backprojection(bscan, x, time, c, h, epsilon_r, image_size, alpha, correlation_threshold)
    % Implements the SBP algorithm from Zhang et al. (2015)
    %
    % Inputs:
    %   bscan                - B-scan matrix (time x antennas)
    %   x                   - Antenna positions
    %   time                - Time vector (s)
    %   c                   - Speed of light in vacuum (m/s)
    %   h                   - Antenna height above ground (m)
    %   epsilon_r           - Relative permittivity of ground
    %   image_size          - Number of pixels per axis in the output image
    %   alpha               - Depth compensation coefficient
    %   correlation_threshold - Threshold ψ for adaptive echo inclusion (typically 0.85–0.95)
    %
    % Output:
    %   reconstructed_image - SBP-reconstructed reflectivity image

    [num_samples, num_antennas] = size(bscan);
    dt = time(2) - time(1);

    % Set up image grid
    x_grid = linspace(min(x), max(x), image_size);
    z_grid = linspace(0, max(time) * c / (2 * sqrt(epsilon_r)), image_size);
    [X, Z] = meshgrid(x_grid, z_grid);
    reconstructed_image = zeros(size(X));

    N = round(length(time) / 40); % Length of echo window (10% of B-scan depth range, can be tuned)

    for i = 1:image_size
        for j = 1:image_size
            x0 = X(i, j);
            z0 = Z(i, j);

            % Find the antenna closest to x0
            [~, idx_center] = min(abs(x - x0));
            u_vals = [];
            C_seqs = {};

            % Compute central response
            [tof_c, ~] = compute_tof(x(idx_center), h, x0, z0, epsilon_r, c, false);
            [~, t_idx_c] = min(abs(time - tof_c));
            u_center = abs(bscan(t_idx_c, idx_center));
            C_center = extract_segment(bscan(:, idx_center), t_idx_c, N);
            u_vals(end+1) = u_center;
            C_seqs{end+1} = C_center;

            % Search left
            for k = idx_center-1:-1:1
                [tof, ~] = compute_tof(x(k), h, x0, z0, epsilon_r, c, false);
                [~, t_idx] = min(abs(time - tof));
                C_k = extract_segment(bscan(:, k), t_idx, N);
                r = corrcoef(C_center, C_k);
                if size(r,1) < 2 || r(1,2) < correlation_threshold
                    break;
                end
                u_vals(end+1) = abs(bscan(t_idx, k));
                C_seqs{end+1} = C_k;
            end

            % Search right
            for k = idx_center+1:1:num_antennas
                [tof, ~] = compute_tof(x(k), h, x0, z0, epsilon_r, c, false);
                [~, t_idx] = min(abs(time - tof));
                C_k = extract_segment(bscan(:, k), t_idx, N);
                r = corrcoef(C_center, C_k);
                if size(r,1) < 2 || r(1,2) < correlation_threshold
                    break;
                end
                u_vals(end+1) = abs(bscan(t_idx, k));
                C_seqs{end+1} = C_k;
            end

            % Cross-correlation sum
            energy = 0;
            for m = 1:length(u_vals)-1
                for n = m+1:length(u_vals)
                    energy = energy + u_vals(m) * u_vals(n);
                end
            end

            % Depth compensation
            compensation = exp(alpha * z0);
            reconstructed_image(i, j) = compensation * energy;
        end
    end

    % Normalize
    reconstructed_image = reconstructed_image / max(reconstructed_image(:));
end

function segment = extract_segment(signal, center_idx, N)
    % Extracts a segment of length 2N centered around center_idx
    L = length(signal);
    idx_start = max(1, center_idx - N);
    idx_end = min(L, center_idx + N);
    segment = signal(idx_start:idx_end);
    % Pad if needed
    if length(segment) < 2*N+1
        segment(end+1:2*N+1) = 0;
    end
end
