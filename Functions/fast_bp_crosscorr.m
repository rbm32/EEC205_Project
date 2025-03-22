function reconstructed_image = fast_bp_crossover(bscan, x, time, c, h, epsilon_r, image_size, alpha, correlation_threshold)
    % Semi-vectorized SBP for faster GPR image reconstruction

    [num_samples, num_antennas] = size(bscan);
    dt = time(2) - time(1);

    % Grid setup
    x_grid = linspace(min(x), max(x), image_size);
    z_grid = linspace(0, max(time) * c / (2 * sqrt(epsilon_r)), image_size);
    [X, Z] = meshgrid(x_grid, z_grid);
    reconstructed_image = zeros(size(X));

    % Pulse window size for extracting segments
    N = round(length(time) / 40);  % ~10% of time axis
    seg_len = 2 * N + 1;

    % Precompute all travel times
    tof_lookup = zeros(image_size, image_size, num_antennas);
    for m = 1:num_antennas
        for i = 1:image_size
            for j = 1:image_size
                x0 = X(i, j);
                z0 = Z(i, j);
                tof_lookup(i, j, m) = compute_fast_tof(x(m), h, x0, z0, epsilon_r, c);
            end
        end
    end

    % Main loop over image pixels
    for i = 1:image_size
        for j = 1:image_size
            x0 = X(i, j);
            z0 = Z(i, j);

            % Center antenna index
            [~, idx_center] = min(abs(x - x0));
            t_c = tof_lookup(i, j, idx_center);
            [~, t_idx_c] = min(abs(time - t_c));

            % Central segment
            C_center = extract_segment(bscan(:, idx_center), t_idx_c, N);
            u_vals = abs(bscan(t_idx_c, idx_center));
            segments = {C_center};

            % Search left
            for k = idx_center-1:-1:1
                t = tof_lookup(i, j, k);
                [~, t_idx] = min(abs(time - t));
                Ck = extract_segment(bscan(:, k), t_idx, N);
                r = corrcoef(C_center, Ck);
                if size(r,1) < 2 || r(1,2) < correlation_threshold
                    break;
                end
                u_vals(end+1) = abs(bscan(t_idx, k));
                segments{end+1} = Ck;
            end

            % Search right
            for k = idx_center+1:num_antennas
                t = tof_lookup(i, j, k);
                [~, t_idx] = min(abs(time - t));
                Ck = extract_segment(bscan(:, k), t_idx, N);
                r = corrcoef(C_center, Ck);
                if size(r,1) < 2 || r(1,2) < correlation_threshold
                    break;
                end
                u_vals(end+1) = abs(bscan(t_idx, k));
                segments{end+1} = Ck;
            end

            % Cross-correlation energy sum
            M = length(u_vals);
            if M > 1
                U = u_vals(:);
                E = sum(triu(U * U', 1), 'all');
            else
                E = 0;
            end

            % Depth compensation
            compensation = exp(alpha * z0);
            reconstructed_image(i, j) = compensation * E;
        end
    end

    % Normalize output
    reconstructed_image = reconstructed_image / max(reconstructed_image(:));
end

function segment = extract_segment(signal, center_idx, N)
    % Extract segment of length 2N+1 around center_idx
    L = length(signal);
    idx_start = max(1, center_idx - N);
    idx_end = min(L, center_idx + N);
    segment = signal(idx_start:idx_end);
    if length(segment) < 2*N+1
        segment(end+1:2*N+1) = 0;
    end
end