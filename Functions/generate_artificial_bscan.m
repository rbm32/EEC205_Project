function [bscan, time, x] = generate_artificial_bscan(image_size, num_time_samples, x, g, pulse, c, t_max, epsilon_r, h, alpha, a)
    % Generates an artificial B-scan of soil with a given reflectivity matrix g using compute_tof
    %
    % Inputs:
    %   image_size       - Size of the reflectivity matrix
    %   num_time_samples - Number of time samples
    %   x                - Linear antenna positions
    %   g                - Reflectivity matrix (image_size x image_size)
    %   pulse            - Custom pulse waveform
    %   c                - Speed of light in vacuum
    %   t_max            - Maximum time value
    %   epsilon_r        - Relative permittivity of the ground
    %   h                - Height of the antenna above ground (m)
    %   alpha            - Soil attenuation factor (Np/m)
    %   a                - (Optional) Noise matrix
    %
    % Outputs:
    %   bscan - Simulated B-scan data (time samples x antenna positions)

    % Set spatial grid
    x_grid = linspace(-1, 1, image_size);  % Cross-range
    z = linspace(0, 2, image_size);        % Range
    [X, Z] = meshgrid(x_grid, z);

    % Time axis
    time = linspace(0, t_max, num_time_samples);

    % Output B-scan
    num_positions = length(x);
    xs = mean(x(2:end) - x(1:end-1)); % Antenna spacing
    bscan = zeros(num_time_samples, num_positions);

    pulse_len = length(pulse);

    % Loop over antenna positions
    for m = 1:num_positions
        x_pos = x(m); % Antenna position
        for i = 1:image_size
            for j = 1:image_size
                if g(i, j) > 0
                    x_b = X(i, j);
                    z_b = Z(i, j);

                    % Use compute_tof to get time of flight and refraction point
                    [tof, x_r] = compute_tof(x_pos, h, x_b, z_b, epsilon_r, c, false);

                    % Distance segments
                    d_air = 2 * sqrt(h^2 + (x_r - x_pos)^2);
                    d_soil = 2 * sqrt(z_b^2 + (x_b - x_r)^2);
                    R_total = d_air + d_soil;

                    % Compute attenuation factors
                    attenuation_R = 1 / (R_total^2);              % 1/R² spreading loss
                    attenuation_alpha = exp(-alpha * d_soil);     % Soil absorption (not air)

                    attenuation = attenuation_R * attenuation_alpha;

                    % Find index in time axis
                    [~, idx] = min(abs(time - tof));

                    if idx > 0 && idx <= num_time_samples
                        end_idx = min(idx + pulse_len - 1, num_time_samples);
                        pulse_crop = pulse(1:(end_idx - idx + 1));
                        bscan(idx:end_idx, m) = bscan(idx:end_idx, m) + ...
                            g(i, j) * attenuation * pulse_crop.';
                    end
                end
            end
        end
    end

    % Add noise if provided
    if nargin == 11 && ~isempty(a)
        if ~isequal(size(a), size(bscan))
            error('Noise matrix "a" must match the size of the B-scan.');
        end
        bscan = bscan + a;
    end
end
