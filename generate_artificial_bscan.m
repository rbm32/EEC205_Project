%% --- Generate Artificial B-scan Function ---
function [bscan, time, x] = generate_artificial_bscan(image_size, num_time_samples, x, g, pulse, c, t_max, a)
    % Generates an artificial B-scan of soil with a given reflectivity matrix g
    %
    % Inputs:
    %   image_size       - Size of the reflectivity matrix
    %   num_time_samples - Number of time samples
    %   x                - Linear antenna positions
    %   g                - Reflectivity matrix (image_size x image_size)
    %   pulse            - Custom pulse waveform
    %   c                - Wave propagation speed
    %   t_max            - Maximum time value
    %   a                - Noise matrix (optional, size must match image_size)
    %
    
    % Outputs:
    %   bscan - Simulated B-scan data (time samples x antenna positions)

    % Set parameters
    x_grid = linspace(-1, 1, image_size); % Cross-range axis
    z = linspace(0, 2, image_size); % Range axis
    [X, Z] = meshgrid(x_grid, z);

    % Define time samples
    time = linspace(0, t_max, num_time_samples); % Time samples

    % Initialize B-scan matrix for time-domain data
    num_positions = length(x);
    xs = mean(x(2:end) - x(1:end-1));
    bscan = zeros(num_time_samples, num_positions);

    % Normalize pulse to prevent amplitude scaling issues
    pulse = pulse / max(abs(pulse));

    % Simulate time-domain data using echo responses
    for m = 1:num_positions
        x_pos = x(m); % Current antenna position
        for i = 1:image_size
            for j = 1:image_size
                if g(i, j) > 0
                    % Edge extension handling
                    if j == 1
                        kidx = 0:-1:-500;
                    elseif j == size(g,2)
                        kidx = 0:500;
                    else
                        kidx = 0;
                    end

                    for k = kidx
                        r = sqrt((X(i, j) - x_pos + xs .* k).^2 + Z(i, j).^2);
                        t = 2 * r / c; % Two-way travel time

                        % Find closest sample in time vector
                        [~, idx] = min(abs(time - t));

                        if idx > 0 && idx <= num_time_samples
                            % Add pulse response at the computed time index
                            pulse_len = length(pulse);
                            end_idx = min(idx + pulse_len - 1, num_time_samples);
                            bscan(idx:end_idx, m) = bscan(idx:end_idx, m) + ...
                                g(i, j) * pulse(1:(end_idx - idx + 1)).';
                        end
                    end
                end
            end
        end
    end

    % Add noise if provided
    if nargin == 8 && ~isempty(a)
        % Ensure noise matrix is the correct size
        if ~isequal(size(a), [num_time_samples, image_size])
            error('Noise matrix "a" must be the same size as the b-scan.');
        end

        % Map noise directly to the B-scan size
        bscan = bscan + a;
    end
end
