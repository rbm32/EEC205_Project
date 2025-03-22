function reconstructed_image = BScanBackprojection(bscan, x, time, c, h, epsilon_r, image_size)
    % Refraction-aware filtered back-projection using compute_tof.
    %
    % Inputs:
    %   bscan       - B-scan matrix (time samples x antenna positions)
    %   x           - Antenna positions (1 x M)
    %   time        - Time vector (1 x T)
    %   c           - Speed of light in vacuum (m/s)
    %   h           - Antenna height above ground (m)
    %   epsilon_r   - Relative permittivity of ground
    %   image_size  - Reconstruction grid size (NxN)
    %
    % Output:
    %   reconstructed_image - Reconstructed reflectivity image (image_size x image_size)

    % Create reconstruction grid
    x_grid = linspace(min(x), max(x), image_size);
    z_max = max(time) * c / (2 * sqrt(epsilon_r));
    z_grid = linspace(h, z_max + h, image_size);  % Depth grid shifted up by antenna height
    [X, Z] = meshgrid(x_grid, z_grid);

    % Initialize image
    reconstructed_image = zeros(image_size, image_size);

    % Loop over all antenna positions and image pixels
    for m = 1:length(x)
        x_pos = x(m);
        for i = 1:image_size
            for j = 1:image_size
                x_p = X(i, j);
                z_p = Z(i, j);

                % Get two-way travel time using refraction-aware method
                [tof, ~] = compute_tof(x_pos, h, x_p, z_p, epsilon_r, c, false);

                % Map ToF to nearest time index
                [~, idx] = min(abs(time - tof));

                % Accumulate if index is valid
                if idx > 0 && idx <= size(bscan, 1)
                    reconstructed_image(i, j) = reconstructed_image(i, j) + abs(bscan(idx, m));
                end
            end
        end
    end

    % Normalize output
    reconstructed_image = reconstructed_image / max(reconstructed_image(:));
end
