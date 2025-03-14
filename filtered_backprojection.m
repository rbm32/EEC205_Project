function reconstructed_image = filtered_backprojection(bscan, x, time, c)
    % Perform filtered back-projection on the B-scan data.
    % 
    % Inputs:
    %   bscan - B-scan matrix (time samples x antenna positions)
    %   x - Antenna positions
    %   time - Time vector
    %   c - Speed of wave (e.g., 0.3 m/ns for EM waves in soil)
    %
    % Outputs:
    %   reconstructed_image - Reconstructed reflectivity image
    
    % Define reconstruction grid
    image_size = 100;
    x_grid = linspace(min(x), max(x), image_size);
    z_grid = linspace(0, max(time) * c / 2, image_size);
    [X, Z] = meshgrid(x_grid, z_grid);
    
    % Initialize reconstructed image
    reconstructed_image = zeros(size(X));
    
    % Back-projection loop
    for m = 1:length(x)
        x_pos = x(m);
        for i = 1:image_size
            for j = 1:image_size
                % Compute range from antenna to point in image
                r = sqrt((X(i, j) - x_pos)^2 + Z(i, j)^2);
                t = 2 * r / c; % Two-way travel time
                
                % Find the closest time index
                [~, idx] = min(abs(time - t));
                
                if idx > 0 && idx <= size(bscan, 1)
                    % Back-project signal to reconstructed image
                    reconstructed_image(i, j) = reconstructed_image(i, j) + abs(bscan(idx, m));
                end
            end
        end
    end
end