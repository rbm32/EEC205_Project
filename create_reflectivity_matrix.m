function g = create_reflectivity_matrix(image_size, targets, ground_z, ground_width, ground_reflectivity)
    % Generates the reflectivity matrix g from given target locations and a ground surface
    % targets: a Nx4 matrix where each row is [x, z, radius, reflectivity]
    % ground_z: normalized z value for the ground surface (e.g., 1.0 for a flat ground at z = 1.0)
    % ground_width: width of the ground surface (range in z axis) for ground reflectivity
    % ground_reflectivity: reflectivity value for the ground surface
    
    % Create a grid of positions
    x_grid = linspace(-1, 1, image_size); % Cross-range axis
    z = linspace(0, 2, image_size); % Range axis
    [X, Z] = meshgrid(x_grid, z);
    
    
    % Initialize reflectivity matrix g
    g = zeros(image_size);
    
    % Add circular targets to g based on input targets
    for i = 1:size(targets, 1)
        dist = sqrt((X - targets(i, 1)).^2 + (Z - targets(i, 2)).^2);
        g(dist <= targets(i, 3)) = targets(i, 4); % Set reflectivity based on target
    end
    
    % Add ground surface to g
    % Define the z-range where the ground is located
    ground_row_start = find(z >= ground_z - ground_width / 2, 1); % Start of ground in z-axis
    ground_row_end = find(z <= ground_z + ground_width / 2, 1, 'last'); % End of ground in z-axis
    
    % Apply ground reflectivity to the specified range in the z-direction
    if ~isempty(ground_row_start) && ~isempty(ground_row_end)
        g(ground_row_start:ground_row_end, :) = ground_reflectivity; % Set the ground reflectivity
    end
end