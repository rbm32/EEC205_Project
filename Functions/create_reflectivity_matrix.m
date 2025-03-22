function g = create_reflectivity_matrix(x, z, targets)
    % Generates the reflectivity matrix g from given target locations
    % Sets reflectivity only at the nearest pixel for each target
    %
    % Inputs:
    %   x       - Vector of horizontal (cross-range) coordinates
    %   z       - Vector of vertical (depth/range) coordinates
    %   targets - Nx3 matrix, each row: [x_center, z_center, reflectivity]
    %
    % Output:
    %   g       - Reflectivity matrix of size length(z) x length(x)

    % Initialize reflectivity matrix
    g = zeros(length(z), length(x));

    % Loop through targets and assign to nearest grid point
    for i = 1:size(targets, 1)
        xc = targets(i, 1);
        zc = targets(i, 2);
        val = targets(i, 3);

        % Find nearest grid indices
        [~, xi] = min(abs(x - xc));
        [~, zi] = min(abs(z - zc));

        % Assign reflectivity
        g(zi, xi) = val;
    end
end
