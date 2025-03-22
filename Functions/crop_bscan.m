function [bscan_cropped, x_cropped, t_cropped] = crop_bscan(raw, x, t, x_range, t_range)
    % CROP_BSCAN crops a B-scan in time and space.
    %
    % Inputs:
    %   raw     - B-scan matrix (size: length(t) x length(x))
    %   x       - Vector of spatial (antenna) positions
    %   t       - Vector of time samples (in seconds)
    %   x_range - 2-element vector [x_min, x_max] for spatial crop
    %   t_range - 2-element vector [t_min, t_max] for time crop
    %
    % Outputs:
    %   bscan_cropped - Cropped B-scan matrix
    %   x_cropped     - Cropped x vector
    %   t_cropped     - Cropped t vector

    % Find indices in range
    x_idx = find(x >= x_range(1) & x <= x_range(2));
    t_idx = find(t >= t_range(1) & t <= t_range(2));

    % Crop B-scan and axes
    bscan_cropped = raw(t_idx, x_idx);
    x_cropped = x(x_idx);
    t_cropped = t(t_idx);
end
