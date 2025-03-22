function [bscan_ds, t_ds, x_ds] = downsample_bscan(bscan, t, x, time_factor, x_factor)
    % DOWNSAMPLE_BSCAN Downsamples the B-scan in time and x (antenna) directions
    %
    % Inputs:
    %   bscan       - Original B-scan matrix (time x antenna positions)
    %   t           - Original time vector (length equal to number of rows in bscan)
    %   x           - Original x (antenna position) vector (length equal to number of columns in bscan)
    %   time_factor - Downsampling factor for time (positive integer)
    %   x_factor    - Downsampling factor for x (positive integer)
    %
    % Outputs:
    %   bscan_ds - Downsampled B-scan
    %   t_ds     - Downsampled time vector
    %   x_ds     - Downsampled x vector (antenna positions)

    % Input validation
    if any([time_factor, x_factor] <= 0) || any(mod([time_factor, x_factor],1) ~= 0)
        error('Both downsampling factors must be positive integers.');
    end

    % Perform downsampling
    bscan_ds = bscan(1:time_factor:end, 1:x_factor:end);
    t_ds = t(1:time_factor:end);
    x_ds = x(1:x_factor:end);
end
