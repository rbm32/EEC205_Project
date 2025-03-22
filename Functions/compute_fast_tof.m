function [tof x_r] = compute_fast_tof(x_a, h, x_b, z_b, epsilon_r, c)
% compute_fast_tof - Fast approximation of time-of-flight (TOF) using Zhou et al. method.
%
% Inputs:
%   x_a        - Antenna x-position
%   h          - Antenna height above ground (m)
%   x_b        - Buried target x-position
%   z_b        - Buried target depth (m)
%   epsilon_r  - Relative permittivity of ground
%   c          - Speed of light in vacuum (m/s)
%
% Output:
%   tof        - Estimated two-way time-of-flight (seconds)

    % Critical delta and threshold from paper
    delta = z_b + h;
    threshold = delta * sqrt(epsilon_r / (epsilon_r - 1));

    % Approximate refraction point
    if abs(x_a - x_b) < threshold
        x1 = x_b + (x_b - x_a);  % Mirror point across ground
        x_r = x_b + (x1 - x_b) / sqrt(epsilon_r);
    elseif x_a > x_b + threshold
        x_r = x_b + z_b / sqrt(epsilon_r - 1);
    else
        x_r = x_b - z_b / sqrt(epsilon_r - 1);
    end

    % Travel distances
    d_air = sqrt(h^2 + (x_r - x_a)^2);
    d_soil = sqrt(z_b^2 + (x_b - x_r)^2);

    % TOF calculation (round trip)
    v = c / sqrt(epsilon_r);
    tof = 2 * (d_air / c + d_soil / v);
end
