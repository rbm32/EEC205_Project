function [bscan_weighted, weights] = compensate_attenuation(bscan, t, c, epsilon_r, alpha)
    % COMPENSATE_ATTENUATION Pre-weights a B-scan to undo exponential attenuation
    %
    % Inputs:
    %   bscan     - Raw B-scan data [T x M] (time x antenna)
    %   t         - Time vector [T x 1]
    %   c         - Speed of light in vacuum
    %   epsilon_r - Relative permittivity of soil
    %   alpha     - Attenuation factor (Np/m)
    %
    % Outputs:
    %   bscan_weighted - Pre-weighted B-scan with attenuation compensation
    %   weights        - Vector of weights applied per time sample

    if size(t,1) == 1
        t = t(:);  % ensure column vector
    end

    v = c / sqrt(epsilon_r);   % wave speed in medium
    z = v * t / 2;             % depth approximation for each time sample (round-trip)
    weights = exp(alpha * z);  % inverse of attenuation exp(-alpha*z)

    % Apply weights to each column (antenna)
    bscan_weighted = bscan .* weights;
end
