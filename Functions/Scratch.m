% --- Parameters ---
c = 3e8;                     % Speed of light
epsilon_r = 9;              % Relative permittivity
h = 1.0;                    % Antenna height
x_ant = linspace(-1, 1, 64); % Antenna positions
x_vec = linspace(-1, 1, 64); % Image lateral axis
z_vec = linspace(0.1, 1.5, 64); % Image depth axis

% --- Select a pixel location to analyze ---
ix = 32; % Mid x
iz = 32; % Mid z
x0 = x_vec(ix);
z0 = z_vec(iz);

% --- Initialize error array ---
tof_error = zeros(size(x_ant));
tof_fast_all = zeros(size(x_ant));
tof_true_all = zeros(size(x_ant));

for m = 1:length(x_ant)
    xk = x_ant(m);

    % --- Fast Approximation ---
    delta = z0 + h;
    threshold = delta * sqrt(epsilon_r / (epsilon_r - 1));
    if abs(xk - x0) < threshold
        x1 = x0 + (x0 - xk);
        xr_fast = x0 + (x1 - x0) / sqrt(epsilon_r);
    elseif xk > x0 + threshold
        xr_fast = x0 + z0 / sqrt(epsilon_r - 1);
    else
        xr_fast = x0 - z0 / sqrt(epsilon_r - 1);
    end

    d_air_fast = sqrt(h^2 + (xr_fast - xk)^2);
    d_soil_fast = sqrt(z0^2 + (x0 - xr_fast)^2);
    tof_fast = 2 * (d_air_fast / c + d_soil_fast / (c / sqrt(epsilon_r)));

    % --- True TOF ---
    [tof_true, ~] = compute_tof(xk, h, x0, z0, epsilon_r, c, false);

    % --- Store results ---
    tof_fast_all(m) = tof_fast;
    tof_true_all(m) = tof_true;
    tof_error(m) = abs(tof_fast - tof_true);
end

% --- Plot TOF Error ---
figure;
plot(x_ant, tof_error * 1e9, 'r', 'LineWidth', 2); % Convert error to ns
xlabel('Antenna Position x (m)');
ylabel('TOF Error (ns)');
title(sprintf('TOF Approximation Error at Point (x=%.2f, z=%.2f)', x0, z0));
grid on;
