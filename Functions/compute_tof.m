function [tof, x_r] = compute_toF(x_a, h, x_b, z, epsilon_r, c, do_plot)
    if nargin < 7
        do_plot = false;
    end

    % Special case: vertical path (no refraction)
    if abs(x_a - x_b) < 1e-10
        x_r = x_a;  % No horizontal offset
        d_air = 2 * h;
        d_ground = 2 * z;
        v_ground = c / sqrt(epsilon_r);

        t_air = d_air / c;
        t_ground = d_ground / v_ground;
        tof = t_air + t_ground;

        if do_plot
            plot_ray(x_a, x_r, x_b, h, z);
        end
        return;
    end

    % Otherwise, use Snell’s law to find x_r
    x_r_guess = (x_a + x_b) / 2;
    refraction_eq = @(x_r) snell_constraint(x_r, x_a, h, x_b, z, epsilon_r);
    options = optimset('Display', 'off');
    x_r = fsolve(refraction_eq, x_r_guess, options);

    % Compute distances
    d_air = 2 * sqrt(h^2 + (x_r - x_a)^2);
    d_ground = 2 * sqrt(z^2 + (x_b - x_r)^2);
    v_ground = c / sqrt(epsilon_r);

    % Compute ToF
    t_air = d_air / c;
    t_ground = d_ground / v_ground;
    tof = t_air + t_ground;

    % Optional plotting
    if do_plot
        plot_ray(x_a, x_r, x_b, h, z);
    end
end

function F = snell_constraint(x_r, x_a, h, x_b, z, epsilon_r)
    delta_x_air = x_r - x_a;
    delta_x_ground = x_b - x_r;

    if abs(delta_x_air) < 1e-12  % Vertical air path
        theta_i = 0;
    else
        theta_i = atan(abs(delta_x_air) / h);
    end

    sin_theta_t = sin(theta_i) / sqrt(epsilon_r);
    sin_theta_t = min(max(sin_theta_t, -1), 1);  % Clamp to valid range
    theta_t = asin(sin_theta_t);

    if abs(tan(theta_t)) < 1e-12
        F = 0;  % Avoid division by zero
    else
        F = abs(delta_x_ground) / tan(theta_t) - z;
    end
end

function plot_ray(x_a, x_r, x_b, h, z)
    clf; hold on; grid on;
    xlim([min([x_a, x_b]) - 1, max([x_a, x_b]) + 1]);
    ylim([-z - 1, h + 1]);
    line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.2); % Ground line
    plot(x_a, h, 'bo', 'MarkerFaceColor', 'b'); % Antenna
    plot(x_r, 0, 'ko', 'MarkerFaceColor', 'k'); % Refraction point
    plot(x_b, -z, 'ro', 'MarkerFaceColor', 'r'); % Target
    line([x_a x_r], [h 0], 'Color', 'b', 'LineWidth', 1.5); % Air path
    line([x_r x_b], [0 -z], 'Color', 'r', 'LineWidth', 1.5); % Ground path
    legend('Ground', 'Antenna', 'Refraction Pt.', 'Target', 'Air Path', 'Ground Path');
    xlabel('x (m)'); ylabel('Depth (m)'); set(gca, 'YDir', 'reverse');
    title('Radar Pulse Path with Refraction');
end
