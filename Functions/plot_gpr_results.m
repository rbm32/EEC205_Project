function plot_gpr_results(g, bscan, reconstructed, x, time, image_size, dB_cutoff)
% PLOT_GPR_RESULTS - Visualizes GPR simulation results.
%
% Inputs:
%   g            - Ground truth reflectivity matrix
%   bscan        - Simulated B-scan data (time x antennas)
%   reconstructed- Reconstructed image (linear scale)
%   x            - Antenna positions
%   time         - Time vector (s)
%   image_size   - Number of image grid points (NxN)
%   dB_cutoff    - Minimum dB display value (e.g., -30)

    if nargin < 7
        dB_cutoff = -3; % Default if not provided
    end

    % Convert reconstruction to dB scale
    reconstructed_dB = 20 * log10(abs(reconstructed) / max(reconstructed(:)));
    reconstructed_dB(reconstructed_dB < dB_cutoff) = dB_cutoff;
    xmin = min(x);
    xmax = max(x);

    % Start plotting
    clf;
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % 1) Ground truth reflectivity matrix
    nexttile;
    imagesc(linspace(xmin, xmax, image_size), linspace(0, 2, image_size), g);
    colormap(gray); colorbar;
    xlabel('x (m)');
    ylabel('z (m)');
    title('Reflectivity Matrix (g)', 'FontSize', 23);
    set(gca, 'YDir', 'reverse');

    % 2) Simulated B-scan
    nexttile;
    imagesc(x, time * 1e9, bscan);
    colormap(gray); colorbar;
    xlabel('Antenna Position (m)');
    ylabel('Time (ns)');
    title('Simulated B-Scan', 'FontSize', 23);

    % 3) Reconstructed image in dB
    nexttile;
    imagesc(linspace(xmin, xmax, image_size), linspace(0, 2, image_size), reconstructed_dB);
    colormap(jet); colorbar;
    xlabel('x (m)');
    ylabel('z (m)');
    title(sprintf('Reconstructed Image (dB, Cutoff: %d)', dB_cutoff), 'FontSize', 23);
    set(gca, 'YDir', 'reverse');
end
