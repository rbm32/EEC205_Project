function plot_g_and_bscan(g, bscan, x, t, z, dB_cutoff, recon_title, bscan_title)
% PLOT_GPR_RESULTS_REALDATA - Visualizes GPR reconstruction and B-scan.
%
% Inputs:
%   reconstructed - Reconstructed image (linear scale)
%   bscan         - Measured B-scan data (time x antennas)
%   x             - Antenna positions
%   t             - Time vector (s)
%   dB_cutoff     - Minimum dB display value (e.g., -30)
%   recon_title   - (Optional) Custom title for reconstructed image
%   bscan_title   - (Optional) Custom title for B-scan

    if nargin < 6 || isempty(dB_cutoff)
        dB_cutoff = -3; % Default cutoff
    end
    if nargin < 7 || isempty(recon_title)
        recon_title = sprintf('Reconstructed Image (dB, Cutoff: %d)', dB_cutoff);
    end
    if nargin < 8 || isempty(bscan_title)
        bscan_title = 'B-Scan';
    end

    % Convert reconstruction to dB scale
    reconstructed_dB = 20 * log10(abs(g) / max(g(:)));
    reconstructed_dB(reconstructed_dB < dB_cutoff) = dB_cutoff;

    % Plot
    clf;
    tiledlayout('horizontal')

    

    % B-scan
    nexttile;
    imagesc(x, t * 1e9, bscan);
    colormap(gray); colorbar;
    xlabel('Antenna Position (m)');
    ylabel('Time (ns)');
    title(bscan_title, FontSize=30);
    
    % Reconstructed image
    nexttile;
    imagesc(x, z, reconstructed_dB);
    colormap(jet); colorbar;
    xlabel('x (m)');
    ylabel('z (m)');
    title(recon_title, FontSize=30);
    set(gca, 'YDir', 'reverse');
    colormap jet
end
