addpath(genpath('utils_Ryan\'));
addpath(genpath('Data\'));
clearvars; close all; clc;

pulsed_filename = "Data\ryanPulse_v2.csv";

S = Scan_V2();
filename = pulsed_filename;
S.read_scope_data(filename)
S.BKGR_PCA(2)
subplot(1,2,1)
S.plot_Ryan("Raw Data", "Raw Data")
subplot(1,2,2)
S.plot_Ryan("BKGR", "BKGR Removed");
%%

bscan = S.D_BKGR;
time = S.t;
x = S.x;
x = linspace(0, 1.5, size(bscan, 2));
c = 3e8;

% Perform Filtered Backprojection to estimate the targets
% Backprojection
reconstructed_image = filtered_backprojection(bscan, x, time, c);

% Convert to dB


%% Plot reconstructed image


reconstructed_image_dB = 20 * log10(abs(reconstructed_image) + eps); % Avoid log(0)

% Apply cutoff
cutoff_dB = 5; % Cutoff 20 dB below the max value
max_dB = max(reconstructed_image_dB(:));
reconstructed_image_dB(reconstructed_image_dB < (max_dB - cutoff_dB)) = NaN;


tiledlayout('horizontal')
nexttile;
imagesc(linspace(min(x), max(x), 100), linspace(0, max(time) * 0.3 / 2, 100), reconstructed_image_dB);
xlabel('Antenna Position (m)');
ylabel('Depth (m)');
title('Reconstructed Image (dB)');
colormap('jet'); % Set reconstructed image colormap to jet
colorbar;
axis tight;


tl = nexttile;
imagesc(x, time, (bscan));
xlabel('Antenna Position (m)');
ylabel('Time (ns)');
title('Artificial B-scan (Custom Pulse)');

colorbar;
axis tight;
colormap(tl, 'gray'); % Set B-scan colormap to gray

