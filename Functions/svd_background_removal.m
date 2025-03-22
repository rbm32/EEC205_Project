function D_BKGR = svd_background_removal(D, L)
%SVD_BACKGROUND_REMOVAL Removes background using SVD truncation.
%   D_BKGR = SVD_BACKGROUND_REMOVAL(D, L) performs background removal on
%   matrix D by zeroing out the first L principal components (largest
%   singular values) and reconstructing the matrix.
%
%   Inputs:
%       D - Input data matrix (e.g., B-scan data)
%       L - Number of leading singular values to zero (background components)
%
%   Output:
%       D_BKGR - Background-removed matrix

    % Compute reduced SVD
    [U, S, V] = svd(D, 'econ');

    % Zero out the first L singular values
    o = diag(S);
    S_filtered = diag([zeros(L,1); o(L+1:end)]);

    % Pad S_filtered if needed (in case it's smaller due to econ mode)
    S_padded = zeros(size(S));
    S_padded(1:size(S_filtered,1), 1:size(S_filtered,2)) = S_filtered;

    % Reconstruct matrix with background removed
    D_BKGR = U * S_padded * V';
end