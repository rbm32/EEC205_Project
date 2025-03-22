function g = reconstruct_g_bayesian(H, bscan, num_scans, num_time_samples, lambda, L)
    % Reconstruct g using Bayesian MAP with quadratic prior
    %
    % Inputs:
    %   H           - Forward model matrix (T*M x N)
    %   bscan       - B-scan matrix (T x M)
    %   image_size  - Scalar image dimension (for reshaping to square)
    %   lambda      - Regularization strength
    %   L           - Regularization operator (defaults to identity if not provided)
    %
    % Output:
    %   g           - Reconstructed reflectivity matrix [image_size x image_size]

    % Flatten B-scan to vector
    b = bscan(:);
    N = size(H, 2);  % Number of image pixels

    % Default L to identity if not specified
    if nargin < 6 || isempty(L)
        L = speye(N);
    end

    % Compute MAP estimate
    A = H' * H + lambda * (L' * L);
    rhs = H' * b;

    % Solve (regularized normal equation)
    g_vec = A \ rhs;

    % Reshape to 2D image
    g = reshape(g_vec, num_time_samples, num_scans);
end