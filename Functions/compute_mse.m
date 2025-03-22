function mse = compute_mse(g, r)
% COMPUTE_MSE - Calculates the mean squared error between two matrices
%
% Inputs:
%   g - Ground truth reflectivity matrix
%   r - Reconstructed reflectivity matrix
%
% Output:
%   mse - Mean squared error value

    % Ensure the matrices are the same size
    if ~isequal(size(g), size(r))
        error('Input matrices must be the same size.');
    end

    % Compute mean squared error
    diff = g - r;
    mse = mean(diff(:).^2);
end
