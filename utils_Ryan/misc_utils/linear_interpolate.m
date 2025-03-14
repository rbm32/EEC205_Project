function M = linear_interpolate(M, dim, interp_ratio)
    % A "wrapper" function to interpolate a matrix along
    % a specific dimension. Linear interpolate.
    % no input checking. we go fast here.
    
    % Change M to have dim on dim1
    M           = permute_matrix(M, dim, size(M));
    size_M      = size(M);
    M           = interp1(...
                        (0:(size_M(1)-1)), ...
                        M,...
                        (0:1/interp_ratio:(size_M(1)-1)),...
                        'linear'...
                    );
    M           = permute_matrix(M, dim, size_M);
end

function M = permute_matrix(M, dim, sz)
    % permutes the matrix such that `dim` is along dim1
    % create the permutation vector
    perm_vec        = 1:numel(sz);
    perm_vec(dim)   = 1;
    perm_vec(1)     = dim;
    
    % permute the matrix
    M = permute(M, perm_vec);
end