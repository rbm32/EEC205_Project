function M = fft_interpolate(M, dim, interp_ratio)
    % no input checking. we go fast here.
    
    % Change M to have dim on dim1
    M           = permute_matrix(M, dim, size(M));
    size_M      = size(M);
    
    M_f         = fftshift(fft(M,[],1),1);
    
    new_size    = size_M(1) * interp_ratio - (interp_ratio - 1);
    extra_size  = new_size - size_M(1);
    T_size      = size_M;
    T_size(1)   = floor(extra_size/2);
    B_size      = size_M;
    B_size(1)   = floor(extra_size/2)+ mod(extra_size,2) ;
    zeropad_T   = zeros(T_size);
    zeropad_B   = zeros(B_size);
    
    M_f         = cat(1, zeropad_T, M_f, zeropad_B);
    
    M           = real(ifft(ifftshift(M_f,1),[], 1));
    size_M(1)   = new_size;

    M           = permute_matrix(M, dim, size_M) .* interp_ratio;
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