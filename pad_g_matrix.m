%% --- Function to Pad the Reflectivity Matrix in y-direction ---
function g_padded = pad_g_matrix(g)
    % Pads the g matrix in the y-direction with the first and last columns to simulate infinite extension
    g_padded = g; % Start with the original g matrix
    n = round(size(g, 2));
    
    % Pad the matrix with the first and last column
    g_padded = [repmat(g(:, 1), 1, n), g, repmat(g(:, end), 1, n)]; % Padding with n columns on both sides
    
    % Ensure the padded size is consistent with the original image_size
    g_padded = g_padded(:, 1:end); % Trim excess columns if necessary
end