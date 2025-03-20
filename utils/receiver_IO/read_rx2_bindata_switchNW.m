function D4 = read_rx2_bindata_switchNW(fn, wfm_len, C)
    % this function will take as an input a filename and return the data
    % -generated data - interleaved as [0101][tx_id][rx_id][data1][data2] in 32bit
    % NOTE: ONLY CURRENTLY WORKS WITH RECEIVERS WHICH FLIP-ORDER
    %   I.E., sampling_reversal is selected on C.setup_order(...)
    % INPUTS _____
    % fn        : filename as an absolute path
    % wfm_len   : number of variables to be read per waveform
    % C         : calibration object
    % OUTPUTS____
    % D4: an N x M x N_tx x N_rx matrix of singles.

    %   read in the data
    assert(isfile(fn), "File " + fn + " does not exist. Make sure that the filename points to the exact path.");
    
    disp("Reading in file: " + fn);
    fp       = fopen(fn,'rb');                    % open the data file
    data_raw = uint32(fread(fp,'uint32'));  % read the entire file as uint32
    fclose(fp);                             % close the data file

    % extract upper & lower 
    % [reserved (4bit)]
    % [tx_id (3bit)]
    % [rx_id (3bit)]
    % [data_upper (11bit)]
    % [data_lower (11bit)]

    data_lower = uint16(bitand(data_raw,               uint32(2^11-1),'uint32'));
    data_upper = uint16(bitand(bitshift(data_raw,-11), uint32(2^11-1),'uint32'));
    rx_id      = uint8( bitand(bitshift(data_raw,-22), uint32(2^3 -1),'uint32'));
    tx_id      = uint8( bitand(bitshift(data_raw,-25), uint32(2^3 -1),'uint32'));
    rx_id      = rx_id(1:wfm_len/2:end) + 1;
    tx_id      = tx_id(1:wfm_len/2:end) + 1;

    D = int16(reshape([data_upper'; data_lower'], wfm_len, []));
    
    N_wfms      = size(D,2);
    N_frames    = N_wfms/C.N_tx/C.N_rx;
    D4          = zeros(wfm_len, N_frames, C.N_tx, C.N_rx, 'single');
    order_inds  = ones(C.N_tx, C.N_rx);
    
    % this is the important part of the data reshaping!
    % in this loop, we will populate the 4-dimensional sample matrix
    % with the correct waveforms.
    % this is coded semi-robustly, in that out-of-order is OK.
    %
    % order_inds(N_tx,N_rx) is a 2D matrix whose elements refer to the
    % next-to-write index of the (i-th TX, j-th RX) pair
    % order_inds starts at one, and the (i,j)th element is incremented by
    % one after the (i,j)th TXRX pair is written to.
    % tx_id and rx_id are used as indices, but incrememnted by one because
    % tx_id spans from 0...N_tx-1 and rx_id spans from 0...N_rx-1
    %
    % D4 (4D sample/data matrix) has dimensions of
    % (Time, SlowTime/Pulse#, Tx_ID, Rx_ID)
    %
    % to gain an understanding of how this resorting works, try the
    % following lines of code 
    % (outside the loop) 
    % figure(); hold on; caxis([1 N_frames]);
    % (inside the loop)
    % imagesc(order_inds); pause(10e-3);
    
    ptr = 1;
%     figure(); hold on; caxis([1 N_frames]);
    while(ptr <= N_wfms)
        D4(:,order_inds(tx_id(ptr), rx_id(ptr)), tx_id(ptr), rx_id(ptr)) = D(:,ptr);
        
        order_inds(tx_id(ptr),rx_id(ptr)) = order_inds(tx_id(ptr),rx_id(ptr)) + 1;
        ptr = ptr + 1;
%         imagesc(order_inds); pause(10e-3);
    end
    
    % subtract DC - required for xcorr
    D4 = D4 - mean(D4, 1);
    
    % This is a loop to create a matrix 'F' which corresponds to
    % the "flipped" state of a TX-RX-Chirp# pair.
    % It loops through all TX/RX/Frames and puts a logical "1" if
    % the waveform is taken with a flipped sampling ramp and a logical
    % "0" if the waveform is taken with a normal sampling ramp
    % flipped-2047,2046,...,0. normal - 0,1,...,2047
    F = zeros(C.N_tx, C.N_rx, N_frames,'logical');
    
    
    is_flipped = false;
    order_inds  = ones(C.N_tx, C.N_rx);
    
    ptr = 1;
    % follow the sample_order from C
    while(ptr <= N_wfms)
        for jj = 1:size(C.sample_order,1)
            tx_id = C.sample_order(jj,1);
            rx_id = C.sample_order(jj,2);
            
            F(tx_id, rx_id, order_inds(tx_id, rx_id)) = is_flipped;
            is_flipped = ~is_flipped;
            order_inds(tx_id, rx_id) = order_inds(tx_id,rx_id)+1;
            ptr = ptr + 1;
        end
    end
    
    
    
    % un-flip the flipped waveforms
    % calculate the "step" (appears in real data) by looking at waveform 1s
    base          = squeeze(D4(:,~squeeze(F(1, 1, :)), 1, 1));
    toflip        = squeeze(D4(:, squeeze(F(1, 1, :)), 1, 1));
    [xC,lags]     = xcorr(mean(base,2), flip(mean(toflip,2),1), 2);
    [~,ind_shift] = max(xC);
    step          = lags(ind_shift);
    
    % iterate over TX-RX pairs and flip what must be flipped.
    for ii = 1:C.N_tx
        for jj = 1:C.N_rx
            flipped_inds = squeeze(F(ii, jj, :));
            Flipped      = squeeze(D4(:, flipped_inds, ii, jj));
            Flipped      = flip(Flipped, 1);
            
            if(step >= 0)
                D4(:, flipped_inds, ii, jj) = ...
                    cat(1, ...
                        zeros(step,size(Flipped,2)),...
                        Flipped(1:end-step,:));
            elseif(step < 0)
                D4(:, flipped_inds, ii, jj) = ...
                    cat(1, ...
                        Flipped(abs(step)+1:end,:), ...
                        zeros(abs(step),size(Flipped,2)));
            end
        end
    end    
end