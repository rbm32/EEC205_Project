function D2 = read_rx2_bindata_2DBscan(fn, wfm_len, has_bit1_indicator)
    % this function will take as an input a filename and return the data
    % -generated data - interleaved as [wfm_id][tx_id][rx_id][data1][data2] in 32bit
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
    if(~exist('has_bit1_indicator','var'))
       has_bit1_indicator = false; 
    end
    
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
    wfm_id     = uint8( bitand(bitshift(data_raw,-28), uint32(2^3 -1),'uint32'));
    wfm_id     = wfm_id(1:wfm_len/2:end);
    
    D = int16(reshape([data_upper'; data_lower'], wfm_len, []));
    
    wfm_mode = mode(wfm_id);
    is_wfm_id_mode = wfm_id == wfm_mode;
    D = D(:,is_wfm_id_mode);
    
    N_wfms      = size(D,2);
    
    % subtract DC - required for xcorr
    D2 = single(D) - mean(single(D), 1);
    
    % This is a loop to create a matrix 'F' which corresponds to
    % the "flipped" state of a TX-RX-Chirp# pair.
    % It loops through all TX/RX/Frames and puts a logical "1" if
    % the waveform is taken with a flipped sampling ramp and a logical
    % "0" if the waveform is taken with a normal sampling ramp
    % flipped-2047,2046,...,0. normal - 0,1,...,2047
    if(has_bit1_indicator)
        is_reversed = uint8(bitand(bitshift(data_raw,-31), uint32(1),'uint32'));

        F = is_reversed(2:wfm_len/2:end) == 1;
    else
        F = zeros(N_wfms, 1,'logical');
        is_flipped = false;
        ptr = 1;
        % follow the sample_order from C
        while(ptr <= N_wfms)
            F(ptr)     = is_flipped;
            is_flipped = ~is_flipped;
            ptr = ptr + 1;
        end
    end
    % un-flip the flipped waveforms
    % calculate the "step" (appears in real data) by looking at waveform 1s
    base          = squeeze(D2(:,~F));
    toflip        = squeeze(D2(:, F));
    [xC,lags]     = xcorr(mean(base,2), flip(mean(toflip,2),1), 2);
    [~,ind_shift] = max(xC);
    step          = lags(ind_shift);
    
    % iterate over TX-RX pairs and flip what must be flipped.
    flipped_inds = F;
    Flipped      = squeeze(D2(:, flipped_inds));
    Flipped      = flip(Flipped, 1);

    if(step >= 0)
        D2(:, flipped_inds) = ...
            cat(1, ...
                zeros(step,size(Flipped,2)),...
                Flipped(1:end-step,:));
    elseif(step < 0)
        D2(:, flipped_inds) = ...
            cat(1, ...
                Flipped(abs(step)+1:end,:), ...
                zeros(abs(step),size(Flipped,2)));
    end

end