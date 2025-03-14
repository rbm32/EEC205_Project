function write_wfm_as_rx2_bindata(fp, wfm, tx_id, rx_id, flip)
    % this function will take an input waveform and write it in the format
    % of rx v2 bindata
    % that is - [0101][tx_id (3bit)][rx_id (3bit)][data_upper][data_lower]

    % inputs____
    % fp    - filePTR (must be fopened) to be written.
    % wfm   - waveform (uint16 vector) to be written
    % tx_id - uint8 tx_id - refers to tx antenna at which wfm was taken
    %         spans from 1...N. -1 will be removed automatically
    % rx_id - uint8 rx_id - refers to rx antenna at which wfm was taken
    %         spans from 1...N. -1 will be removed automatically 
    % flip  - logical - should the waveform be flipped (i.e., written reverse)?
    % outputs___
    % this will output a file at path fn 
    assert(min(size(wfm))==1,"wfm must be a vector.");
    assert(strcmp(string(class(wfm)),"uint16"),"wfm must be a vector of uint16 scaled from 0 to 2047.");
    assert(mod(numel(wfm),2)==0,"wfm must have an even number of elements.");
    
    N = numel(wfm);
    header    = bitshift(bitor(uint32(rx_id-1),uint32((tx_id-1))*2^3),22);
    word      = zeros(N/2,1,'uint32');
    
    wordctr=1;
    if(~flip)
        for kk = 1:2:(N-1)
            %select the lower 11 bits from the uint16s
            upper_bits      = uint32(bitand(wfm(kk),uint16(2047)))*uint32(2048);
            lower_bits      = bitand(wfm(kk+1),uint16(2047));
            UL_bits         = bitor(upper_bits,uint32(lower_bits),'uint32');
            word(wordctr)   = bitor(UL_bits,header,'uint32');
            wordctr=wordctr+1;
        end
    else
        for kk = N:-2:2
            %select the lower 11 bits from the uint16s
            upper_bits      = uint32(bitand(wfm(kk),uint16(2047)))*uint32(2048);
            lower_bits      = bitand(wfm(kk-1),uint16(2047));
            UL_bits         = bitor(upper_bits,uint32(lower_bits),'uint32');
            word(wordctr)   = bitor(UL_bits,header,'uint32');
            wordctr=wordctr+1;
        end
    end
    fwrite(fp, word, 'uint32');
end