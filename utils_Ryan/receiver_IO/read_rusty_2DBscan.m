function [txrx, D2,C,B2,x] = read_rusty_2DBscan(fn,is_v2)
    % Samuel Wagner
    % Feb 11, 2022
    % read_rusty_2DBscan(filename)
    %
    % This is a function intended to read a "rusty" B-scan from a file
    % By "rusty", I am saying that the file was dumped via a Rust program
    % using the special format I developed 
    % 
    % description of the format:
    % first byte u8 in the file refers to the contents of the rest of the file
    % the first byte (u8) is called `content_info`
    % if `content_info` == 0
    %    => the file contains `metadata` and `data` (in that order)
    % if `content_info` == 1
    %    => the file contains `metadata`, `data`, and `cal` (in that order)
    % if `content_info` == 2
    %    => the file contains `metadata`, `data`, and `bkgr` (in that order)
    % if `content_info` == 3
    %    => the file contains `metadata`, `data`, `cal`, and `bkgr` (in that order)
    %
    % note: all data is stored as big endian (BE) or 'b'
    %
    % `metadata` is two u32  immediately after `content_info`
    %    1st u32 of `metadata` => N_t or number of fast-time samples
    %    2nd u32 of `metadata` => N_x or number of slow-time samples
    %    
    % `txrx` is a N_x*1 vector of uint8 following `metadata`
    % 
    % `data` is the raw 2D array of the B-scan
    %   `data` may or may not be time-aligned with a calibration waveform
    %   internal data is stored in row-major order as f32. the total number of
    %   samples in data will be N_t * N_x.
    %
    % `cal` is the calibration waveform vector, if it exists
    %   `cal` will be size N_t * 1 and stored as f32. 
    %
    % `bkgr` is the background-removed version of data, if it exists
    %   If `bkgr` exists, it will be time-aligned with `data` and `cal`. 
    %   It will be the same size as `data` (N_t*N_x) and stored as f32
    %   in the same method as data.

    
    % ensure file exists
    assert(isfile(fn), "File " + fn + " does not exist. Make sure that the filename points to the exact path.");

    if(~exist('is_v2','var'))
       is_v2 = false; 
    end
    
    endianness = 'b';
    
    % open the file and let them know
    fp = fopen(fn,'rb');   
    
    % read the first byte `content_info`
    content_info = fread(fp,[1 1], 'uint8',endianness);
    
    % read in the next u32 `metadata`
    N_t = fread(fp, [1 1], 'uint32',endianness);
    N_x = fread(fp, [1 1], 'uint32',endianness);
    
    txrx_u8   = uint8(fread(fp, [N_x, 1], 'uint8', endianness));
    txrx      = zeros(N_x,2);
    txrx(:,1) = bitand(bitshift(txrx_u8,-4),uint8(15),'uint8');  
    txrx(:,2) = bitand(txrx_u8,uint8(15),'uint8');  
    
    % now, depending on content_info, read in different pieces of data.
    D2 = []; % container for `data`
    C  = []; % container for `cal`
    B2 = []; % container for `bkgr`
    x  = []; % container for x positions
    file_info_str = fn + "found containing: ";
    
    % for some reason, there is a junk vector at the start.?
    %fread(fp,[N_t, 1],'float32',endianness);
    
    switch content_info
        case uint8(0)  % contains only data
            D2 = fread(fp,[N_t, N_x], 'float32', endianness);
            if(is_v2)
                x  = fread(fp,[N_x, 1], 'float32', endianness);
            end
            file_info_str = file_info_str + "meta["+num2str(N_t)+"*"+num2str(N_x)+"],data";
        case uint8(1)  % contains data,cal
            D2 = fread(fp,[N_t, N_x], 'float32', endianness);
            C  = fread(fp,[N_t, 1],   'float32', endianness);
            if(is_v2)
                x  = fread(fp,[N_x, 1], 'float32', endianness);
            end
            file_info_str = file_info_str + "meta["+num2str(N_t)+"*"+num2str(N_x)+"],data,cal";
        case uint8(2)  % contains data,bkgr
            D2 = fread(fp,[N_t, N_x], 'float32', endianness);
            B2 = fread(fp,[N_t, N_x], 'float32', endianness);
            if(is_v2)
                x  = fread(fp,[N_x, 1], 'float32', endianness);
            end
            file_info_str = file_info_str + "meta["+num2str(N_t)+"*"+num2str(N_x)+"],data,bkgr";
        case uint8(3)  % contains data,cal,bkgr
            D2 = fread(fp,[N_t, N_x], 'float32', endianness);
            C  = fread(fp,[N_t, 1],   'float32', endianness);
            B2 = fread(fp,[N_t, N_x], 'float32', endianness);
            if(is_v2)
                x  = fread(fp,[N_x, 1], 'float32', endianness);
            end
            file_info_str = file_info_str + "meta["+num2str(N_t)+"*"+num2str(N_x)+"],data,cal,bkgr";            
        otherwise
            error('content_info was not an expected format. This data file may not be rusty?');
    end
    
    disp(file_info_str);
    
    fclose(fp);
end