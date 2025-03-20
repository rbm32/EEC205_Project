function SS = remove_GPR_Bscan_jitter(D, is_interp_enable, interp_ratio)
% remove_GPR_skyshot_jitter.m
% Samuel Wagner, UCD ECE MML, Jan. 16, 2022
%
% This is a function intended to "de-jitter" the waveforms in the
% matrix Skyshot
%
% inputs ______
% Skyshot - an N_t x N_tx x N_rx matrix of singles or doubles
%           containing the to-be-calibrated coupling signal of the array
% is_interp_enable - a boolean to whether the data in D4(:,) can be
%                    interpolated for better de-jittering accuracy
% interp_multiple  - the ratio of interpolated points to normal points
%                    for example, interp_multiple = 10 turns a time
%                    sampling vector from 0:Ts:(N_t-1)*Ts to 
%                    0:(Ts/10):(N_t-1)*Ts. Does not need to be supplied
%                    if is_interp_enable is false


    %% input checking ____________________________________________________
    % Skyshot checks
    assert(numel(size(D))==2, ...
        "D must be a 2-D matrix with dim1: time, dim2: frames");
    assert(strcmp(string(class(D)),"double")||strcmp(string(class(D)),"single"),...
        "D4 must be of type single or double.");
    
    size_SS = size(D);

    % is_interp_enable checks
    if(~strcmp(class(is_interp_enable),"logical"))
        assert(ismember(is_interp_enable,[0 1]),...
            "is_interp_enable should be a logical 0 or 1");
        is_interp_enable = logical(is_interp_enable);
    end

    % interp_ratio checks
    if(~is_interp_enable)
        assert(round(interp_ratio)==interp_ratio && interp_ratio >= 1, ...
            "interp_multiple must be a positive integer great than or equal to 1.")
    else
       interp_ratio = 1;
    end

    %% perform an interpolation
    maxlag      = 200*interp_ratio;
    N_frames    = size_SS(2);

    % interpolate the 4-D matrix
    D = linear_interpolate(D, 1, interp_ratio);

    % find the reference "waveform" foreach.
    env_max_inds = zeros(N_frames, 1);

    for kk = 1:N_frames
        [~,env_max_inds(kk)] = max(abs(hilbert(squeeze(D(:,kk)))));
    end
    [~, refwfm_frame_ind] = min(abs(env_max_inds - median(env_max_inds)));

    
    % loop through to do a proper interpolation for EACH wfm.
    ref_wfm = D(:, refwfm_frame_ind);

    for kk = 1:N_frames
        wfm                 = squeeze(D(:,kk));
        [xC, lags]          = xcorr(ref_wfm, wfm, maxlag);
        [~,max_xc_ind]      = max(xC);
        wfm                 = circshift(wfm, lags(max_xc_ind));
        D(:,kk) = wfm;
    end
    
    SS = D(1:interp_ratio:end,:,:,:);
end