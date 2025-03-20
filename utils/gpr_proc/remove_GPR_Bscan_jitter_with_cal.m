function D2 = remove_GPR_Bscan_jitter_with_cal(D2, cal, is_interp_enable, interp_ratio)
% remove_4d_GPR_jitter.m
% Samuel Wagner, UCD ECE MML, Jan. 6, 2022
%
% This is a function intended to "de-jitter" the waveforms in the
% matrix D4 with respect to the matrix Skyshot
%
% inputs ______
% D2 - an N_t x N_x  matrix of singles or doubles
%      corresponding to the scan data to be de-jittered
% Cal - an N_t x 1 matrix of singles or doubles
%           containing the calibrated coupling signal of the array
% is_interp_enable - a boolean to whether the data in D4(:,) can be
%                    interpolated for better de-jittering accuracy
% interp_multiple  - the ratio of interpolated points to normal points
%                    for example, interp_multiple = 10 turns a time
%                    sampling vector from 0:Ts:(N_t-1)*Ts to 
%                    0:(Ts/10):(N_t-1)*Ts. Does not need to be supplied
%                    if is_interp_enable is false


    %% input checking ____________________________________________________
    % D4 checks
    assert(is_double_single_2d_matrix(D2), "`D2` must be a 2-D single or double matrix with dim1: time, dim2: frames");
    assert(is_double_single_vector(cal), "`cal` must be a 1-D (vector) single or double with dim1: time");
    

    size_D2 = size(D2);
    size_cal = size(cal);

    assert(size_D2(1) == size_cal(1),...
        "`D2` of size [N_t, N_x] must have matching dimensions with `cal` of size " + ...
        "[N_t, 1]. Check that the dimensions of each `D2` and `cal` fit accordingly.");

    % is_interp_enable checks
    assert(is_double_single_logical_scalar(is_interp_enable),"`is_interp_enable` must be ");

    % interp_ratio checks
    if(~is_interp_enable)
        assert(round(interp_ratio)==interp_ratio && interp_ratio >= 1, ...
            "interp_multiple must be a positive integer great than or equal to 1.")
    else
       interp_ratio = 1;
    end

    %% perform an interpolation
    maxlag      = 200*interp_ratio;
    N_frames    = size_D2(2);


    % interpolate the 4-D matrix
    D2      = linear_interpolate(D2,      1, interp_ratio);
    cal     = linear_interpolate(cal, 1, interp_ratio);

    % loop through to do a proper interpolation for EACH wfm.

    for kk = 1:N_frames
        wfm            = squeeze(D2(:,kk));
        [xC, lags]     = xcorr(cal, wfm, maxlag);
        [~,max_xc_ind] = max(xC);
        wfm            = circshift(wfm, lags(max_xc_ind));
        D2(:,kk) = wfm;
    end


    D2 = D2(1:interp_ratio:end,:,:,:);

end