classdef Scan < handle
    properties
        % Main data containers
        D        (:,:) single; % Raw Data of time by x (rows - time, columns- x)
        D_BKGR   (:,:) single; % BKGR data
        D_MF     (:,:) single; % matched filter data
        I        (:,:) single; % migrated image (from D_BKGR)
        I_MF     (:,:) single; % migrated image (from D_MF)
        MF_WFM   (1,1) Pulse;  % contains the matched wfm pulse
        C        (:,1) single; % calibration waveform
        
        % target params
        N_targets (1,1) single;
        Targets   (:,1) Target2D;
        
        % Imaging parameters
        xC          (:,1) single; % calculation points (x) for image
        zC          (:,1) single; % calculation points (z) for image
        h           (1,1) single; % height of scan above ground
        er          (1,1) single; % relative permittivity of ground
        tx_pos      (1,1) single; % tx distance (relative to x = 0)
        rx_pos      (1,1) single; % rx distance (relative to x = 0)
        pixel_width (1,1) single; % pixel width for a single image
        
        % Scan parameters
        x   (:,1) single; %x vector
        t   (:,1) single; %time vector
        f   (:,1) single; %frequency vector
        N_x (1,1) single;
        N_t (1,1) single;
        Ts  (1,1) single; %time sampling period
        Fs  (1,1) single; %frequency sampling rate
        Kx  (1,1) single; %spatial sampling rate
        kx  (:,1) single; %spatial frequency vector
        dx  (1,1) single; %spatial sampling period
        has_bit1_flip_indicator (1,1) logical;
        has_bkgr (1,1) logical;
        has_cal (1,1) logical;
        
        name string;      % title of scan
        log  string;        % processing log
        has_parallel_proc_toolbox (1,1) logical;
    end
    
    methods
        % methods
        % S = Scan()
        %   | constructor
        % S.read_raw_bin_data(fn)
        %   | reads raw bin data at fn
        % S.setup_time_axis(Ts, N_t)
        %   | creates the time/freq vectors
        % S.setup_x_axis(x_max)
        %   | creates the x vector based on max distance travelled
        % S.filter_F(use_BKGR, f_min, f_max, N_FIR)
        %   | filters D (use_BKGR=0) or D_BKGR(use_BKGR=1) using an
        %   | N_FIR fir filter. can drop N_FIR as a variable default 50
        % S.remove_jitter(is_interp_enable, interp_ratio)
        %   | will remove most ToA jitter
        % S.BKGR_PCA(L)
        %   | performs principal component analysis to remove
        %   | bacgrkound with parameter L
        % S.setup_imaging_params(pixel_width,zlims,xlims,h,er,tx_pos,rx_pos)
        %   | will set the correct imaging params based on pixel_width
        % S.matched_filter(P, take_positive_only)
        %   | performs a matched-filter on D_BKGR
        % S.migrate(progress)
        %   | migrates D_BKGR. progress = optional progress bar
        % S.migrate_MF(progress)
        %   | migrates D_MF. progress = optional progress bar
        % S.plot(enable_vec)
        %   | plots the Scan object, with optional parameters
        % S.find_targets(target_threshold_dB)
        %
        % S.filter_X(use_BKGR, kx_max, N_FIR)
        
        % constructor
        function obj = Scan()
            obj.log = "";
            obj.N_x = 0;
            obj.N_t = 0;
            obj.Fs  = 0;
            obj.Ts  = 0;
            obj.Kx  = 0;
            obj.dx  = 0;
            obj.has_bit1_flip_indicator = false;
            obj.has_cal = false;
            obj.has_bkgr = false;
        end
        
        function obj = read_raw_bin_data(obj, fn, has_bit1_flip_indicator)
            assert(isfile(fn), fn + " does not exist. Check paths.");
            assert(obj.N_t > 0, "setup_time_axis(Ts, N_t) must be called first.");
            [~,obj.name,~]  = fileparts(fn);
            obj.name  = strrep(obj.name, '_','-');

            if(~exist('has_bit1_flip_indicator','var'))
                obj.has_bit1_flip_indicator = false;
            else
                obj.has_bit1_flip_indicator = has_bit1_flip_indicator;
            end
            
            obj.D     = read_rx2_bindata_2DBscan(fn, obj.N_t, obj.has_bit1_flip_indicator);
            obj.N_x   = size(obj.D,2);
            
            obj.log   = append(obj.log, "read; ");
        end
        
        function obj = read_rusty_bin_data(obj, fn,is_v2)
            assert(isfile(fn), fn + " does not exist. Check paths.");
            if(~exist('is_v2','var'))
                is_v2 = false; 
            end
            [~,D_in,C_in,B_in,x_vec] = read_rusty_2DBscan(fn,is_v2);
            
            if(obj.N_t > 0)
                assert(size(D_in,1) == obj.N_t, "Error: Mismatch between N_t and read Data (" + num2str(obj.N_t) + ") vs (" +num2str(size(D_in,1)) +")");
            else
                obj.N_t = size(D_in,1);
            end
            
            obj.N_x = size(D_in,2);
            obj.D   = D_in;
            
            if(numel(C_in)>0)
                obj.C = C_in;
                obj.has_cal = true;
            end
            if(numel(B_in)>0)
                obj.D_BKGR=B_in;
                obj.has_bkgr = true;
            end
            if(numel(x_vec)>0)
               obj.x = x_vec;
               % future self - write smoothing implementation here!!
            end

            
            [~,obj.name,~] = fileparts(fn);
            obj.name  = strrep(obj.name, '_','-');
            obj.log   = append(obj.log, "read_rusty; ");
        end
        
        
        function obj = setup_time_axis(obj, Ts, N_t)
            obj.N_t = N_t;
            obj.Ts  = Ts;
            obj.Fs  = 1/obj.Ts;
            obj.t   = (0:(obj.N_t - 1)).*obj.Ts;
            obj.f   = linspace(-obj.Fs/2, obj.Fs/2, obj.N_t);
        end
        
        function obj = setup_x_axis(obj, x_max)
            obj.x   = linspace(0,x_max, obj.N_x);
            obj.dx  = obj.x(2) - obj.x(1);
            obj.Kx  = 1/obj.dx;
            obj.kx  = linspace(-obj.Kx/2,obj.Kx/2,obj.N_x);
        end
        
        function obj = filter_F(obj, use_BKGR, f_min, f_max, N_FIR)
            % input checking
            assert(f_max>f_min,      "f_max must be greater than f_min");
            assert(obj.Fs > 0,       "must call setup_time_axis(Ts, N_t) first");
            assert(f_max < obj.Fs/2, "f_max must be less than Fs/2");
            
            if(~exist('N_FIR','var'))
                N_FIR = 50;
            end
            assert(mod(N_FIR,2)==0,"N_FIR must be even.");
            
            
            TF = freqz(fir1(N_FIR, [f_min f_max]./(obj.Fs/2)),1, obj.f, obj.Fs);
            if(size(TF,1)==1)
                TF = TF.';
            end
            
            if(use_BKGR)
                assert(size(obj.D_BKGR,1)>0,"must call a BKGR function first before filter_F if use_BKGR is enabled.")
                Df = fftshift(fft(obj.D_BKGR,[],1),1).*repmat(TF,1,obj.N_x);
                obj.D_BKGR = real(ifft(ifftshift(Df,1),[],1));
            else
                Df = fftshift(fft(obj.D,[],1),1).*repmat(TF,1,obj.N_x);
                obj.D = real(ifft(ifftshift(Df,1),[],1));
            end
            
            obj.D = circshift(obj.D, -N_FIR/2, 1);
            obj.log = append(obj.log, "filteredF:"+f_min/1e9+"GHz-"+f_max/1e9+"GHz|N_FIR="+N_FIR+"; ");
        end
        
        function obj = filter_X(obj, use_BKGR, kx_max, N_FIR)
            % input checking
            assert(obj.Kx > 0,       "must call setup_x_axis(x_max) first");
            assert(kx_max < obj.Kx/2, "kx_max must be less than Ks/2");
            
            if(~exist('N_FIR','var'))
                N_FIR = 50;
            end
            assert(mod(N_FIR,2)==0,"N_FIR must be even.");
            
            
            TF = freqz(fir1(N_FIR, kx_max./(obj.Kx/2)),1, obj.kx, obj.Kx);
            if(size(TF,2)==1)
                TF = TF.';
            end
            
            if(use_BKGR)
                assert(size(obj.D_BKGR,1)>0,"must call a BKGR function first before filter_F if use_BKGR is enabled.")
                Dx = fftshift(fft(obj.D_BKGR,[],2),2).*repmat(TF,obj.N_t,1);
                obj.D_BKGR = real(ifft(ifftshift(Dx,2),[],2));
            else
                Dx = fftshift(fft(obj.D,[],2),2).*repmat(TF,obj.N_t,1);
                obj.D = real(ifft(ifftshift(Dx,2),[],2));
            end
            
            obj.D = circshift(obj.D, -N_FIR/2, 2);
            obj.log = append(obj.log, "filteredX:"+kx_max+"m^-1|N_FIR="+N_FIR+"; ");
        end
        
        function obj = BKGR_cal(obj, cal, samples_pm_to_search)
            sample_vecs = -samples_pm_to_search:samples_pm_to_search;
            N_energies  = numel(sample_vecs);
            energies    = zeros(N_energies,1);
            
            D_BKG = mean(obj.D,2);
            
            if(obj.has_parallel_proc_toolbox)
                parfor ii = 1:N_energies
                    energies(ii) = sum((D_BKG - circshift(cal,sample_vecs(ii))).^2); 
                end
                
            else
                for ii = 1:N_energies
                    energies(ii) = sum((D_BKG - circshift(cal,sample_vecs(ii))).^2); 
                end 
            end
            
            [~,ind_best_shift] = min(energies);
            cal_wfm            = circshift(cal,sample_vecs(ind_best_shift));
            
            obj.D_BKGR = obj.D - cal_wfm;
            
%             figure()
%             hold on;
%             plot(mean(obj.D,2));
%             plot(cal_wfm,'r');
        end
        
        function obj = remove_jitter(obj, is_interp_enable, interp_ratio)
            % wrapper function to remove the ToA jitter from a scan
            
            % input checking
            if(~exist('is_interp_enable','var'))
                is_interp_enable = 0;
                interp_ratio     = 1;
            end
            if(isinteger(interp_ratio))
               interp_ratio = single(interp_ratio); 
            end
            
            assert(interp_ratio >= 1, "interp_ratio must be greater or equal to 1");
            assert(is_double_single_logical_scalar(is_interp_enable), "is_interp_enable must be a double/single/logical scalar");
            assert(is_double_single_2d_matrix(obj.D), "obj.D must be a 2-D matrix of singles or doubles");
            assert(is_double_single_scalar(interp_ratio), "interp_ratio must be a single or double scalar");
            assert(is_integer_scalar(interp_ratio),"interp_ratio must be a single or double scalar integer");
            
            is_interp_enable = logical(is_interp_enable);
            obj.D   = remove_GPR_Bscan_jitter(obj.D, is_interp_enable, interp_ratio);
            obj.log = append(obj.log, "dejittered interp?="+is_interp_enable+"; ");
        end
        
        function obj = BKGR_PCA(obj, L)
            % function to remove the background from the scan using
            % principal component analysis
            % inputs_______________________________________________
            %   L: an integer - the number of principal components to remove 
            % outputs______________________________________________
            %   none
            
            % input checking
            assert(is_double_single_2d_matrix(obj.D),       "The data matrix (obj.D) must be populated as a single or double 2-D matrix");
            assert(is_integer_scalar(L) && L>0,"Input L must be an scalar integer greater than zero.");
            
            % perform a singular value decomposition
            [U, S, V]     = svd(obj.D,'econ');
            
            % find o - the relative strength of the principal components
            o             = diag(S);
            
            % reconstruct the matrix, inserting zeros for the first `L`
            % singular values. save as D_BKGR.
            obj.D_BKGR    = U*diag([zeros(L,1); o(L+1:end)])*V';
            
            obj.log       = append(obj.log,"BKGR_PCA with L="+num2str(L)+"; ");
        end
        
        function obj = chop_time_axis(obj, chop_time_env_threshold, ignore_time_lessthan)
            % function to chop the "before" time out of the time axis
            % inputs_______________________________________________
            %   chop_time_env_threshold : a threshold number from 0 to 1 at
            %        which to call the pulse t = 0 "start". typ approx 0.075. needs adjustment
            %        for different pulses
            %   ignore_time_lessthan    : if there are glitches before the coupling signal
            %                             then use this parameter. if there
            %                             is an event that triggers the
            %                             chop_time_env_threshold before
            %                             ignore_time_lessthan, this
            %                             function will ignore the event.
            % outputs______________________________________________
            %   none
            
            if(~exist('ignore_time_lessthan','var'))
               ignore_time_lessthan = 0e-9; 
            end
            % input checking  
            assert(is_double_single_scalar(chop_time_env_threshold), "chop_time_env_threshold must be a double/single scalar");
            assert(chop_time_env_threshold > 0 && chop_time_env_threshold <1, "chop_time_env_threshold must be between 0 and 1, typ 0.05 - 0.1");
            assert(is_double_single_scalar(ignore_time_lessthan), "ignore_time_lessthan must be a double/single scalar");
            assert(ignore_time_lessthan >= 0 && ignore_time_lessthan <= max(obj.t), "ignore_time_lessthan must be between 0 and max(t)");
            assert(is_double_single_scalar(obj.tx_pos),"tx_pos must be a single or double vector. It is not. Perhaps it is not constructed (setup_imaging_params()).");
            assert(is_double_single_vector(obj.t),"The time vector must be a single or double vector. It is not. Perhaps it is not constructed (setup_time_axis()).");
            assert(is_double_single_2d_matrix(obj.D),"The data matrix (obj.D) must be populated as a single or double 2-D matrix");
            assert(is_double_single_scalar(obj.Ts),"Sampling period (obj.Ts) must be a single/double scalar");
            
            % reduce the verbosity...
            env_thresh  = chop_time_env_threshold;
            t_ignore_lt = ignore_time_lessthan;
            
            % Estimate the "background" (BK)
            BK       = normalize(mean(obj.D,2));
            
            % To perform t=0 detection, we take the cumulative integral of the
            % background energy (cumtrapz), and then take the derivative
            % when the derivative starts, that means that we have a pulse
            % starting. this is the coupling waveform. by setting
            % chop_time_env_threshold, we can control when on the slope we
            % decide to call near "t=0"
            BKE_grad = normalize(gradient(cumtrapz(BK.^2)));
            
            % find the FIRST index where the cumulative-energy slope is
            % greater than the threshold AND time is greater than
            % ignore-time
            t0_ind = find(BKE_grad>env_thresh & (obj.t > t_ignore_lt), 1, 'first');

            % control for the distance between tx and rx at t=0, as the
            % direct coupling wave must travel some distance.
            t0_ind = t0_ind + round((abs(obj.tx_pos - obj.rx_pos)./3e8)./obj.Ts);
            
            % slice the obj.D matrix based on the new t = 0
            obj.D = obj.D(t0_ind:end,:);
            
            % if obj.D_BKGR exists, slice it too.
            if(size(obj.D_BKGR,1)>0)
                obj.D_BKGR = obj.D_BKGR(t0_ind:end,:);
            end
            
            obj.setup_time_axis(obj.Ts, size(obj.D,1));

            obj.log = append(obj.log,"Timechop at time: " + obj.t(t0_ind).*1e9 + "ns; ");
        end
        
        function obj = setup_imaging_params(obj,pixel_width,zlims,xlims,h,er,tx_pos,rx_pos)
            % function to validate and set the Scan object parameters for imaging
            % inputs_______________________________________________
            %   pixel_width: width of a pixel in the image.
            %   zlims      : [z_min, z_max] imaging domain in depth
            %   xlims      : [x_min, x_max] imaging domain in depth
            %   h          : antenna phase center height above ground
            %   er         : relative permittivity (est) of soil  / gnd
            %   tx_pos     : tx antenna position relative to scan axis center
            %   rx_pos     : rx antenna position relative to scan axis center
            % outputs______________________________________________
            %   none

            %input validation
            is_double_single_2vec    = @(x) (strcmp(string(class(x)),"double")||strcmp(string(class(x)),"single")) && numel(x) == 2;
            
            assert((h>0)  && is_double_single_scalar(h),                    "h must be a double/single scalar greater than 0");
            assert((er>1) && is_double_single_scalar(er),                   "er must be a double/single scalar greater than 1");
            assert((pixel_width>0) && is_double_single_scalar(pixel_width), "pixel_width must be a double/single scalar greater than 0");
            assert(is_double_single_2vec(zlims),                            "zlims must be a two-element vector containing the upper and lower limits for depth calculation");
            assert(min(zlims)>=0,                                           "min(zlims) must be greater than 0");
            assert(is_double_single_2vec(xlims),                            "xlims must be a two-element vector containing the upper and lower limits for scan axis calculation");
            assert(is_double_single_scalar(tx_pos),                         "tx_pos must a double/single scalar");
            assert(is_double_single_scalar(rx_pos),                         "rx_pos must a double/single scalar");
            
            % set up the imaging variables
            obj.xC          = min(xlims):pixel_width:max(xlims);
            obj.zC          = min(zlims):pixel_width:max(zlims);
            obj.pixel_width = pixel_width;
            obj.h           = h;
            obj.er          = er;
            obj.tx_pos      = tx_pos;
            obj.rx_pos      = rx_pos;
        end
        
        function obj = migrate(obj,progress)
            if(~exist('progress','var'))
                progress = 0;
            end
            obj.I   = BScan_kirchoff_migration3(obj.D_BKGR,...
                obj.h, obj.er, obj.t, obj.x, obj.tx_pos, obj.rx_pos, obj.xC, obj.zC, obj.has_parallel_proc_toolbox, progress);
            obj.log = append(obj.log,"migrated (normal); ");
        end
        
        function obj = matched_filter(obj, P, take_positive_only)
            % wfm - waveform to use as reference
            % positive - do you want to throw out negative xcorrs?
            
            % input checking
            is_one_pulse                    = @(x)  strcmp(string(class(x)),"Pulse")  && (numel(x) == 1);

            assert(is_one_pulse(P), "P must be a single Pulse object.");
            assert(is_double_single_scalar(P.Ts), "P.Ts must be a single or double scalar");
            assert(is_double_single_vector(P.x) && is_double_single_vector(P.t), "P.x and P.t must be filled with a double/single vector.");
            assert(is_double_single_logical_scalar(take_positive_only), "take_positive_only must be a logical or double/single scalar");
            
            % to do a matched-filter, flip the input waveform
            to_conv    = flip(P.x);
            
            % check that the sampling period is the same. if not,
            % interpolate to_conv to the correct sampling period
            if(P.Ts ~= obj.Ts)
                to_conv = interp1(P.t - min(P.t), to_conv, obj.t);
            end
            
            % create a new pulse object based on to_conv to save the MFWFM
            obj.MF_WFM = Pulse(flip(to_conv), obj.t);
            
            % perform a 2-D convolution and save the central portion
            obj.D_MF = conv2(obj.D_BKGR,to_conv,'same');
            
            % if the take_positive_only flag is set, undo the negative
            % xcorr by setting it to zero.
            logmsg = "matched-filtered; ";
            if(take_positive_only)
                [~,max_ind] = max(abs(obj.D_MF(:)));
                if(obj.D_MF(max_ind) < 0)
                    obj.D_MF(obj.D_MF>0)=0;
                    obj.D_MF = obj.D_MF.*-1;
                else
                    obj.D_MF(obj.D_MF<0)=0;
                end
                logmsg = "matched-filtered (positive only); ";
            end
            
            % update log
            
            obj.log = append(obj.log,logmsg);
        end
        
        function obj = migrate_MF(obj,progress)
            if(~exist('progress','var'))
                progress = 0;
            end
            obj.I_MF   = BScan_kirchoff_migration3(obj.D_MF,...
                obj.h, obj.er, obj.t, obj.x, obj.tx_pos, obj.rx_pos, obj.xC, obj.zC, obj.has_parallel_proc_toolbox, progress);
            obj.log = append(obj.log,"migrated (matched-filter); ");
        end
        
        function obj = find_targets(obj, target_threshold_dB, max_Z, max_X, target_MPP, target_MPS)
            assert(target_threshold_dB < 0, "target_threshold_dB should be in dB and less than zero. typ -3.")
            if(numel(obj.I_MF)<1)
                IMG = obj.I;
            else
                IMG = obj.I_MF;
            end
            
            obj.Targets   = find_2D_targets(obj, IMG, target_threshold_dB, max_X, max_Z, target_MPP, target_MPS);
            obj.N_targets = numel(obj.Targets);
            obj.log = append(obj.log,"found-targets; ");
        end
        
        function obj = discriminate_targets(obj,method,arg)
            implemented_methods = [...
                "snr",...
                "maxarea"...
                ];
            method = lower(method);
            assert(ismember(method,implemented_methods),"method must be one of the following: " + strjoin(implemented_methods, ", "));
            
            if(strcmp(method, "SNR"))
                if(arg>0)
                    arg=-arg;
                end
                min_allowable_SNR = max([obj.Targets.SNR]) + arg;
                obj.Targets([obj.Targets.SNR] < min_allowable_SNR) = [];
            end
            
            if(strcmp(method, "MaxArea"))
                max_allowable_area = min([obj.Targets.area].*arg);
                obj.Targets([obj.Targets.area] > max_allowable_area) = [];
            end
            
            obj.N_targets = numel(obj.Targets);
        end
        
        function obj = update_target_SNRs(obj, target_estimation_method, noise_estimation_method)
            if(numel(obj.I_MF)<1)
                IMG = obj.I;
            else
                IMG = obj.I_MF;
            end
            target_estimation_method = lower(target_estimation_method);
            noise_estimation_method  = lower(noise_estimation_method);
            
            % implemented method definitions
            implemented_noise_methods =  [...
                "removeall", ...
                "removeselfonly"...
                ];
            implemented_target_methods = [...
                "peak",...
                "avgenergy",...
                "totalenergy"...
                ];
            
            % input checking
            assert(ismember(noise_estimation_method,implemented_noise_methods), ...
                "noise_estimation_method must be one of the following: " + strjoin(implemented_noise_methods, ", "));
            assert(ismember(target_estimation_method,implemented_target_methods), ...
                "target_estimation_method must be one of the following: " + strjoin(implemented_target_methods, ", "));
            
            
            % estimate the targets
            if(strcmp(target_estimation_method, "totalenergy"))
                Target_signal = [obj.Targets.energy];
            elseif(strcmp(target_estimation_method, "avgenergy"))
                Target_signal = [obj.Targets.energy]./[obj.Targets.area];
            elseif(strcmp(target_estimation_method, "peak"))
                Target_signal = [obj.Targets.peak];
            end
            
            
            % estimate the noises
            if(strcmp(noise_estimation_method, "removeall"))
                for ii = 1:obj.N_targets
                    [~,x_ind] = min(abs(obj.xC-obj.Targets(ii).upper_left_x));
                    [~,z_ind] = min(abs(obj.zC-obj.Targets(ii).upper_left_z));
                    Nz        = obj.Targets(ii).Nz;
                    Nx        = obj.Targets(ii).Nx;
                    IMG(z_ind + (0:Nz-1), x_ind + (0:Nx-1)) = nan;
                end
                IMG             = IMG(:);
                IMG(isnan(IMG)) = [];
                Target_noise    = std(IMG).*ones(obj.N_targets,1);
            elseif(strcmp(method, "removeselfonly"))
                Target_noise = zeros(obj.N_targets, 1);
                for ii = 1:obj.N_targets
                    IMG2 = IMG;
                    % target under test
                    [~,x_ind] = min(abs(obj.xC-obj.Targets(ii).upper_left_x));
                    [~,z_ind] = min(abs(obj.zC-obj.Targets(ii).upper_left_z));
                    Nz        = obj.Targets(ii).Nz;
                    Nx        = obj.Targets(ii).Nx;
                    IMG2(z_ind + (0:Nz-1), x_ind + (0:Nx-1)) = nan;
                    IMG2             = IMG2(:);
                    IMG2(isnan(IMG2)) = [];
                    Target_noise(ii)    = std(IMG2).*ones(obj.N_targets,1);
                end
            end
            
            % update the SNRs
            for ii = 1:obj.N_targets
                obj.Targets(ii).SNR = Target_signal(ii)/Target_noise(ii);
            end
            
        end
        
        function obj = remove_glitches(obj,threshold)
%             progressbar('glitch removal');
%             for ii = 1:obj.N_x
%                wfm = obj.D(:,ii);
%                diff_wfm = abs(diff(wfm)) > threshold;
%                
%                for jj = 1:obj.N_t - 1
%                    if(diff_wfm(jj))
%                         %replace the larger abs one
%                         left  = wfm(jj);
%                         right = wfm(jj+1);
%                         
%                         if(abs(left)>abs(right))
%                             wfm(jj)   = right;
%                         else
%                             wfm(jj+1) = left;
%                         end    
% 
%                    end
%                end
%                progressbar(ii/obj.N_x);
%                obj.D(:,ii) = wfm;
%                 obj.D(abs(obj.D) > threshold) = 0; 
%             end
            
                obj.D(1:10,:) = 0;
                obj.D(obj.N_t - 10:end,:) = 0;
        end
        
        function plot(obj, enable_vec)
            enable_vec = logical(enable_vec);
            N_plots = sum(enable_vec);
            N_max_plots = numel(enable_vec);
            assert(N_plots>0, "enable_vec must be a logical vector corresponding to whether or not you want to plot.. [raw bkgr mf image image_mf]");
            % raw bkgr mf image image_mf
            plot_ptr = 0;
            
            figure();
            for ii = 1:N_max_plots
                if(enable_vec(ii))
                    plot_ptr = plot_ptr + 1;
                    ax(plot_ptr) = subplot(1, N_plots, plot_ptr);
                    plot_targs = false;
                    switch ii
                        case 1
                            x_vec = obj.x;
                            y_vec = obj.t.*1e9;
                            I_mat = obj.D;
                            cmap  = gray(2^12);
                            tit   = obj.name + ":Raw";
                            xL    = "x (m)";
                            yL    = "t (ns)";
                            ydir  = 'reverse';
                            clims = [min(I_mat(:)) max(I_mat(:))].*0.975;
                        case 2
                            x_vec = obj.x;
                            y_vec = obj.t.*1e9;
                            I_mat = obj.D_BKGR;
                            cmap  = gray(2^12);
                            tit   = obj.name + ":BKGR";
                            xL    = "x (m)";
                            yL    = "t (ns)";
                            ydir  = 'reverse';
                            clims = [min(I_mat(:)) max(I_mat(:))].*0.975;
                        case 3
                            x_vec = obj.x;
                            y_vec = obj.t.*1e9;
                            I_mat = obj.D_MF;
                            cmap  = gray(2^12);
                            tit   = obj.name + ":MF";
                            xL    = "x (m)";
                            yL    = "t (ns)";
                            ydir  = 'reverse';
                            clims = [min(I_mat(:)) max(I_mat(:))].*0.975;
                        case 4
                            x_vec = obj.xC;
                            y_vec = obj.zC.*1e2;
                            I_mat = 10.*log10(normalize(obj.I));
                            cmap  = jet;
                            tit   = obj.name + ":Image";
                            xL    = "x (m)";
                            yL    = "z (cm)";
                            ydir  = 'reverse';
                            clims = [-10 0];
                        case 5
                            x_vec = obj.xC;
                            y_vec = obj.zC.*1e2;
                            I_mat = 10.*log10(normalize(obj.I_MF));
                            cmap  = jet;
                            tit   = obj.name + ":ImageMF";
                            xL    = "x (m)";
                            yL    = "z (cm)";
                            ydir  = 'reverse';
                            clims  = [-6 0];
                        case 6
                            x_vec = obj.xC;
                            y_vec = obj.zC.*1e2;
                            if(size(obj.I_MF,1)~= numel(obj.zC))
                                I_mat = 10.*log10(normalize(obj.I));
                            else
                                I_mat = 10.*log10(normalize(obj.I_MF));
                            end
                            cmap  = jet;
                            tit   = obj.name + ":ImageAnnotate";
                            xL    = "x (m)";
                            yL    = "z (cm)";
                            ydir  = 'reverse';
                            clims  = [-10 0];
                            plot_targs = true;
                    end
                    
                    colormap(ax(plot_ptr), cmap); hold on;
                    imagesc(x_vec, y_vec, I_mat);
                    ax(plot_ptr).YDir = ydir;
                    title(tit);
                    xlabel(xL);
                    ylabel(yL);
                    colorbar;
                    ylim([min(y_vec) max(y_vec)]);
                    xlim([min(x_vec) max(x_vec)]);
                    caxis(clims);
                    if(plot_targs)
                        for jj = 1:obj.N_targets
                            obj.Targets(jj).draw(ax(plot_ptr),1,1e2);
                        end
                    end
                end
            end
        end
        
    end
end