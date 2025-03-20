classdef Pulse < handle
    % Samuel Wagner, UCD ECE / MML, 2020
    % Object which holds a "Pulse"
    % (any time domain waveform)
    % useful for keeping pulses in check.

    properties
        t   (:,1) double;    % time vector 
        x   (:,1) double;    % time-domain signal
        f   (:,1) double;    % frequency vector
        X   (:,1) double;    % frequency-domain signal
        N   {mustBeNumeric};    % number of samples
        Ts  {mustBeNumeric};    % sampling period
        Fs  {mustBeNumeric};    % sampling rate
        PSD (:,1) double;    % power spectral density
        label {string};         % title of pulse
        vpp (1,1) double;       % something to keep track of...
        ref (1,1) double;       % 
    end

    methods
        % method signatures
        % P = Pulse()
        %     | constructor - don't use much
        % P.normalizeTtoVoltage(Vmax)
        %     | normalizes peak voltage to Vmax
        % P.convertFtoTDomain()
        %     | have F domain, want T domain. must have
        %     | frequency vector f & X defined.
        % P.constructFVector()
        %     | given x & t, create f, Fs
        % P.convertTtoFDomain()
        %     | have T domain, want F domain (& PSD)
        % P.ReadFromFile(fn)
        %     | reads from a .csv file where t is in the 
        %     | first column and x is in the second column
        % P.ReadFromFile_noTime(fn, Fs)
        %     | reads from a .csv file where x is the only column  
        %     | and will generate t,f,etc.. from Fs
        % [flow,fhi] = P.findBandwidth(threshold_dB20)
        %     | will find the bandwidth specified by threshold_db20 in PSD
        % P.chopPulse(threshold)
        %     | will chop the pulse to only the important times based on
        %       its envelope's magnitude (threshold)
        % vpp = P.get_vpp()
        %     | returns the vpp and sets the obj.vpp variable.
        % P.padToGranularity(gran)
        %     | adds zeros at the end of the pulse to be a length of an
        %       integer multiple of gran
        % P.delayWRTRef(Ref)
        %     | delays a pulse such that its peak lines up with another's
        %       peak. Uses envelopes. Ref is input pulse.
        % P.plot()
        %     | plots the time and frequency representations
        
        function obj = Pulse(x,t)
            if(nargin>0)
                obj.N   = numel(x);
                obj.x   = x;
                obj.t   = t;
                obj.Ts  = t(2)-t(1);
                obj.Fs  = 1./obj.Ts;
                obj.f   = linspace(-obj.Fs/2,obj.Fs/2,obj.N);
                obj.X   = fftshift(fft(x));
                obj.PSD = 10.*log10(1./obj.Fs./obj.N.*abs(obj.X).^2);
            end
            obj.label="";
        end

        function obj = normalizeTtoVoltage(obj,Vmax)
            %must have t-domain defined already
            obj.x = obj.x./max(abs((obj.x))).*Vmax;
            obj.convertTtoFDomain();
        end

        function obj = convertFtoTDomain(obj) %have F domain, want T domain
            obj.x  = real(ifft(ifftshift(obj.X)));
            obj.Fs = 2*max(obj.f);
            obj.Ts = 1./obj.Fs;
            obj.N  = numel(obj.x);
            obj.t  = (0:(obj.N-1)).*obj.Ts;
        end

        function obj = constructFVector(obj)
            obj.N  = numel(obj.x);
            obj.Ts = obj.t(2) - obj.t(1);
            obj.Fs = 1./obj.Ts;
            obj.f  = linspace(-obj.Fs/2, obj.Fs/2, obj.N);
        end
        
        function obj = constructTVector(obj)
            obj.N   = numel(obj.x);
            obj.Ts  = 1/obj.Fs;
            obj.t   = (0:(obj.N-1)).*obj.Ts;
        end

        function obj=convertTtoFDomain(obj) %have T domain, want F domain
            if(isempty(obj.f))
                obj.constructFVector();
            end
            obj.X   = fftshift(fft(obj.x));
            obj.PSD = 10.*log10(1./obj.Fs./obj.N.*abs(obj.X).^2);
        end

        function obj =ReadFromFile(obj, fn)
            %time must be in row 1, data in row 2
            M = csvread(fn);
            obj.t = M(:,1);
            obj.x = M(:,2);
            obj.Ts = obj.t(2)-obj.t(1);
            obj.constructFVector();
            obj.convertTtoFDomain();
        end
        
        function ReadFromFile_noTime(obj, fn, Fs)
            obj.x   = csvread(fn);
            obj.Fs  = Fs;
            obj.constructTVector();
            obj.convertTtoFDomain();
        end
        
        % a function to find the bandwidth of a pulse given its PSD
        function [f_low, f_high] = findBandwidth(obj, threshold_dB20)
            % note: if the pulse has ripples of more than threshold_dB20,
            % this will ignore the ripples and give the FIRST and LAST 
            % times that the pulse transitions below threshold_dB20
            if(threshold_dB20 > 0)
                warning('Expected a negative threshold_dB20. Negating input threshold.');
                threshold_dB20 = - threshold_dB20;
            end
            
            PSD_gt0     = obj.PSD(obj.f>0);                             % PSD (already in dBs) where f>0
            f_gt0       = obj.f(obj.f>0);                               % f where f > 0
            PSD_gt0     = (PSD_gt0 - max(PSD_gt0)) >= threshold_dB20;   % now logical - is PSD greater than threshold?
            left_ind    = find(PSD_gt0,1,'first');                      % left index where PSD < threshold
            right_ind   = find(PSD_gt0,1,'last');                       % right index where PSD < threshold
            f_low       = f_gt0(left_ind);                              % low frequency
            f_high      = f_gt0(right_ind);                             % high frequency
        end
        
        function obj = chopPulse(obj, threshold)
            env         = normalize(abs(hilbert(obj.x)));
            [~,ind_max] = max(env);
            env_log     = env > threshold;
            left        = find(env_log(1:ind_max),1,'first');
            right       = find(env_log(ind_max:end),1,'last');
            right       = right + ind_max - 1;

            obj.x = obj.x(left:right);
            obj.t = obj.t(left:right);
            obj.N = numel(obj.t);
            obj.constructFVector();
            obj.convertTtoFDomain();
        end
        
        function vpp = get_vpp(obj)
            vpp = max(obj.x)-min(obj.x);
            obj.vpp = vpp;
        end
        
        function obj = padToGranularity(obj, gran)
            pad = gran - mod(obj.N, gran);
            if(isrow(obj.x))
                zero_pad = zeros(1,pad);
                obj.x = cat(2, obj.x, zero_pad);
                obj.t = cat(2, obj.t, max(obj.t)+(1:pad).*obj.Ts);
            else
                zero_pad = zeros(pad,1);
                obj.x = cat(1, obj.x, zero_pad);
                obj.t = cat(1, obj.t, max(obj.t)+(1:pad)'.*obj.Ts);
            end
            obj.N = numel(obj.x);
            obj.constructFVector();
            obj.convertTtoFDomain();
        end
        
        function obj = delayWRTRef(obj,R)
            % delayWRTRef - delays current pulse with respect to
            % a reference pulse R
            
            % find the peak of R.
            env_ref = abs(hilbert(R.x));
            [~,max_ref_ind] = max(env_ref);
            max_ref_t = R.t(max_ref_ind);
            
            assert(max_ref_t > min(obj.t),...
                "Reference t0 must be greater than the minimum Pulse time");

            % find the peak of obj.
            env_obj = abs(hilbert(obj.x));
            [~,max_obj_ind] = max(env_obj);
            max_obj_t = obj.t(max_obj_ind);
            
            peak_difference = round((max_ref_t - max_obj_t)./obj.Ts);
            
            % how many samples must we shift?
            if(peak_difference > 0) % shift obj.x to the right (pad left)
                if(isrow(obj.x))
                    zero_pad = zeros(1, peak_difference);
                    obj.x = cat(2, zero_pad, obj.x);
                    obj.t = cat(2, obj.t, max(obj.t)+(1:peak_difference).*obj.Ts);
                else
                    zero_pad = zeros(peak_difference,1);
                    obj.x = cat(1, zero_pad, obj.x);
                    obj.t = cat(1, obj.t, max(obj.t)+(1:peak_difference)'.*obj.Ts);
                end    
            else % shift obj.x to the left (eat into the right)
                obj.x = obj.x((peak_difference + 1):end);
                obj.t = obj.t(1:(end-peak_difference));
            end
            obj.N = numel(obj.x);
            obj.constructFVector();
            obj.convertTtoFDomain();
        end
        
        function obj = plot(obj)
            [~,fhi] = obj.findBandwidth(-20);
            figure();
            subplot(211); hold on; grid on;
            xlabel("Time (ns)");
            ylabel("Mag");
            plot(obj.t.*1e9, obj.x,'k');
            title(obj.label + " Time domain");
            
            subplot(212); hold on; grid on;
            xlabel("Freq (GHz)");
            ylabel("PSD (dB/Hz)");
            title(obj.label + " Frequency domain");
            plot(obj.f./1e9, obj.PSD, 'k');
            xlim([0, min(2*fhi,obj.Fs/2)]./1e9);
            yL = get(gca,'YLim');
            ylim([yL(2)-80 yL(2)]);
        end
    end
end