% Target2D Class - a class containing a 2-Dimensional B-Scan target
% meant to be populated by a function similar to find_2D_targets.m
%
% Samuel Wagner
% UCD ECE/MML Jan 21,2022
% rev1
classdef Target2D < handle
    properties
        % horizontal (x) properties
        upper_left_x    (1,1) single; % x-position (m) of lower left corner
        x_centroid      (1,1) single; % x position (m) of target centroid
        x_width         (1,1) single; % width of target rectangle (x, m)
        x               (:,1) single; % the x vector (in meters)
        Nx              (:,1) single; % the number of x elements
        dx              (1,1) single; % the spacing in the x-direction (should be = pixel_width)

        % vertical (z) properties
        upper_left_z    (1,1) single; % z-position (m) of lower left corner
        z_centroid      (1,1) single; % z position (m) of target centroid
        z_height        (1,1) single; % height of target rectangle (z, m)
        dz              (1,1) single; % the spacing in the z-directoin (should be = pixel_width)
        z               (:,1) single; % the z vector (in meters)
        Nz              (:,1) single; % the number of z elements
        
        % target or derived properties
        area            (1,1) single; % x_width * z_height (m^2)
        Image           (:,:) single; % target rectangle, un-normalized.
        SNR             (1,1) single; % signal-to-noise ratio of this target
        peak            (1,1) single; % peak power of the target
        energy          (1,1) single; % total energy of this target
        WL_ratio        (1,1) single; % the width/length ratio of the target
    end
    
    methods
        % Target2D methods
        % typically, the only method that you want to really call is
        % T.fill_target(). The constructor will call that if correctly
        % supplied.
        %
        % T = Target2D(S, I, upper_left_x, upper_left_z, x_width, z_height)
        %   | constructor method. will call T.fill_target() if possible.
        % T.fill_target(S, I,upper_left_x, upper_left_z, x_width, z_height)
        %   | main population method. calculates and fills most properties in the class 
        % T.calculate_centroid()
        %   | calculates the centroid based on Image and the x,z axes.
        %     called by fill_target.
        % T.draw(ax, scalex, scalez)
        %   | will draw a white target onto the axes given by ax
        %     with the x- and z-scales scalex and scalez
        %     for example, if you want to plot with a x-scale in meters
        %     and a z-scale in cm, call: T.draw(ax, 1, 1e2);
        
        
        % Constructor
        function obj = Target2D(S, I, upper_left_x, upper_left_z, x_width, z_height)
            % inputs __________
            % S            - Scan object from which the target originates
            % I            - Image (2D double/single) from which the target originiates
            % upper_left_x - left x-value of target start (m)
            % upper_left_z - top (near-ground) z-value of target start (m)
            % x_width      - how long / wide the target is (m)
            % z_height     - how tall / deep the target is (m)
            
            if(nargin == 0)
                obj.upper_left_x = 0;
                obj.upper_left_z = 0;
                obj.x_width      = 0;
                obj.z_height     = 0;
            else
                % input checking
                % anonymous functions for abstraction
                is_double_single_scalar = @(x) (strcmp(string(class(x)),"double")||strcmp(string(class(x)),"single")) && numel(x) == 1;
                is_double_single_vector = @(x) (strcmp(string(class(x)),"double")||strcmp(string(class(x)),"single")) && numel(x) >= 1;
                is_double_single_matrix  = @(x) (strcmp(string(class(x)),"double")||strcmp(string(class(x)),"single")) && numel(x) ~= max(size(x)) && numel(size(x)) == 2;
                is_one_scan             = @(x)  strcmp(string(class(x)),"Scan")  && (numel(x) == 1);
                isaninteger             = @(x) (isfinite(x) & x==floor(x));    
                                
                % assertions - throw an error if not true
                assert(is_one_scan(S),                        "S must be a 1x1 Scan object.");
                assert(is_double_single_vector(S.xC),         "S.xC (imaging x-axis) must be filled.");
                assert(is_double_single_vector(S.zC),         "S.zC (imaging z-axis) must be filled.");
                assert(is_double_single_matrix(I),            "I must be a 2-dimensional double or single matrix.");  
                assert(is_double_single_scalar(upper_left_x), "upper_left_x must be a scalar double or single.");
                assert(is_double_single_scalar(upper_left_z), "upper_left_z must be a scalar double or single.");  
                assert(is_double_single_scalar(x_width),      "x_width must be a scalar double or single.");                      
                assert(is_double_single_scalar(z_height),     "z_height must be a scalar double or single.");  
                
                if(sum(isaninteger([x_width, z_height, upper_left_x, upper_left_z])) > 0)
                    warning("Target2D expects x_width, z_height, upper_left_x, and upper_left_z to be in units of meters. "+...
                        "The inputs may have been in units of samples. Double check the inputs");
                end
                
                obj.fill_target(S, I, upper_left_x, upper_left_z, x_width, z_height);               
            end
        end
        
        function obj = fill_target(obj,S, I,upper_left_x, upper_left_z, x_width, z_height)
            % at this point, inputs have been checked. all input defs are
            % the same as in the constructor.

            % fill x parameters
            obj.upper_left_x = upper_left_x;
            obj.x_width      = x_width;
            obj.dx           = S.xC(2)-S.xC(1);
            obj.x            = obj.upper_left_x:obj.dx:(obj.upper_left_x + obj.x_width);
            obj.Nx           = numel(obj.x);
            
            % fill z parameters
            obj.upper_left_z = upper_left_z;
            obj.z_height     = z_height;
            obj.dz           = S.zC(2)-S.zC(1);
            obj.z            = obj.upper_left_z:obj.dz:(obj.upper_left_z + obj.z_height);
            obj.Nz           = numel(obj.z);
            
            % obtain target area, W/L ratio
            obj.area         = obj.x_width * obj.z_height;
            obj.WL_ratio     = obj.x_width/obj.z_height;
            
            % extract the image and related parameters
            [~,xl]  = min(abs(S.xC-obj.upper_left_x));      
            [~,zl]  = min(abs(S.zC-obj.upper_left_z));
            xu      = xl + obj.Nx - 1;
            zu      = zl + obj.Nz - 1;
            
            assert(xu <= size(I,2) && zu <= size(I,1), "target upper limits exceed image dimensions");
            
            obj.Image       = I(zl:zu,xl:xu);
            obj.energy      = trapz(obj.dz,trapz(obj.dx,obj.Image));
            obj.peak        = max(obj.Image(:));
            
            % calculate the energy centroid (x,z)
            obj.calculate_centroid();
        end
        
        function obj = calculate_centroid(obj)
            I_renorm = normalize(obj.Image);
            M = sum(I_renorm(:));
            obj.x_centroid = sum(sum(I_renorm.*repmat(obj.x', obj.Nz, 1)))./M;
            obj.z_centroid = sum(sum(I_renorm.*repmat(obj.z, 1, obj.Nx)))./M;
        end
        
        function draw(obj,ax,scalex,scalez)
           if(~exist('scalex','var'))
              scalex=1; 
           end
           if(~exist('scalez','var'))
              scalez=1; 
           end
           text_z = 0.012;
           axes(ax);
           plot(obj.upper_left_x .*[1 1].*scalex,obj.upper_left_z.*scalez+[0 obj.z_height].*scalez,'w-');
           plot((obj.upper_left_x + obj.x_width) .*[1 1].*scalex, obj.upper_left_z.*scalez+[0 obj.z_height].*scalez,'w-');
           plot(obj.upper_left_x.*scalex+[0 obj.x_width].*scalex, obj.upper_left_z .*[1 1].*scalez,'w-');
           plot(obj.upper_left_x.*scalex+[0 obj.x_width].*scalex, (obj.upper_left_z + obj.z_height) .*[1 1].*scalez,'w-');
           plot(obj.x_centroid.*scalex, obj.z_centroid.*scalez, 'wx');
           text((obj.upper_left_x + obj.x_width).*scalex, obj.upper_left_z.*scalez, "SNR="+sprintf("%.1f",obj.SNR)+"dB",'Color','white');
           text((obj.upper_left_x + obj.x_width).*scalex, (obj.upper_left_z + text_z).*scalez, "width="+sprintf("%.1f",obj.x_width*1e2)+"cm",'Color','white');
           text((obj.upper_left_x + obj.x_width).*scalex, (obj.upper_left_z + 2*text_z).*scalez, "height="+sprintf("%.1f",obj.z_height*1e2)+"cm",'Color','white');
           text((obj.upper_left_x + obj.x_width).*scalex, (obj.upper_left_z + 3*text_z).*scalez, "(x,z)=("+sprintf("%.2f",obj.x_centroid)+","+sprintf("%.2f",obj.z_centroid)+")",'Color','white');
        end
        
    end
end

