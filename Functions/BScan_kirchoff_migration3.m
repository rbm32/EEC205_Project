function I = BScan_kirchoff_migration3(C, h, er, t, x, tx_pos, rx_pos, xQ, zQ, parallellize, progress)
% Sam W. UCD MML
% This is a function that will convert an input B-scan (GPR) into
% an "image" using a Kirchhoff-Migration inspired algorithm
%
% The new, updated & speedy algorithm has 2 parts:
%   part 1 - pre-calculate all indices for what hyperbolas would look like
%            for an antenna centered at x = 0
%   part 2 - iterate through xQ, fit the pre-calculated indices to the
%            data, and integrate to get an image pixel.
%
% INPUTS
% C         - pre-processed B-scan. should be Nt * Nx * 1
% h         - antenna height off of ground
% er        - relative permittivity of ground
% t         - time vector (sampling rate)
% x         - spatial sample vector of array - where each pulse is tx/rx
% tx_pos    - the x/y position of each tx, 1x1  with  row1=x
% rx_pos    - the x/y position of each rx, 1x1  with  row1=x
% xQ        - calculation domain for x (xQ = z QUERY)
% zQ        - calculation domain for z (zQ = z QUERY)
% progress - display progress bar? yes or no (1 or 0)
%
% OUTPUTS
% I - migrated data for [xQ zQ] calculation domain

c   = 3e8;          % speed of light
NxQ = numel(xQ);    % number of x query points
NzQ = numel(zQ);    % number of z query points

% parameters to make sure we don't take ALL x points into integration when
% calculating. - no need, just look at typ. +/- 0.5 m ahead/behind of x pos
x_lims = 0.5;                   % how much +/- in scan x of the tx/rx midpoint to include in calculation.
dx     = x(2)-x(1);             % distance travelled by array inbetween waveforms
xc     = -x_lims:dx:x_lims;     % the true x "calculation" vector - limited to +/- x_lims
Nxc    = numel(xc);             % size of xc


if(mod(Nxc,2)==1)
    Nxc_was_odd_flag = 1;
    xc(end + 1) = x_lims + dx;
    Nxc = Nxc + 1;
else
    Nxc_was_odd_flag = 0;
end

% !Important!
% TTD_inds - True Time Delay indices. Each row is an expected "hyperbola"
% at a given depth. Part 1 of algorithm is tasked with filling this matrix.
% The expected "hyperbola" is traced out by the time-of-arrival of what a point target
% at a given (xq, zq) would produce.
% TTD_inds is calculated by positioning % an antenna pair at x = 0 
% and z = -h and calculating the time-of-arrival % for all of xc and zQ.
% With this matrix, we can drastically speed up the migration algorithm.

TTD_inds = zeros(Nxc,NzQ,'uint32');


% iterate over xc (x calculation domain in antenna scan-x units, limited in
% extend)
for ii = 1:Nxc
    % select xq
    xc_q = xc(ii);  

        for jj = 1:NzQ
            % select zq
            zQ_q = zQ(jj);

            % calculate the GPR refraction angles for the Tx antenna
            % GPR_transmission_angles4 is an optimized method to calculate
            % the angles by minimizing time-of-arrival between two points
            % in a dielectric half-space
            [theta_a_tx, theta_g_tx] = GPR_transmission_angles4(...
                er,...              % dielectric permittivity
                h,...               % height of antenna (z-value)
                tx_pos + xc_q,...   % x-value of antenna
                0,...               % y-value of antenna
                0,...               % query x-point
                0,...               % query y-point
                zQ_q);              % query z-point

            % True Time Delay due to TX - line / ray approximation
            TTD_tx = (h*sec(theta_a_tx) + zQ_q*sec(theta_g_tx)*sqrt(er))/c; 

            % calculate the GPR angles for the Rx - can use the same
            % exact method because source and receiver are interchangeable
            [theta_a_rx, theta_g_rx] = GPR_transmission_angles4(...
                er,...              % dielectric permittivity
                h,...               % height of antenna (z-value)
                rx_pos + xc_q,...   % x-value of antenna
                0,...               % y-value of antenna
                0,...               % query x-point
                0,...               % query y-point
                zQ_q);              % query z-point

            % True Time Delay (seconds) due to RX and TX
            TTD = (h*sec(theta_a_rx) + zQ_q*sec(theta_g_rx)*sqrt(er))/c + TTD_tx; 
            
            % Convert the TTD (seconds) into time indices (1 ... Nt)
            % min1 is a binary search method, as t is a sorted time vector
            % it will return the time index at which TTD occurs
            TTD_inds(ii,jj) = min1(t,TTD);                    
        end         %z-iter
end                 %x-iter

% !Important!
% Part 2 of the algorithm - sliding TTD_inds over the data window
% to recreate migration integration
% 
% Iterate over xQ and zQ - the query pixel points at which to calculate
% the image magnitude
% 
% The majority of this loop is to calculate where the TTD_inds window
% lines up and the extents of its applicability

% Empty image matrix NzQ by NxQ - rows are depth, columns are x-extent
I = zeros(NzQ, NxQ, 'single');

%iterate over xQ
if(parallellize)
parfor ii = 1:NxQ
     % goal: find the columns in the data matrix C that will be
    %       integrated over. and align with TTD_inds.
    
    xq          = xQ(ii);           % what xQ are we looking at? xq
    [~,x0_ind]  = min(abs(x-xq));   % find the index of xq in antenna-scan-x units
    
    %determine the limits of antenna-scan-x that we will look at.
    xa_lower = max(min(x),x(x0_ind)-x_lims);    % minimum x to consider in integration
    xa_upper = min(max(x),x(x0_ind)+x_lims + Nxc_was_odd_flag*dx);    % maximum x to consider in integration

    % find out what indices in x to subset
    size_left  = floor((x(x0_ind)-xa_lower) ./ dx); % on the left of xq
    size_right = floor((xa_upper-x(x0_ind)) ./ dx); % on the right of xq
 
    % important! center xq onto the "middle" (Nxc/2, x = 0) of TTD_inds
    % so that the pre-calculation (algo. step 1) is valid
    ind_left_TTD    = floor(Nxc/2)-size_left;   % left TTD_ind index to consider
    ind_right_TTD   = floor(Nxc/2)+size_right;  % right TTD_ind index to consider
    
    % how much do we shift from x = 0?
    [~,xa_lower_ind]=min(abs(x-xa_lower));
    
    % subset the data matrix to select the correct data columns
    % after this step, C11 contains ONLY columns that will be integrated
    % over. By performing all of the above steps in this loop, we
    % avoid integrating over the entire matrix. much faster.
    C11 = C(:,xa_lower_ind + (0:(size_left+size_right)));
    
    % create a dummy vector based on C11's size to help with sub2ind later
    I2  = (1:size(C11,2))';
    
    % iterate over zQ
    for jj = 1:NzQ
        % select a TTD_vector (called TTD_mat from C-scan migration)
        % that corresponds to the zq. also select the columns derived
        % from the size_left and size_right
        TTD_mat     = TTD_inds(ind_left_TTD:ind_right_TTD,jj);
        
        % convert TTD_mat into linear indices (sub2ind) so we can subset
        % the data matrix in a single shot - no loops! loops, slow.
        linInds     = sub2ind(size(C11),TTD_mat(:),I2);
        
        % integrate the data matrix, subsetted by linInds. linInds
        % describes a hyperbola at (xq,zq).
        I(jj,ii)    = sum(C11(linInds)); % integration is basically summation, right?
    end
end
else % no parallellization. sad!
if(progress)
    progressbar('Sliding B-scan migration')
end
for ii = 1:NxQ
     % goal: find the columns in the data matrix C that will be
    %       integrated over. and align with TTD_inds.
    
    xq          = xQ(ii);           % what xQ are we looking at? xq
    [~,x0_ind]  = min(abs(x-xq));   % find the index of xq in antenna-scan-x units
    
    %determine the limits of antenna-scan-x that we will look at.
    xa_lower = max(min(x),x(x0_ind)-x_lims);    % minimum x to consider in integration
    xa_upper = min(max(x),x(x0_ind)+x_lims + Nxc_was_odd_flag*dx);    % maximum x to consider in integration

    % find out what indices in x to subset
    size_left  = floor((x(x0_ind)-xa_lower) ./ dx); % on the left of xq
    size_right = floor((xa_upper-x(x0_ind)) ./ dx); % on the right of xq
 
    % important! center xq onto the "middle" (Nxc/2, x = 0) of TTD_inds
    % so that the pre-calculation (algo. step 1) is valid
    ind_left_TTD    = floor(Nxc/2)-size_left;   % left TTD_ind index to consider
    ind_right_TTD   = floor(Nxc/2)+size_right;  % right TTD_ind index to consider
    
    % how much do we shift from x = 0?
    [~,xa_lower_ind]=min(abs(x-xa_lower));
    
    % subset the data matrix to select the correct data columns
    % after this step, C11 contains ONLY columns that will be integrated
    % over. By performing all of the above steps in this loop, we
    % avoid integrating over the entire matrix. much faster.
    C11 = C(:,xa_lower_ind + (0:(size_left+size_right)));
    
    % create a dummy vector based on C11's size to help with sub2ind later
    I2  = (1:size(C11,2))';
    
    % iterate over zQ
    for jj = 1:NzQ
        % select a TTD_vector (called TTD_mat from C-scan migration)
        % that corresponds to the zq. also select the columns derived
        % from the size_left and size_right
        TTD_mat     = TTD_inds(ind_left_TTD:ind_right_TTD,jj);
        
        % convert TTD_mat into linear indices (sub2ind) so we can subset
        % the data matrix in a single shot - no loops! loops, slow.
        linInds     = sub2ind(size(C11),TTD_mat(:),I2);
        
        % integrate the data matrix, subsetted by linInds. linInds
        % describes a hyperbola at (xq,zq).
        I(jj,ii)    = sum(C11(linInds)); % integration is basically summation, right?
    end    
    progressbar(ii/NxQ);
end    
    
end

I = abs(I); % take the absolute value.