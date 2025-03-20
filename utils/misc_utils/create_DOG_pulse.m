function P = create_DOG_pulse(derivative_number, sigma, time_vec, delay, amplitude, polarity)

% function to create a Pulse() object with a corresponding DOG 
% (Derivative Of Gaussian) pulse.
%
% inputs: 
%   derivative_number: determines the shape of pulse
%                      0 - gaussian
%                      1 - differential
%                      2 - ricker
%                      3 - 3rd derivative ... etc
%   sigma:      the width parameter (o) of the base gaussian pulse
%   time_vec:   the time vector over which to calculate the pulse
%   delay:      where the envelope of the output pulse will peak in time
%   amplitude:  the peak value of the pulse DEFAULT 1
%   polarity:   if 0, do not flip pulse. if 1, flip pulse DEFAULT 0
%
%
% outputs:
%   P: pulse object with corresponding parameters
%
% example:
%   P = create_DOG_pulse(1, 50e-12, linspace(-5e-9,5e-9,1024), 2e-9, 4, 0)
%
%   will create a pulse object P with a differential sigma=50ps pulse over
%   the time scale -5 to 5 ns delayed by 2 ns with a peak amplitude of 4.

if (~exist('amplitude','var'))
    amplitude = 1;
end

if(~exist('polarity','var'))
    polarity = 0;
end

% create the base gaussian pulse
DoG = exp(-1/2 .* ((time_vec-delay)./sigma).^2);

% take the derivative of DoG derivative_number of times
for ii = 1:derivative_number
    DoG = gradient(DoG);
end

% flip if requested.
if(polarity)
    DoG = -DoG;
end

% instantiate the pulse object
P = Pulse();

% fill out relevant parameters
P.x     = DoG;
P.Ts    = time_vec(2)-time_vec(1);
P.t     = time_vec;
P.N     = numel(time_vec);
P.normalizeTtoVoltage(amplitude);








