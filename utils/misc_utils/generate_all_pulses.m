addpath(genpath('.\utils\'));
fixit;
prettygraphs;
DEFS;

ref_fn = "_delayed_Differential_o_100ps.csv";
Fs = 64e9;
output_fn = INPUT_PULSE_DIR + "base_ricker_pulses\";

% get reference pulse
R = Pulse();
R.ReadFromFile_noTime(INPUT_PULSE_DIR + ref_fn, Fs);
t = R.t;

o_v = 20e-12:5e-12:300e-12;
N_o = numel(o_v);

for ii = 1:N_o
    P = create_DOG_pulse(0, o_v(ii), t, 3e-9, 1, 0);
    P.delayWRTRef(R);
    P.x = P.x(1:numel(R.x));
    csvwrite(output_fn + "ricker_o_"+num2str(o_v(ii)*1e12)+"ps_rxv2.csv", P.x);
    
end

%%


