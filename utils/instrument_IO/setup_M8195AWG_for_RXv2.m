function awg = setup_M8195AWG_for_RXv2(IP, set_reset, voltset)        
    awg = visa('agilent', 'TCPIP0::localhost::hislip0::INSTR');
    awg.OutputBufferSize = 2^20;
    awg.InputBufferSize = 2^20;
    
    fopen(awg); %open AWG
    fprintf(awg,":ABOR");
    if(set_reset)
        fprintf(awg,"*RST"); %reset settings to default
    end
    fprintf(awg,":ABOR");

    % clock tab settings
    fprintf(awg,":SOUR:ROSC:RANG RANG1");       % range 10 to 300 MHz (range1)
    fprintf(awg,":SOUR:ROSC:SOUR EXT");         % source: external
    fprintf(awg,":OUTP:ROSC:SOUR EXT");         % output ref clk is external
    fprintf(awg,":SOUR:ROSC:FREQ 10000000");    % set freq to 10 MHz

    % output tab settings
    fprintf(awg, ":VOLT1 "+num2str(voltset));                   % sets voltage on channel 1
    fprintf(awg,":SOUR:FREQ:RAST 64000000000"); 

    % trigger tab settings
    fprintf(awg, ":ARM:TRIG:OPER SYNC");        % synchronous operation
    fprintf(awg, ":ARM:TRIG:SOUR TRIG");        % source "Trigger" from "Trig in"
    fprintf(awg, ":ARM:TRIG:SLOP POS");         % Trigger on rising edge
    fprintf(awg, ":TRIG:SOUR:ADV INT");         % source "Advance" from "Event In"
    fprintf(awg, ":TRIG:SOUR:ENAB EVEN");       % source "Enable" from "Trig In"
    fprintf(awg, ":ARM:EVEN:SLOP POS");         % set level event to 250mV

    fprintf(awg, ":ARM:TRIG:LEV 1.0");         % set level trigger to 250mV
    fprintf(awg, ":ARM:EVEN:LEV 1.5");         % set level event to 250mV

    fprintf(awg, ":INIT:CONT:ENAB SELF");        % Armed mode
    fprintf(awg, ":INIT:CONT:STAT OFF");

    % import the first file into TRAC1
    fprintf(awg,":TRACe1:DEL:ALL");
    fprintf(awg,":INST:DACM DUAL");
    fprintf(awg,":INST:MEM:EXT:RDIV DIV1");
    fprintf(awg,":TRACe1:MMODe EXTended");

    fprintf(awg, ":TRACe1:DEF 1, "     + num2str(numel(IP.x))); 
    fprintf(awg, ":TRACe1:DATA 1, 0, " + join(string(IP.x),', '));
    fprintf(awg, ":OUTP1:FILT:FRAT:SCAL 1"); 

    fprintf(awg,":FUNCtion:MODE ARBitrary");

    % turn on the output
    fprintf(awg, ":OUTP1 ON");                  
    fprintf(awg, ":TRAC:SEL 1");
    fprintf(awg, ":INIT:IMM"); 
end