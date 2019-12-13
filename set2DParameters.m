function gc = set2DParameters
%% all given constants

gc.up_sample = 8; %?
gc.c_light = 3e8;
gc.fc = 2.4e9;               % 2.4 GHz ISM Band
gc.lambda = gc.c_light/gc.fc;
gc.nT = 290;       % Noise Temp in deg K
gc.scan_az = -90:90;

% modulation
gc.M = 16;
gc.N_fft = 64;
gc.pilot_pos = [12;26;40;54];
gc.pil = 3+3j;
gc.guard_bands = [6,5];
gc.N_bit_block = (gc.N_fft - max(size(gc.pilot_pos)) - sum(gc.guard_bands)) * log2(gc.M); %block is 1 ofdm symbol
gc.N_cyclepref = 16;
gc.N_carried_data = gc.N_fft - numel(gc.pilot_pos) - sum(gc.guard_bands);
