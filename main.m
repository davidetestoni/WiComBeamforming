%% Initialize system constants
clear all,close all,clc;
rng(2014);
%given constants
gc = set2DParameters();

% Tunable parameters
tp.snr = 10;

tp.tgt(1).range = 2750;    % m
tp.tgt(1).angle = 0;       % degrees
tp.tgt(2).range = 2750;    % m
tp.tgt(2).angle = 10;       % degrees

tp.intf(1).range = 9000;    % m
tp.intf(1).angle =   26;    % degrees
tp.intf(2).range = 9000;    % m
tp.intf(2).angle =   -50;    % degrees

tp.num_tx_el = 16;       
tp.steer_angle = 0;     % degrees

tp.N_bits = 1e4;
gc.N_blocks = ceil(tp.N_bits/gc.N_bit_block);


%% antenna definition

ula = phased.ULA( tp.num_tx_el, ...
        'ElementSpacing', 0.5*gc.lambda, ...
        'Element', phased.IsotropicAntennaElement('BackBaffled', true));

arrayresp = phased.ArrayResponse('SensorArray',ula,'WeightsInputPort',true);

%visualize2D(gc, tp,ula);

%% Signal generation
[data_bit,data_sample,data_qam,cp_sig] = modulateData(gc,tp);


%% DOA estimation
% target and interferer send signaling data to allow discovery

tp.tgt(1).tx_sig = cp_sig(:,:,1);
tp.tgt(1).data_bit = data_bit(:,:,1);


tp.tgt(2).tx_sig = cp_sig(:,:,11);
tp.tgt(2).data_bit = data_bit(:,:,11);

tp.intf(1).tx_sig = cp_sig(:,:,22);
tp.intf(1).data_bit = data_bit(:,:,22);

tp.intf(2).tx_sig = cp_sig(:,:,33);
tp.intf(2).data_bit = data_bit(:,:,33);

%free space loss
tp.tgt(1).rx_sig = tp.tgt(1).tx_sig.* 10^-(fspl(tp.tgt(1).range,gc.lambda)/10);
tp.tgt(2).rx_sig = tp.tgt(2).tx_sig.* 10^-(fspl(tp.tgt(2).range,gc.lambda)/10);
tp.intf(1).rx_sig = tp.intf(1).tx_sig.* 10^-(fspl(tp.intf(1).range,gc.lambda)/10);
tp.intf(2).rx_sig = tp.intf(2).tx_sig.* 10^-(fspl(tp.intf(2).range,gc.lambda)/10);

all_sig = [tp.tgt(1).rx_sig tp.tgt(2).rx_sig tp.intf(1).rx_sig tp.intf(2).rx_sig];
rx = collectPlaneWave(ula,all_sig,[tp.tgt(1).angle tp.tgt(2).angle  tp.intf(1).angle tp.intf(2).angle; 0 0 0 0],gc.fc);
%rx = collectPlaneWave(ula,[tp.tgt(1).rx_sig tp.tgt(2).rx_sig],[tp.tgt(1).angle tp.tgt(2).angle;0 0],gc.fc);

rx = awgn(rx,tp.snr,'measured');

musicEstimator = phased.MUSICEstimator('SensorArray',ula,...
        'OperatingFrequency',gc.fc,'ScanAngles',gc.scan_az,...
        'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',4);

    
[~,doas] = musicEstimator(rx);

%figure
%plotSpectrum(musicEstimator,'NormalizeResponse',true);
disp(doas);
disp([tp.tgt(1).angle tp.tgt(2).angle  tp.intf(1).angle tp.intf(2).angle]);

%% Steering
num_tgt = length(tp.tgt);
steer_vec = steer_vec_ula(ula,gc.lambda,doas(1:num_tgt));
% weighted response
y = zeros(size(rx,1),num_tgt);
for j = 1 :size(steer_vec,2)
    s0 = steer_vec(:,j);
    w = s0 / tp.num_tx_el;
    for i = 1:size(rx,1)
        y(i,j) = w' * transpose(rx(i,:));
    end
end
%% demodulate ofdm symbol
y_cap = zeros(size(data_qam,1),num_tgt);
y_sample = zeros(size(data_qam,1),num_tgt);
for i = 1:size(y_cap,2)
    y_cap(:,i) = ofdm_demod(y(:,i),gc.N_fft,gc.N_cyclepref,gc.pilot_pos,gc.guard_bands);
    y_sample(:,i) = qamdemod(y_cap(:,i),gc.M);
end

[~,ber(1)] = biterr(de2bi(y_sample(:,1)),tp.tgt(1).data_bit);
[~,ber(2)] = biterr(de2bi(y_sample(:,2)),tp.tgt(2).data_bit);
[~,ber_c(1)] = biterr(de2bi(y_sample(:,2)),tp.tgt(1).data_bit);
[~,ber_c(2)] = biterr(de2bi(y_sample(:,1)),tp.tgt(2).data_bit);
if ber(1) < ber_c(1)
    clear ber_c
else
    ber = ber_c;
    clear ber_c
end
