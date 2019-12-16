%% Initialize system constants
clear all,close all,clc;
rng(2014);
%given constants
gc = set2DParameters();

% Tunable parameters
tp.snr = 10;
tp.txGain = 0;           % dB
tp.tgt1_range = 2750;    % m
tp.tgt1_angle = 40;       % degrees
tp.tgt2_range = 2750;    % m
tp.tgt2_angle = 10;       % degrees

tp.intf1_range = 9000;    % m
tp.intf1_angle =   26;    % degrees
tp.intf2_range = 9000;    % m
tp.intf2_angle =   -50;    % degrees

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
%free space loss
rx_tgt1 = cp_sig(:,:,1).* 10^-(fspl(tp.tgt1_range,gc.lambda)/10);
rx_tgt2 = cp_sig(:,:,11).* 10^-(fspl(tp.tgt2_range,gc.lambda)/10);
rx_intf1 = cp_sig(:,:,22).* 10^-(fspl(tp.intf1_range,gc.lambda)/10);
rx_intf2 = cp_sig(:,:,33).* 10^-(fspl(tp.intf2_range,gc.lambda)/10);

rx = collectPlaneWave(ula,[rx_tgt1 rx_tgt2 rx_intf1 rx_intf2],[tp.tgt1_angle tp.tgt2_angle  tp.intf1_angle tp.intf2_angle; 0 0 0 0],gc.fc);
rx = awgn(rx,tp.snr,'measured');

musicEstimator = phased.MUSICEstimator('SensorArray',ula,...
        'OperatingFrequency',gc.fc,'ScanAngles',gc.scan_az,...
        'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',4);

    
[~,doas] = musicEstimator(rx);

%figure
%plotSpectrum(musicEstimator,'NormalizeResponse',true);
disp(doas);
disp([tp.tgt1_angle tp.tgt2_angle  tp.intf1_angle tp.intf2_angle]);

%% Steering
num_tgt = 2;
steer_vec = steer_vec_ula(ula,gc.lambda,doas(1:num_tgt));
y = zeros(size(rx,1),num_tgt);
for j = 1 :size(steer_vec,2)
    s0 = steer_vec(:,j);
    w = s0 / tp.num_tx_el;
    for i = 1:size(rx,1)
        %y(i,j) = w' * reshape(rx(i,:),[tp.num_tx_el 1]);
        y(i,j) = w' * transpose(rx(i,:));
    end
end

% demodulate ofdm symbol
y_cap = zeros(size(data_qam,1),num_tgt);
y_qam = zeros(size(data_qam,1),num_tgt);
for i = 1:size(y_cap,2)
    y_cap(:,i)= ofdm_demod(y(:,i),gc.N_fft,gc.N_cyclepref,gc.pilot_pos,gc.guard_bands);
    y_qam(:,i) = qamdemod(y_cap(:,i),gc.M);
end