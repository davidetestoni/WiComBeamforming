%% Initialize system constants
clear all,close all,clc;
rng(2014);
%given constants
gc = set2DParameters();

% Tunable parameters
tp.snr = 10;
tp.txGain = 0;           % dB
tp.tgt1_range = 2750;    % m
tp.tgt1_angle = -50;       % degrees
tp.tgt2_range = 2750;    % m
tp.tgt2_angle = 30;       % degrees

tp.intf1_range = 9000;    % m
tp.intf1_angle =   50;    % degrees
tp.intf2_range = 9000;    % m
tp.intf2_angle =   29;    % degrees

tp.numTXElements = 8;       
tp.steeringAngle = 0;     % degrees
tp.rxGain = 100; % dB ???

tp.N_bits = 1e4;
gc.N_blocks = ceil(tp.N_bits/gc.N_bit_block);

numTx= tp.numTXElements;

%% antenna definition

ula = phased.ULA( tp.numTXElements, ...
        'ElementSpacing', 0.5*gc.lambda, ...
        'Element', phased.IsotropicAntennaElement('BackBaffled', true));

arrayresp = phased.ArrayResponse('SensorArray',ula,'WeightsInputPort',true);

visualize2D(gc, tp,ula);

%% Signal generation
[data_bit,data_sample,data_qam,cp_sig] = modulateData(gc,tp);


%% DOA estimation
% target and interferer send signaling data to allow discovery
%free space loss
rx_tgt1 = cp_sig(:,:,1).* 10^-(fspl(tp.tgt1_range,gc.lambda)/10);
rx_tgt2 = cp_sig(:,:,2).* 10^-(fspl(tp.tgt2_range,gc.lambda)/10);
rx_intf1 = cp_sig(:,:,3).* 10^-(fspl(tp.intf1_range,gc.lambda)/10);
rx_intf2 = cp_sig(:,:,4).* 10^-(fspl(tp.intf2_range,gc.lambda)/10);

rx = collectPlaneWave(ula,[rx_tgt1 rx_tgt2 rx_intf1 rx_intf2],[tp.tgt1_angle tp.tgt2_angle  tp.intf1_angle tp.intf2_angle; 0 0 0 0],gc.fc);
rx = awgn(rx,tp.snr,'measured');

musicEstimator = phased.MUSICEstimator('SensorArray',ula,...
        'OperatingFrequency',gc.fc,'ScanAngles',gc.scan_az,...
        'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',4);

    
[y,doas] = musicEstimator(rx);
figure
plotSpectrum(musicEstimator,'NormalizeResponse',true);