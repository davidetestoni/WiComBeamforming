%% Initialize system constants
clear all,close all,clc;
rng(2014);
%given constants
gc = set2DParameters();

% Tunable parameters
tp.txGain = 0;           % dB
tp.mobileRange = 2750;    % m
tp.mobileAngle = 3;       % degrees
tp.interfGain = -15;      % dB
tp.interfRange = 9000;    % m
tp.interfAngle =   20;    % degrees
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

musicEstimator = phased.MUSICEstimator('SensorArray',ula,...
        'OperatingFrequency',gc.fc,'ScanAngles',gc.scan_az,...
        'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',2);
    return
[~,ang] = musicEstimator(signal)
ymvdr = mvdrspatialspect(signal);
ymusic = musicEstimator(signal);
helperPlotDOASpectra(mvdrspatialspect.ScanAngles,...
  musicEstimator.ScanAngles,ymvdr,ymusic,'ULA')