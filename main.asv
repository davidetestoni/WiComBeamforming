%% Initialize system constants
rng(2014);
gc = set2DParameters();

% Tunable parameters
tp.txPower = 9;           % watt
tp.txGain = 0;           % dB
tp.mobileRange = 2750;    % m
tp.mobileAngle = 3;       % degrees
tp.interfPower = 1;       % watt
tp.interfGain = -15;      % dB
tp.interfRange = 9000;    % m
tp.interfAngle =   20;    % degrees
tp.numTXElements = 8;       
tp.steeringAngle = 0;     % degrees
tp.rxGain = 100; % dB ???

numTx= tp.numTXElements;

%% antenna definition

ula = phased.ULA( tp.numTXElements, ...
        'ElementSpacing', 0.5*gc.lambda, ...
        'Element', phased.IsotropicAntennaElement('BackBaffled', true));

steeringvec = ...
    phased.SteeringVector('SensorArray',ula,'PropagationSpeed',gc.cLight);

arrayresp = phased.ArrayResponse('SensorArray',ula,'WeightsInputPort',true);

