clear all;
close all;
clc

fc = 1e9;
cLight = 3e8;
lambda = cLight/fc;
angle = 55;
txEl = 8;
k = 2*pi / lambda;
d = 0.5*lambda;

ula = phased.ULA( txEl, ...
        'ElementSpacing', d, ...
        'Element', phased.IsotropicAntennaElement('BackBaffled', false));

steeringvec = ...
    phased.SteeringVector('SensorArray',ula,'PropagationSpeed',cLight);

s = zeros(txEl,1);
for n = 1:txEl
    s(txEl - n+1) = exp(-1i*k*(n-1) * d * sin(degtorad(angle)) );
end

arrayresp = phased.ArrayResponse('SensorArray',ula,'WeightsInputPort',true);


% Steer the transmitter main lobe
wT = steeringvec(fc,[angle;0]);

% Get array response for visualization
arrayResponse = arrayresp(fc,-180:180,wT);
aR = arrayresp(fc,-180:180,s);
polarplot(degtorad(-180:180),abs(arrayResponse(:))*5e3/4, '-m');
figure
polarplot(degtorad(-180:180),abs(aR(:))*5e3/4, '-r');



