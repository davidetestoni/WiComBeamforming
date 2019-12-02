function gc = set2DParameters
%% all given constants

gc.modOrder = 4;
gc.modMode = gc.modOrder^2;
gc.upSample = 8; %?
gc.cLight = 3e8;
gc.fc = 2.4e9;               % 2.4 GHz ISM Band
gc.lambda = gc.cLight/gc.fc;
gc.nT = 290;       % Noise Temp in deg K
gc.scanAz = -180:180;

sz = get(0, 'ScreenSize');

sz2 = repmat(sz(3:4), 1,2);
gc.envPlotPosition = round([0, 0.25, 0.25, 0.25] .* sz2);
gc.constPlotPosition = round([0.0125, 0.52, 0.2381, 0.38] .* sz2);
