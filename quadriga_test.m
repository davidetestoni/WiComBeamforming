clear all 
close all
clc

%% UE track
t = qd_track('linear',200,0);
t.name = 'Terminal';
t.initial_position = [100;100;2];      
t.calc_orientation;
t.movement_profile = [0, 50; 0, 50];
%t.scenario = '3GPP_38.881_Urban_LOS'; 

%% apparently it is necessary to create a track also for the stationary BS
t2 = qd_track('linear',0,0);
t2.name = 'Base';
t2.initial_position = [0;0;15];      
t2.calc_orientation;
t2.scenario = '3GPP_38.881_Urban_LOS'; 

%% quadriga layout
l = qd_layout;   
l.tx_track = t;
l.rx_track = t2;
l.simpar.center_frequency = 26e9; 
l.rx_array = qd_arrayant( '3gpp-3d', 4, 4, l.simpar.center_frequency); % 4x4 array BS
l.rx_name{1} = 'BS';

l.tx_array = qd_arrayant('patch'); % UE antenna is a patch                  
l.tx_array.center_frequency = l.simpar.center_frequency;

l.visualize([],[],0);
axis equal
title('Track layout');



%%
c = l.get_channels(0.01);                               % Generate channels
% it took 10 minutes to generate the channels (200s observation time)
% 2 minutes with a 50s observation time
% in both cases generated 34 paths

pow  = 10*log10( reshape( sum(abs(c.coeff(:,:,:,:)).^2,3) ,2,[] ) );    % Calculate the power
time = (0:c.no_snap-1)*0.01;   