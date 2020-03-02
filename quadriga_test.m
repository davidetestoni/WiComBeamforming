clear all 
close all
clc

t = qd_track('linear',200,0);
t.name = 'Terminal';
t.initial_position = [100;100;2];      
t.calc_orientation;
t.movement_profile = [0,100;0,200];
l = qd_layout;   

%[~,l.tx_track] = interpolate( t.copy,'distance',0.1 );
l.tx_track = t;

t = qd_track('linear',0,0);
t.initial_position = [0;0;15];
l.rx_track = t;

l.visualize([],[],0);    
axis equal
title('Track layout');

%%
l.simpar.center_frequency = 26e9; 
l.rx_array = qd_arrayant( '3gpp-3d', 4, 4, l.simpar.center_frequency);
l.set_scenario( '3GPP_38.881_Urban_LOS' );


l.rx_name{1} = 'BS'; 
l.tx_name{1} = 'MS';
l.tx_array = qd_arrayant('patch');                     
l.tx_array.center_frequency = l.simpar.center_frequency;

%%
c = l.get_channels();                               % Generate channels

pow  = 10*log10( reshape( sum(abs(c.coeff(:,:,:,:)).^2,3) ,2,[] ) );    % Calculate the power
time = (0:c.no_snap-1)*0.01;   