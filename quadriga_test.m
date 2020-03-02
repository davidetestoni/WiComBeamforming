clear all 
close all
clc

t = qd_track('linear',200,0);
t.name = 'Terminal';
t.initial_position = [100;100;2];      
t.calc_orientation;
l = qd_layout;   

[~,l.tx_track] = interpolate( t.copy,'distance',0.1 );
l.rx_position = [0;0;15];

l.visualize([],[],0);    
axis equal
title('Track layout');

%%
l.simpar.center_frequency = 26e9; 
l.rx_array = qd_arrayant( '3gpp-3d', 4, 4, l.simpar.center_frequency);
l.rx_track.orientation = [ 0 ; 0 ; 0]*pi/180;
l.rx_name{1} = 'BS'; 

l.tx_array = qd_arrayant('patch');                     
l.tx_array.center_frequency = l.simpar.center_frequency;

%%
c = l.get_channels(0.01);                               % Generate channels

pow  = 10*log10( reshape( sum(abs(c.coeff(:,:,:,:)).^2,3) ,2,[] ) );    % Calculate the power
time = (0:c.no_snap-1)*0.01;   