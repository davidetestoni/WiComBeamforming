clear all 
close all
clc

t = qd_track('linear',200,0);
t.name = 'Terminal';
t.initial_position(3,1) = 2;      
t.calc_orientation;
l = qd_layout;   

[~,l.rx_track] = interpolate( t.copy,'distance',0.1 );
l.visualize([],[],0);    
axis equal
title('Track layout');