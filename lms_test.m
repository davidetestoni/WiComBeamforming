clear all, close all, clc;
M = 16;
N_bit = 1e4;
N_tx_el = [4 4];
fc = 26e9;               % 26 GHz 
lambda = physconst('LightSpeed')/fc;

x_bit = randi([0 1],N_bit,1);
x_bit = reshape(x_bit,ceil(length(x_bit)/log2(M)),log2(M));
x_sample = bi2de(x_bit);

x = qammod(x_sample,M);

intf_sample = randi([0 M-1],N_bit/log2(M),3);
intf = qammod(intf_sample,M);

x_elev = rand()*180-90;
x_azim = rand()*360-180;


ura = phased.URA( N_tx_el, ...
        'ElementSpacing', 0.5*lambda, ...
        'Element', phased.IsotropicAntennaElement('BackBaffled', false));

all_sig = [x intf(:,1) intf(:,2) intf(:,3)];

real_angles = [x_azim rand()*360-180 rand()*360-180 rand()*360-180];
real_angles = [real_angles;x_elev rand()*180-90 rand()*180-90 rand()*180-90];

rx = collectPlaneWave(ura,all_sig,real_angles,fc);    

rx_n = awgn(rx,10,'measured');



%steeringvec = phased.SteeringVector('SensorArray',ura,'PropagationSpeed',physconst('LightSpeed'));

%S = steeringvec(fc,real_angles);
S = steer_vec_ura(ura,lambda,real_angles);

g_1 = [1 0 0 0];
S_inv = S' / (S * S');

w_h = g_1 * S_inv;

