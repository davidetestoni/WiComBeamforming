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


real_angles = [x_azim rand()*360-180 rand()*360-180 rand()*360-180];
real_angles = [real_angles;x_elev rand()*180-90 rand()*180-90 rand()*180-90];


maxSnr = 20;
snr = 0:maxSnr;
ber = zeros(maxSnr + 1,1);
ber_conv = zeros(maxSnr + 1,1);
snr_out = zeros(maxSnr + 1,1);
%steeringvec = phased.SteeringVector('SensorArray',ura,'PropagationSpeed',physconst('LightSpeed'));

%S = steeringvec(fc,real_angles);
S = steer_vec_ura(ura,lambda,real_angles);
% Null-beamforming
g_1 = [1 0 0 0];
S_inv = pinv(S);%S' / (S * S');
w_h = g_1 * S_inv;

% Conventional beamforming
s_0 = S(:,1);
w_h_conv = (s_0/ura.getNumElements)';

for i = 1 :21
    
    all_sig = [x intf(:,1) intf(:,2) intf(:,3)];
    rx = collectPlaneWave(ura,all_sig,real_angles,fc);    
    rx_n = awgn(rx,snr(i),mean(abs(x).^2));
    
    all_noise_in = all_sig - x;
    
    y =  rx_n * transpose(w_h);
    noise_out = y - x;
    
    y_s = qamdemod(y,M);
    y_b = de2bi(y_s);
    [ ~,ber(i) ] = biterr(x_bit,y_b);
    
    
    gain = 10*log10(mean(mean(abs(rx_n - rx).^2)) / mean(abs(noise_out).^2));
    snr_out(i) = gain + snr(i);
    %gain
    
    
    y_conv = rx_n * transpose(w_h_conv);
    
    n_out_conv = y_conv - x;
    
    gain_conv = 10*log10( mean(mean(abs(rx_n - rx).^2)) / mean(abs(n_out_conv).^2) );
    snr_out_conv(i) = gain_conv + snr(i);
    
    y_conv = de2bi(qamdemod(y_conv,M));
    [ ~,ber_conv(i) ] = biterr(x_bit,y_conv);
end
figure
plot(snr,snr_out,snr,snr_out_conv);
title("Consider SNR input")

legend("Null","Conventional");
xlabel("SINR input");
ylabel("SINR output");


% figure
% semilogy(snr,ber,snr,ber_conv)
% legend("Null","Conventional");