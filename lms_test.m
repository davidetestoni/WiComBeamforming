clear all, close all, clc;
M = 16;
N_bit = 1e4;
N_tx_el = [4 4];
fc = 26e9;               % 26 GHz 
lambda = physconst('LightSpeed')/fc;

x_bit = randi([0 1],N_bit,1);
x_bit = reshape(x_bit,ceil(length(x_bit)/log2(M)),log2(M));
x_sample = bi2de(x_bit);

N_sample = length(x_sample);

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
ber_opt = zeros(maxSnr + 1,1);
ber_conv = zeros(maxSnr + 1,1);
ber_null = zeros(maxSnr + 1,1);
%steeringvec = phased.SteeringVector('SensorArray',ura,'PropagationSpeed',physconst('LightSpeed'));

%S = steeringvec(fc,real_angles);
S = steer_vec_ura(ura,lambda,real_angles);
% Null-beamforming
g_1 = [1 0 0 0];
S_inv = pinv(S);
w_h = g_1 * S_inv;

% Conventional beamforming
s_0 = S(:,1);
w_h_conv = (s_0/ura.getNumElements)';

Ru = zeros(prod(N_tx_el),prod(N_tx_el),N_sample);


for i = 1 :21
    
    all_sig = [x intf(:,1) intf(:,2) intf(:,3)];
    rx = collectPlaneWave(ura,all_sig,real_angles,fc);    
    rx_n = awgn(rx,snr(i),mean(abs(x).^2));
    
    all_noise_in = all_sig - x;
    
    
    %% MVDR

    % Autocorrelation of the signals rapresent the correlation between all
    % the antennas

    Ru = transpose(rx_n) * transpose(rx_n)'./N_sample;

    w_mvdr = inv(Ru) * s_0 / (s_0' * inv(Ru) * s_0);
    %w_mvdr = w_mvdr / sum(abs(w_mvdr));

    y_mvdr(:,1) = w_mvdr' * transpose(rx_n);
    
    U = transpose(all_sig) * transpose(all_sig)';
    n_pow = mean(mean(abs(rx_n - rx).^2));
    
    Ru_ = S * U * S' + n_pow * eye(prod(N_tx_el)) /N_sample; % alternative formulation. Ru == Ru_
    
    noise_out = y_mvdr - x;
    
    gain = 10*log10( mean(mean(abs(rx_n - rx).^2)) / mean(abs(noise_out).^2));
    snr_mvdr(i) = gain + snr(i);
    
    [ ~,ber_mvdr(i) ] = biterr(x_bit,de2bi(qamdemod(y_mvdr,M)) );
    
    %% Wiener and MSE
    dn = transpose(x);
    un = transpose(rx_n);
    p = (un * dn') /N_sample; % cross correlation beteween rx data and 
    mse = (dn * dn')/N_sample + w_mvdr' * Ru * w_mvdr - 2* w_mvdr' * p;
    w_opt = inv(Ru) * p;
    y_opt(:,1) = w_opt' * un;
    
    noise_out = y_opt - x;
    
    gain = 10*log10( mean(mean(abs(rx_n - rx).^2)) / mean(abs(noise_out).^2));
    snr_opt(i) = gain + snr(i);
    
    [ ~,ber_opt(i) ] = biterr(x_bit,de2bi(qamdemod(y_opt,M)) );
    
    %% MMSE
    mmse = (dn * dn')/N_sample + w_opt' * Ru * w_opt;
    
    %% LMS
    mu = 1/trace(Ru);
    wn = rand(prod(N_tx_el),1) + 1i*rand(prod(N_tx_el),1);
    wn = wn ./ sum(abs(wn));
    for iter = 1:1000
        en = un' * wn - dn'; % error of actual iteration
        mse_grad = 2*un * en /N_sample; % gradient of the MSE with the actual wn
        wn = wn - 0.5 * mu * mse_grad;
    end
    y_lms(:,1) = wn' * un;
    
    noise_out = y_lms - x;
    
    gain = 10*log10( mean(mean(abs(rx_n - rx).^2)) / mean(abs(noise_out).^2));
    snr_lms(i) = gain + snr(i);
    
    
    [ ~,ber_lms(i) ] = biterr(x_bit,de2bi(qamdemod(y_lms,M)) );
    
    %% NULL
    y_null =  rx_n * transpose(w_h);
    noise_out = y_null - x;
    
    y_s = qamdemod(y_null,M);
    y_b = de2bi(y_s);
    [ ~,ber_null(i) ] = biterr(x_bit,y_b);
    
    
    gain = 10*log10(mean(mean(abs(rx_n - rx).^2)) / mean(abs(noise_out).^2));
    snr_null(i) = gain + snr(i);

    %% CONVENTIONAL
    y_conv = rx_n * transpose(w_h_conv);
    
    n_out_conv = y_conv - x;
    
    gain_conv = 10*log10( mean(mean(abs(rx_n - rx).^2)) / mean(abs(n_out_conv).^2) );
    snr_conv(i) = gain_conv + snr(i);
    
    [ ~,ber_conv(i) ] = biterr(x_bit,de2bi(qamdemod(y_conv,M)));
end

figure
plot(snr,snr_null,'gs-')
hold on
plot(snr,snr_conv,'gx--')
plot(snr,snr_mvdr,'bo-')
plot(snr,snr_opt,'cd-')
plot(snr,snr_lms,'rx-')

title("SNR Performance")

legend("Null","Conventional","MVDR","Wiener","LMS");
xlabel("SINR input");
ylabel("SINR output");
hold off;

ber_opt(ber_opt == 0) = 1e-5;
ber_lms(ber_lms == 0) = 1e-5;
ber_null(ber_null == 0) = 1e-5;
ber_mvdr(ber_mvdr == 0) = 2e-5;

figure
semilogy(snr,ber_null,'gs-')
hold on
semilogy(snr,ber_conv,'gx--')
semilogy(snr,ber_mvdr,'bo-')
semilogy(snr,ber_opt,'cd-')
semilogy(snr,ber_lms,'rx-')

title("BER")

legend("Null","Conventional","MVDR","Wiener","LMS");
xlabel("SINR input");
ylabel("BER");

real_angles



% figure
% semilogy(snr,ber,snr,ber_conv)
% legend("Null","Conventional");