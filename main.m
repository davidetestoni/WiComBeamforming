%% Initialize system constants
clear all,close all,clc;
rng(2014);
%given constants
gc = set2DParameters();

% Tunable parameters
tp.tx_pow = 20; %dbm

tp.tgt(1).range = 500;    % m
tp.tgt(1).angle = 7;       % degrees
tp.tgt(2).range = 500;    % m
tp.tgt(2).angle = 10;       % degrees

tp.intf(1).range = 1000;    % m
tp.intf(1).angle =   26;    % degrees
tp.intf(2).range = 1000;    % m
tp.intf(2).angle =   -50;    % degrees

tp.num_tx_el = 16;       
tp.steer_angle = 0;     % degrees

tp.N_bits = 1e4;
gc.N_blocks = ceil(tp.N_bits/gc.N_bit_block);


%% antenna definition

ula = phased.ULA( tp.num_tx_el, ...
        'ElementSpacing', 0.5*gc.lambda, ...
        'Element', phased.IsotropicAntennaElement('BackBaffled', true));

arrayresp = phased.ArrayResponse('SensorArray',ula,'WeightsInputPort',true);

%visualize2D(gc, tp,ula);

%% Signal generation
[data_bit,data_sample,data_qam,cp_sig] = modulateData(gc,tp);
% Normalization
norm_factor = mean(abs(cp_sig));
cp_sig = cp_sig ./ norm_factor;

% Amplification
cp_sig = cp_sig .* 10^(0.5 * tp.tx_pow/10);

%% Signal Collection
% target and interferer send signaling data to allow discovery

tp.tgt(1).tx_sig = cp_sig(:,:,1);
tp.tgt(1).data_bit = data_bit(:,:,1);


tp.tgt(2).tx_sig = cp_sig(:,:,11);
tp.tgt(2).data_bit = data_bit(:,:,11);

tp.intf(1).tx_sig = cp_sig(:,:,22);
tp.intf(1).data_bit = data_bit(:,:,22);

tp.intf(2).tx_sig = cp_sig(:,:,33);
tp.intf(2).data_bit = data_bit(:,:,33);

%free space loss
tp.tgt(1).rx_sig = tp.tgt(1).tx_sig.* 10^(-0.5*(fspl(tp.tgt(1).range,gc.lambda)/10));
tp.tgt(2).rx_sig = tp.tgt(2).tx_sig.* 10^(-0.5*(fspl(tp.tgt(2).range,gc.lambda)/10));
tp.intf(1).rx_sig = tp.intf(1).tx_sig.* 10^(-0.5*(fspl(tp.intf(1).range,gc.lambda)/10));
tp.intf(2).rx_sig = tp.intf(2).tx_sig.* 10^(-0.5*(fspl(tp.intf(2).range,gc.lambda)/10));

all_sig = [tp.tgt(1).rx_sig tp.tgt(2).rx_sig tp.intf(1).rx_sig tp.intf(2).rx_sig];
rx = collectPlaneWave(ula,all_sig,[tp.tgt(1).angle tp.tgt(2).angle  tp.intf(1).angle tp.intf(2).angle; 0 0 0 0],gc.fc);
%rx = collectPlaneWave(ula,[tp.tgt(1).rx_sig tp.tgt(2).rx_sig],[tp.tgt(1).angle tp.tgt(2).angle;0 0],gc.fc);

noise = sqrt(physconst('Boltzmann') * gc.nT * gc.B / 2) * randn(size(rx));
rx = rx + noise;
total_noise = mean(mean(abs(noise).^2))+ mean(abs(tp.tgt(2).rx_sig).^2) + mean(abs(tp.intf(1).rx_sig).^2)+ mean(abs(tp.intf(2).rx_sig).^2);
snr = 10*log10(mean(abs(tp.tgt(1).rx_sig).^2)/ total_noise);
%rx = awgn(rx,tp.snr,'measured');

%% DOA estimation
musicEstimator = phased.MUSICEstimator('SensorArray',ula,...
        'OperatingFrequency',gc.fc,'ScanAngles',gc.scan_az,...
        'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',4);

    
[~,doas] = musicEstimator(rx);

%figure
%plotSpectrum(musicEstimator,'NormalizeResponse',true);
disp(doas);
disp([tp.tgt(1).angle tp.tgt(2).angle  tp.intf(1).angle tp.intf(2).angle]);

%% Steering
num_tgt = length(tp.tgt);
steer_vec = steer_vec_ula(ula,gc.lambda,doas(1:num_tgt));
% weighted response
y = zeros(size(rx,1),num_tgt);
for j = 1 :size(steer_vec,2)
    s0 = steer_vec(:,j);
    w = s0 / tp.num_tx_el;
    for i = 1:size(rx,1)
        y(i,j) = w' * transpose(rx(i,:));
    end
end

%% amplify data
ampl_gain(1) = 10*log10(mean(abs(tp.tgt(1).tx_sig).^2) / mean(abs(y(:,1).^2)));
ampl_gain(2) = 10*log10(mean(abs(tp.tgt(2).tx_sig).^2)/ mean(abs(y(:,2).^2)));

y_amp(:,1) = y(:,1) .* 10^( 0.5 * ampl_gain(1)/10);
y_amp(:,2) = y(:,2) .* 10^( 0.5 * ampl_gain(2)/10);
% re-normalize
y_norm(:,1) = y_amp(:,1) ./ mean(abs(y_amp(:,1)));
y_norm(:,2) = y_amp(:,2) ./ mean(abs(y_amp(:,2)));
y_amp = y_norm .* norm_factor;



%% demodulate ofdm symbol

y_cap = zeros(size(data_qam,1),num_tgt);
y_sample = zeros(size(data_qam,1),num_tgt);
for i = 1:size(y_cap,2)
    y_cap(:,i) = ofdm_demod(y_amp(:,i),gc.N_fft,gc.N_cyclepref,gc.pilot_pos,gc.guard_bands);
    y_sample(:,i) = qamdemod(y_cap(:,i),gc.M);
end


%% BER
[~,ber(1)] = biterr(de2bi(y_sample(:,1)),tp.tgt(1).data_bit);
[~,ber(2)] = biterr(de2bi(y_sample(:,2)),tp.tgt(2).data_bit);
[~,ber_c(1)] = biterr(de2bi(y_sample(:,2)),tp.tgt(1).data_bit);
[~,ber_c(2)] = biterr(de2bi(y_sample(:,1)),tp.tgt(2).data_bit);
if ber(1) < ber_c(1)
    clear ber_c
else
    ber = ber_c;
    clear ber_c
end
