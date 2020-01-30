%% Initialize system constants
clear all,close all,clc;
%rng(2014);


% Tunable parameters
tp.tx_pow = 20; %dbm
tp.demo_area = [-1000 1000; -1000 1000; -20 0]; %the limits of demo area, in x;y;z 
tp.num_tx_el = 16;       
%given constants
tp.N_bits = 1e4;
gc = set2DParameters(tp);
gc.N_mob_stat = 4;
gc.N_tgt = 2;
gc.N_intf = 2;
%% MS POS

for i = 1:gc.N_tgt
    x = randi(tp.demo_area(1,:));
    y = randi(tp.demo_area(2,:));
    z = randi(tp.demo_area(3,:));
    tp.tgt(i).pos = [x;y;z];
    [az ,el, tp.tgt(i).range] = cart2sph(x,y,z);
    tp.tgt(i).azimuth = rad2deg(az);
    tp.tgt(i).elev = rad2deg(el);
end
for i = 1:gc.N_intf
    x = randi(tp.demo_area(1,:));
    y = randi(tp.demo_area(2,:));
    z = randi(tp.demo_area(3,:));
    tp.intf(i).pos = [x;y;z];
    [az ,el, tp.intf(i).range] = cart2sph(x,y,z);
    tp.intf(i).azimuth = rad2deg(az);
    tp.intf(i).elev = rad2deg(el);
end
tp.real_angles = [tp.tgt(1).azimuth tp.tgt(2).azimuth tp.intf(1).azimuth tp.intf(2).azimuth];
tp.real_angles = [tp.real_angles ;tp.tgt(1).elev tp.tgt(2).elev tp.intf(1).elev tp.intf(2).elev];
disp(tp.real_angles);



%% antenna definition

ula = phased.ULA( tp.num_tx_el, ...
        'ElementSpacing', 0.5*gc.lambda, ...
        'Element', phased.IsotropicAntennaElement('BackBaffled', false));

arrayresp = phased.ArrayResponse('SensorArray',ula,'WeightsInputPort',true);

%visualize2D(gc, tp,ula);

%% Signal generation
% modulate data for target
n_factors = zeros(gc.N_mob_stat,1);
for a = 1:length(tp.tgt)
    [data_bit,data_sample,data_qam,cp_sig] = modulateData(gc,tp);
    % Normalization
    tp.tgt(a).orig_sig_vec = cp_sig;
    norm_factor = mean(mean(abs(cp_sig)));
    cp_sig = cp_sig ./ norm_factor;
    n_factors(a) = norm_factor;
    % Amplification
    cp_sig = cp_sig .* 10^(0.5 * tp.tx_pow/10);
    tp.tgt(a).tx_sig_vec = cp_sig;
    tp.tgt(a).data_bit_vec = data_bit;
end
% modulate data for interferer
for a = 1:length(tp.intf)
    [data_bit,data_sample,data_qam,cp_sig] = modulateData(gc,tp);
    % Normalization
    tp.intf(a).orig_sig_vec = cp_sig;
    norm_factor = mean(mean(abs(cp_sig)));
    cp_sig = cp_sig ./ norm_factor;
    n_factors(a+2) = norm_factor;
    % Amplification
    cp_sig = cp_sig .* 10^(0.5 * tp.tx_pow/10);
    tp.intf(a).tx_sig_vec = cp_sig;
    tp.intf(a).data_bit_vec = data_bit;
end
norm_factor = mean(n_factors);
clear a;
%% Signal Collection
% target and interferer send signaling data to allow discovery
%free space loss
for ii = 1:size(tp.tgt(1).tx_sig_vec,3)

tp.tgt(1).tx_sig = tp.tgt(1).tx_sig_vec(:,:,ii);
tp.tgt(2).tx_sig = tp.tgt(2).tx_sig_vec(:,:,ii);
tp.intf(1).tx_sig = tp.intf(1).tx_sig_vec(:,:,ii);
tp.intf(2).tx_sig = tp.intf(2).tx_sig_vec(:,:,ii);

tp.tgt(1).rx_sig = tp.tgt(1).tx_sig.* 10^(-0.5*(fspl(tp.tgt(1).range,gc.lambda)/10));
tp.tgt(2).rx_sig = tp.tgt(2).tx_sig.* 10^(-0.5*(fspl(tp.tgt(2).range,gc.lambda)/10));
tp.intf(1).rx_sig = tp.intf(1).tx_sig.* 10^(-0.5*(fspl(tp.intf(1).range,gc.lambda)/10));
tp.intf(2).rx_sig = tp.intf(2).tx_sig.* 10^(-0.5*(fspl(tp.intf(2).range,gc.lambda)/10));

all_sig = [tp.tgt(1).rx_sig tp.tgt(2).rx_sig tp.intf(1).rx_sig tp.intf(2).rx_sig];
rx = collectPlaneWave(ula,all_sig,[tp.real_angles],gc.fc);
%rx = collectPlaneWave(ula,[tp.tgt(1).rx_sig tp.tgt(2).rx_sig],[tp.tgt(1).angle tp.tgt(2).angle;0 0],gc.fc);

%% AWGN
noise = sqrt(physconst('Boltzmann') * gc.nT * gc.B  *1e3) * randn(size(rx)); %noise in dBm

snr = 10*log10(mean(abs(tp.tgt(1).rx_sig).^2)/ mean(mean(abs(noise).^2)));
total_noise = mean(mean(abs(noise).^2))+ mean(abs(tp.tgt(2).rx_sig).^2) + mean(abs(tp.intf(1).rx_sig).^2)+ mean(abs(tp.intf(2).rx_sig).^2);
sinr = 10*log10(mean(abs(tp.tgt(1).rx_sig).^2)/ total_noise);

rx = rx + noise;

%% DOA estimation
musicEstimator = phased.MUSICEstimator('SensorArray',ula,...
        'OperatingFrequency',gc.fc,'ScanAngles',gc.scan_az,...
        'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',gc.N_mob_stat);

    
[~,doas] = musicEstimator(rx);

%figure
%plotSpectrum(musicEstimator,'NormalizeResponse',true);
%disp(doas);
%disp([tp.tgt(1).angle tp.tgt(2).angle  tp.intf(1).angle tp.intf(2).angle]);

%% Steering
num_tgt = length(tp.tgt);
steer_vec = steer_vec_ula(ula,gc.lambda,doas);
% weighted response
y = zeros(size(rx,1),length(doas));
for j = 1 :size(steer_vec,2)
    s0 = steer_vec(:,j);
    w = s0 / tp.num_tx_el;
    for i = 1:size(rx,1)
        y(i,j) = w' * transpose(rx(i,:));
    end
end

%% amplify data

% re-normalize
y_norm = y;
for aa = 1:length(doas)
    y_norm(:,aa) = y(:,aa) ./ mean(abs(y(:,aa)));
end
y_amp = y_norm .* norm_factor;
clear aa;


%% demodulate ofdm symbol

y_cap = zeros(size(data_qam,1),length(doas));
y_sample = zeros(size(data_qam,1),length(doas));
for i = 1:size(y_cap,2)
    y_cap(:,i) = ofdm_demod(y_amp(:,i),gc.N_fft,gc.N_cyclepref,gc.pilot_pos,gc.guard_bands);
    y_sample(:,i) = qamdemod(y_cap(:,i),gc.M);
end
%% find targets

tp.tgt(1).data_bit = tp.tgt(1).data_bit_vec(:,:,ii);
tp.tgt(2).data_bit = tp.tgt(2).data_bit_vec(:,:,ii);

tp.intf(1).data_bit = tp.intf(1).data_bit_vec(:,:,ii);
tp.intf(2).data_bit = tp.intf(2).data_bit_vec(:,:,ii);

tgt_err = zeros(gc.N_mob_stat,length(tp.tgt));
tgt_index = zeros(length(tp.tgt),1);
% check each signal with the right target
for j = 1:length(tp.tgt)
    for i = 1:gc.N_mob_stat
        [~,tgt_err(i,j)] = biterr(de2bi(y_sample(:,i),sqrt(gc.M)),tp.tgt(j).data_bit);
    end
    [~,tgt_index(j)] = min(tgt_err(:,j));
end


%% BER
for j = 1:length(tp.tgt)
    [~,ber(j,ii)] = biterr(de2bi(y_sample(:,tgt_index(j)),sqrt(gc.M)),tp.tgt(j).data_bit);
end
end
clear i j
