close all
clear all
clc
%%
N_bits = 1e4;
M = 16;
N_fft = 64;
pilot_pos = [12;26;40;54];
guard_bands = [6,5];
N_bit_block = (N_fft - max(size(pilot_pos)) - sum(guard_bands)) * log2(M); %block is 1 ofdm symbol
N_blocks = ceil(N_bits/N_bit_block);
N_cyclepref = 16;
N_carried_data = N_fft - numel(pilot_pos) - sum(guard_bands);

fc = 2.4e9;
lambda = physconst('LightSpeed')/fc;
delta_f = 15; %kHz
B = delta_f * N_fft;
Ts = 1/B;
%% Random data and Modulation
data_sample = zeros(N_bit_block/log2(M),1,N_blocks);
data_bit = randi([0 1],N_bits,1);
% zero padding
if mod(N_bits,N_bit_block) > 0
    data_bit(end+1 :  N_blocks * N_bit_block, 1) = 0;
end
% prepare data for qam
data_bit = reshape(data_bit,N_bit_block/log2(M),log2(M),N_blocks);

for i = 1 : size(data_bit,3)
    data_sample(:,:,i) = bi2de(data_bit(:,:,i),'left-msb');
end
clear i;

data_qam = qammod(data_sample,M);

%% OFDM
data_ofdm =  zeros(N_fft,1,size(data_qam,3));
% insert the pilots equi-spaced and padding
pil = 3 + 3j;
N_pil = max(size(pilot_pos));
N_q_data = size(data_qam,1); 
cp_sig = ofdm_mod(data_qam,pilot_pos,pil,guard_bands,N_cyclepref,N_fft,N_blocks);
%
t = 0:1/B *1e3:((1/delta_f)-1/B) *1e3;
%bar(t,abs(ifft_sig(:,:,5)));
%xlabel('time(us)')
%ylabel('Amplitude')




%% RX
% channel
rx_sig = awgn(cp_sig,60,'measured');
% remove cp
rx_sig = rx_sig(N_cyclepref+1:end,:,:);
% move to freq
rx_data = fft(rx_sig,N_fft);
% remove guard bands and pilots
rx_data(pilot_pos,:,:) = [];
rx_data = rx_data(guard_bands(1)+1 : end - guard_bands(2),:,:);

demod_sample = qamdemod(rx_data,M);
demod_bit = zeros(size(data_bit));
numErr = zeros(size(demod_sample,3),1);
ber = zeros(size(demod_sample,3),1);
for i = 1: size(demod_sample,3)
    demod_bit(:,:,i) = de2bi(demod_sample(:,:,i),'left-msb');
    [numErr(i),ber(i)]= biterr(data_bit(:,:,i),demod_bit(:,:,i));    
end
totErr = sum(numErr);
compBER =totErr /  N_bits
avgBER = mean(ber)




%% reference modulator
modulatorOFDM = comm.OFDMModulator(...
    'FFTLength' ,                      N_fft,...
    'NumGuardBandCarriers', guard_bands',...
    'InsertDCNull',                   false, ...
    'PilotInputPort',                 true,...
    'PilotCarrierIndices',         pilot_pos,...
    'CyclicPrefixLength',         N_cyclepref);
    %'NumSymbols',                  numDataSymbols,...
    %'NumTransmitAntennas',  numTx);
    
    pils = repmat(pil,4,1);
    ref_ofdm = modulatorOFDM(data_qam(:,:,1),pils);
    
    demodulatorOFDM = comm.OFDMDemodulator(...
     'FFTLength' ,                      N_fft,...
     'NumGuardBandCarriers', guard_bands',...
     'RemoveDCCarrier',           false, ...
     'PilotOutputPort',           true, ...
     'PilotCarrierIndices',        pilot_pos,...
     'CyclicPrefixLength',         N_cyclepref);
     %'NumSymbols',                  numDataSymbols);
     
     ref_rx = demodulatorOFDM(ref_ofdm);
    

%% Collect
N_ant = 4;
ula =phased.ULA( N_ant, ...
        'ElementSpacing', 0.5*lambda, ...
        'Element', phased.IsotropicAntennaElement('BackBaffled', true));

y = zeros(size(cp_sig,1),N_ant,size(cp_sig,3));
for i = 1:size(cp_sig,3)
    y(:,:,i) = collectPlaneWave(ula,cp_sig(:,:,i),[2; 0],fc,physconst('LightSpeed'));
end

y_freq = zeros(N_carried_data,N_ant,size(y,3));
y_sample = zeros(N_carried_data,N_ant,size(y,3));
for i = 1 : N_ant
    y_freq (:,i,:)= ofdm_demod(y(:,i,:),N_fft,N_cyclepref,pilot_pos,guard_bands);
    y_sample(:,i,:)= qamdemod(y_freq(:,i,:),M);
end

y_ber = zeros(size(y_sample));
for i = 1:size(y_sample,2)
    y_ber(:,i,:) = (y_sample(:,i,:) == data_sample);
    
end



%% Freq Response
over_samp = 1e2;
df = delta_f/over_samp;
dt = 1/df;
t = 1:dt:Ts;
f = -20*delta_f:df:20*delta_f;
% we convolve int the freq with sinc 
s = sinc(f / delta_f);

% make ofdm continuos signal
ofdm_cont = zeros(N_fft * over_samp,1,size(data_ofdm,3));
ofdm_freq = zeros(size(ofdm_cont,1) + max(size(s)) -1,size(ofdm_cont,2),size(ofdm_cont,2)) ;
j = 1;
for k = 1 : size(data_ofdm,3)
    for i = 1:size(data_ofdm(:,:,k),1)
        ofdm_cont(j) = data_ofdm(i,1,k);
        j = j + over_samp;
    end
    ofdm_freq(:,1,k) = conv(s,ofdm_cont(:,1,k));
end
clear i j k;
%stem(ofdm_cont)
%plot(f,s)


x_siz = size(ofdm_freq,1);
ff = floor(-x_siz/2) * df : df: (floor(x_siz/2)-1)*df;
plot(ff,log(abs(ofdm_freq(:,:,1)).^2)*10);
xlabel('freq')
ylabel('PSD (dB)')
ylim([-100 100]);
xlim([-N_fft/2*delta_f N_fft/2*delta_f])
%stem(data_ofdm(:,:,1));
%plot(real(ifft_sig(:,:,1)));

%stem(real(data_ofdm(:,:,1)));


%% averge on all symbols

x_siz = size(ofdm_freq,1);
ff = floor(-x_siz/2) * df : df: (floor(x_siz/2)-1)*df;
plot(ff,log(abs(mean(ofdm_freq,3)).^2)*10);
xlabel('freq')
ylabel('PSD (dB)')
ylim([-100 100]);
xlim([-N_fft/2*delta_f N_fft/2*delta_f])









function ins = insert(x,a,n)
    x(n+1:end) = x(n:end-1);
    x(n) = a;
    ins = x;
end