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
for i = 1 : size(data_qam,3)
    data_ofdm(:,:,i) = cat(1,zeros(guard_bands(1),1),data_qam(:,:,i),zeros(guard_bands(2) + N_pil,1));
    for j = 1: max(size(pilot_pos))
        data_ofdm(:,:,i) = insert(data_ofdm(:,:,i),pil,pilot_pos(j));
    end
end
clear t i j;

%% IFFT
ifft_sig = ifft(data_ofdm,N_fft);

% add cycle prefix
cp_sig = zeros(N_fft + N_cyclepref,1,N_blocks);
cp_sig(1:N_cyclepref,:,:) = ifft_sig(end-N_cyclepref+1:end,:,:);
cp_sig(N_cyclepref + 1:end,:,:) = ifft_sig(:,:,:);

%stem(real(data_ofdm(:,:,5)));
%scatterplot(data_ofdm(:,:,1));


%% Freq Response
over_samp = 1e2;
df = delta_f/over_samp;
dt = 1/df;
t = 1:dt:Ts;
f = -20*delta_f:df:20*delta_f;
% we convolve int the freq with sinc 
s = sinc(f / delta_f);

% make ofdm continuos signal
ofdm_cont = zeros(N_fft * over_samp,1);
j = 1;
step = size(ofdm_cont,1)/N_fft;
for i = 1:size(data_ofdm(:,:,1),1)
    ofdm_cont(j) = data_ofdm(i,1,1);
    j = j + step;
end
clear j;
%stem(ofdm_cont)
%plot(f,s)

ofdm_freq = conv(s,ofdm_cont);
x_siz = size(ofdm_freq,1);
ff = floor(-x_siz/2) * df : df: (floor(x_siz/2)-1)*df;
plot(ff,log(abs(ofdm_freq).^2)*10);
xlabel('freq')
ylabel('PSD (dB)')
ylim([-100 100]);
%stem(data_ofdm(:,:,1));
%plot(real(ifft_sig(:,:,1)));

%stem(real(data_ofdm(:,:,1)));






%% 
a = [ 1 1 1 1 1 0 0 0 0];
a = insert(a,2,3)
a = insert(a,2,5)
a = insert(a,2,6)
a = insert(a,2,7)








function ins = insert(x,a,n)
    x(n+1:end) = x(n:end-1);
    x(n) = a;
    ins = x;
end