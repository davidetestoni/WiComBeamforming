close all
clear all
clc
%%
N_bits = 1e4;
N_bit_block = 96 * 2; %doublesize for FEC (hypotize already in the data)
N_fft = 64; % has to match with the block size
N_blocks = ceil(N_bits/N_bit_block);
N_cyclepref = 16;
M = 16;
pilot_pos = [12;26;40;54];
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

data_qam = qammod(data_sample,M);

%% OFDM
data_ofdm =  zeros(N_fft,1,size(data_qam,3));
% insert the pilots equi-spaced and padding
pil = 3 + 3j;
N_pil = max(size(pilot_pos));
N_q_data = size(data_qam,1); 
for i = 1 : size(data_qam,3)
    t = cat(1,data_qam(:,:,i),zeros(N_fft-N_q_data,1));
    for j = pilot_pos
        data_ofdm(:,:,i) = insert(t,pil,j);
    end
end
clear t i j;

%% IFFT
ifft_sig = ifft(data_ofdm,N_fft);

% add cycle prefix
cp_sig = zeros(N_fft + N_cyclepref,1,N_blocks);
cp_sig(1:N_cyclepref,:,:) = ifft_sig(end-N_cyclepref+1:end,:,:);
cp_sig(N_cyclepref + 1:end,:,:) = ifft_sig(:,:,:);

%stem(real(data_ofdm(:,:,1)));
%scatterplot(data_ofdm(:,:,1));


%% Freq Response
df = delta_f*1e-2;
dt = 1/df;
t = 1:dt:Ts;
f = -20*delta_f:df:20*delta_f;
s = sinc(f / delta_f);
ofdm_cont = zeros(N_fft * 1e2,1);
j = 1;
step = size(ofdm_cont,1)/N_fft;
for i = 1:size(data_ofdm(:,:,1),1)
    ofdm_cont(j) = data_ofdm(i,1,1);
    j = j + step;
end
%stem(ofdm_cont)
%plot(f,s)
grid

r1 = cat(1,ones(100,1),zeros(50,1));
r2 = ones(200,1);

xx = conv(r1,r2);
plot(xx);

%%

x = conv(s,ofdm_cont);
x_siz = size(x,1);
ff = floor(-x_siz/2) * df : df: (floor(x_siz/2)-1)*df;
semilogy(ff,abs(x).^2);
%stem(data_ofdm(:,:,1));
%plot(real(ifft_sig(:,:,1)));

%stem(real(data_ofdm(:,:,1)));

%% y = sinc(f);
for x=data_ofdm(:,1,1)
    x = abs(x);
    yy = yy + conv(y,x);
    plot(yy);
    pause
    clear x;
end














function ins = insert(x,a,n)
    x(n+1:end) = x(n:end-1);
    x(n) = a;
    ins = x;
end