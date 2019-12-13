function cp_sig = ofdm_mod(data_qam,pilot_pos,pil,guard_bands,N_cyclepref,N_fft,N_blocks)
%OFDM_MOD Summary of this function goes here
%   Detailed explanation goes here

data_ofdm =  zeros(N_fft,1,size(data_qam,3));
% insert the pilots equi-spaced and padding
N_pil = max(size(pilot_pos));
for i = 1 : size(data_qam,3)
    data_ofdm(:,:,i) = cat(1,zeros(guard_bands(1),1),data_qam(:,:,i),zeros(guard_bands(2) + N_pil,1));
    for j = 1: max(size(pilot_pos))
        data_ofdm(:,:,i) = insert(data_ofdm(:,:,i),pil,pilot_pos(j));
    end
end

ifft_sig = ifft(data_ofdm,N_fft);

% add cycle prefix
cp_sig = zeros(N_fft + N_cyclepref,1,N_blocks);
cp_sig(1:N_cyclepref,:,:) = ifft_sig(end-N_cyclepref+1:end,:,:);
cp_sig(N_cyclepref + 1:end,:,:) = ifft_sig(:,:,:);

end

function ins = insert(x,a,n)
    x(n+1:end) = x(n:end-1);
    x(n) = a;
    ins = x;
end

