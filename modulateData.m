function [data_bit,data_sample,data_qam,cp_sig] = modulateData(gc,tp)
%MODAULATEDATA Summary of this function goes here
%   Detailed explanation goes here
data_bit = randi([0 1],tp.N_bits,1);
data_sample = zeros(gc.N_bit_block/log2(gc.M),1,gc.N_blocks);
% zero padding
if mod(tp.N_bits,gc.N_bit_block) > 0
    data_bit(end+1 :  gc.N_blocks * gc.N_bit_block, 1) = 0;
end
% prepare data for qam
data_bit = reshape(data_bit,gc.N_bit_block/log2(gc.M),log2(gc.M),gc.N_blocks);

for i = 1 : size(data_bit,3)
    data_sample(:,:,i) = bi2de(data_bit(:,:,i),'left-msb');
end

data_qam = qammod(data_sample,gc.M);
% OFDM
cp_sig = ofdm_mod(data_qam,gc.pilot_pos,gc.pil,gc.guard_bands,gc.N_cyclepref,gc.N_fft,gc.N_blocks);

end

