function [rx_demod] = ofdm_demod(rx_sig,N_fft,N_cyclepref,pilot_pos,guard_bands)
%OFDM_DEMOD Summary of this function goes here
%   Detailed explanation goes here
rx_sig = rx_sig(N_cyclepref+1:end,:,:);
% move to freq
rx_data = fft(rx_sig,N_fft);
% remove guard bands and pilots
rx_data(pilot_pos,:,:) = [];
rx_demod= rx_data(guard_bands(1)+1 : end - guard_bands(2),:,:);
end

