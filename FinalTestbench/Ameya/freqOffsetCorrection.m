function [corrected_sig, detected_freq] = freqOffsetCorrection(rx_signal, sampling_freq)
%This function removes the frequency offset from the received signal

L = length(rx_signal);
t = 0:1/sampling_freq:(L-1)/sampling_freq;
F = (-L/2 : L/2 - 1) * sampling_freq/(L);

checker = rx_signal.^4;
dopp_shift = fft(double(checker));
dopp_shift = fftshift(dopp_shift);
mag = 10*log10(abs(dopp_shift));
plot(F,mag);
[~,I] = max(mag);

detected_freq = F(I)/4;
shift_back = exp(1i*2*pi*detected_freq*t);
corrected_sig = rx_signal' .* shift_back;


end