function [out_sig , offset_freq] = FreqOffsetCorrector(in_sig,M,Fs)
%   This function corrects the frequency offset in a given signal
%   in_sig : Input signal with frequency offset
%   M  : No. of constellation points in MPSK modulation
%   Fs : Sampling frequency

L = length(in_sig);
t = 0:1/Fs:(L-1)/Fs;
F = (-L/2 : L/2 - 1) * Fs/(L);

checker = in_sig.^M;
dopp_shift = fft(double(checker));
dopp_shift = fftshift(dopp_shift);
mag = 10*log10(abs(dopp_shift));
[~,I] = max(mag);

offset_freq = F(I)/4;
shift_back = exp(1i*2*pi*offset_freq*t);
out_sig = in_sig' .* shift_back;

end