function shifted_sig = addFreqOffset(noisy_signal, freq_offset, sampling_rate)
% Adds Frequency offset to the input signal

t = 0:1/sampling_rate:(length(noisy_signal)-1)/sampling_rate;
shift = exp(-1i*2*pi*freq_offset*(t))';
shifted_sig = noisy_signal .* shift;

end