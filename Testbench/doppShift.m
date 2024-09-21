function  [shifted_signal,freq_offset] = doppShift(transmitted_sig, sampling_rate, oversampling_factor, N, M, wordLength, fractionLength)
%doppSHift This function shifts the input signal by a random frequency. It
%          also provides a random offset in time.

freq_offset = (sampling_rate/(M*oversampling_factor))*(1 - 2*rand);
t = 0:1/sampling_rate:(length(transmitted_sig)-1)/sampling_rate;
shift = fi(exp(-1i*2*pi*freq_offset*(t))', 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
shifted_signal = transmitted_sig .* shift;

end