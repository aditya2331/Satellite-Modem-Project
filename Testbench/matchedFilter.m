function received_message = matchedFilter(received_signal,oversampling_factor,filter,wordLength,fractionLength)
% This function performs matched filtering on the received signal, and
% gives the output as received symbols
% filter : filter coefficients for the matched filter

filt_out = fi(conv(received_signal,filter,'same'), 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
received_message = downsample(filt_out,oversampling_factor);

end