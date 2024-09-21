function [received_message,detected_freq] = matchedFilter(received_signal,oversampling_factor,filter,M,symbol_rate,wordLength,fractionLength)
% This function performs matched filtering on the received signal, and
% gives the output as received symbols
% filter : filter coefficients for the matched filter

filt_out = fi(conv(received_signal,filter,'same'), 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
[corrected_sig,detected_freq] = FreqOffsetCorrector(filt_out,M,symbol_rate);
received_message = downsample(corrected_sig,oversampling_factor);
received_message = received_message';


end