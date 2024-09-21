function transmitted_signal = pulseShaper(modulated_signal,oversampling_factor,pulse,wordLength,fractionLength)
% This function performs pulse shaping on a given modulated signal, with a
% given pulse.
% Oversampling_factor : Oversampling factor of the pulse
% Pulse               : Coeffiecients of the pulse shaping filter

filt_inp = upsample(modulated_signal,oversampling_factor);
transmitted_signal = fi(conv(filt_inp,pulse,'same'), 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);

end