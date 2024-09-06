function out = awgn_channel(transmitted_sig,EbN0,signal_power,oversampling_factor,wordLength,fractionLength)
%  This function simulates an AWGN channel by adding noise 
%  to the transmitted signal, EbN0 is required Eb/N0

beta = oversampling_factor*signal_power/EbN0;    % calculating noise variance
noise = fi(sqrt(beta/2)*(randn(size(transmitted_sig)) + 1i*randn(size(transmitted_sig))), 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
out = transmitted_sig + noise;

end