function out = awgn_channel(transmitted_sig,EbN0,oversampling_factor,wordLength,fractionLength)
%  This function simulates an AWGN channel by adding noise to the transmitted signal
%  All the operations are in fixed point
%  EbN0 : Eb/N0 for the channel
%  tranmitted_sig: The transmitted signal, in fixed point
%  wordLength and fractionLength are fixed point parameters

% Calculating transmitted signal power---------------------------------------------
signal_power = mean(double(transmitted_sig).*conj(double(transmitted_sig)))/2;
% Calculating noise variance for required SNR--------------------------------------
beta = oversampling_factor*signal_power/EbN0;
% Generating noise ----------------------------------------------------------------
noise = fi(sqrt(beta/2)*(randn(size(transmitted_sig)) + 1i*randn(size(transmitted_sig))), 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
out = transmitted_sig + noise;

end