function modulatedSig = qpskMod1(in_signal,wordLength,fractionLength)
% This function takes a sequence of bits as input and does QPSK modulation
% on the sequemce. Modulation is done in fixed point. The output sequence
% is half the length of input sequence.

modulatedSig = zeros(length(in_signal)/2,1);
% Constellation map for QPSK
symbol_map = [1+1i, -1+1i, 1-1i, -1-1i] / sqrt(2);
symbol_map = fi(symbol_map,1,wordLength,fractionLength,'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);

% Modulation
for i = 1:2:length(in_signal)
    if in_signal(i) == 0 && in_signal(i+1) == 0
        modulatedSig((i+1)/2) = symbol_map(1);
    elseif in_signal(i) == 0 && in_signal(i+1) == 1
        modulatedSig((i+1)/2) = symbol_map(2);
    elseif in_signal(i) == 1 && in_signal(i+1) == 0
        modulatedSig((i+1)/2) = symbol_map(3);
    else
        modulatedSig((i+1)/2) = symbol_map(4);
    end
end


end