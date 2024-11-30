function output_sig = polyphaseRRCfilter_quantized(framed_symbols, upsampling_factor, RRC_FIR_coefficients, wordLength, fractionLength)
%Polyphase implementation of the Pulse Shaping filter
%  'input_vars' = ...
%  'output_vars' = ...
%
%References: 
%[1] Digital Communications, John G Proakis and Masoud Salehi

output_len = length(framed_symbols)*upsampling_factor + length(RRC_FIR_coefficients) - 1;
output_sig = fi(zeros(output_len,1), 1, wordLength, fractionLength);

result = zeros(length(framed_symbols) + ceil(length(RRC_FIR_coefficients) / upsampling_factor) - 1, upsampling_factor);

result = fi(result, 1, wordLength, fractionLength);

for i = 1:upsampling_factor
    polyphase_filt = RRC_FIR_coefficients(i:upsampling_factor:end);
    temp = fi(conv(polyphase_filt, framed_symbols), 1, 12, 9);
    result(1:length(temp), i) = temp;
end

for i = 1:output_len
    p = mod((i-1),upsampling_factor) + 1;
    k = ceil(i/8);
    output_sig(i) = fi(result(k,p), 1, wordLength, fractionLength);
end

end
