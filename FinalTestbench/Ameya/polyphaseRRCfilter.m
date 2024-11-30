function output_sig = polyphaseRRCfilter(framed_symbols, upsampling_factor, RRC_FIR_coefficients)
%Polyphase implementation of the Pulse Shaping filter
%  'input_vars' = ...
%  'output_vars' = ...
%
%References: 
%[1] Digital Communications, John G Proakis and Masoud Salehi

output_len = length(framed_symbols)*upsampling_factor + length(RRC_FIR_coefficients) - 1;
output_sig = zeros(output_len,1);

result = zeros(length(framed_symbols) + ceil(length(RRC_FIR_coefficients) / upsampling_factor) - 1, upsampling_factor);

for i = 1:upsampling_factor
    polyphase_filt = RRC_FIR_coefficients(i:upsampling_factor:end);
    temp = conv(polyphase_filt, framed_symbols);
    result(1:length(temp), i) = temp;
end

for i = 1:output_len
    p = mod((i-1),upsampling_factor) + 1;
    k = ceil(i/upsampling_factor);
    output_sig(i) = result(k,p);
end

end

