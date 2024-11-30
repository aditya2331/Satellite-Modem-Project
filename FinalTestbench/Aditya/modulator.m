function modulated_symbols = modulator(CC_encoded_bits, bits_per_symbol, constellation)
    
    if iscolumn(CC_encoded_bits)
        % Transpose to make it a row vector
        CC_encoded_bits = CC_encoded_bits.';
    end

    % Reshape the bit sequence into pairs
    symbols = reshape(CC_encoded_bits, bits_per_symbol, [])';
    % Map bit pairs to QPSK symbols
    modulated_symbols = constellation(bi2de(symbols, 'left-msb') + 1);

    % Ensure the output is a column vector
    modulated_symbols = modulated_symbols.';
end